"""
modules/scp_module.py

Streamlit module for Single-Cell Proteomics (SCP) analysis in Pro-Visualize.

Tab layout:
  1. Upload & QC       — load DIA-NN PG matrix + report.stats + annotation;
                         QC violin plots, scatter, protein detection histograms
  2. Preprocessing     — normalize → log1p → regress covariates → scale;
                         library-size plots, covariate-correlation heatmap
  3. Embedding         — PCA elbow, PCA scatter, UMAP, Leiden clustering,
                         cluster composition bar charts
  4. Diff. Expression  — group selection, Wilcoxon DE, volcano plot, heatmap,
                         protein class breakdown table
  5. Activity Scoring  — paste / define gene sets; UMAP + violin panels;
                         per-score summary stats
"""

import pickle
import streamlit as st
import pandas as pd
import logging
import plotly.express as px
from plotly import colors as pc

from visualizations.scp_visualizer import SCPVisualizer
from utils.plot_manager import PlotManager, MplPlotManager

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Session-state helpers
# ─────────────────────────────────────────────────────────────────────────────

def _viz() -> SCPVisualizer:
    return st.session_state.scp_visualizer


def _has_viz() -> bool:
    return (
        "scp_visualizer" in st.session_state
        and st.session_state.scp_visualizer is not None
    )


def _state(key: str, default=None):
    if key not in st.session_state:
        st.session_state[key] = default
    return st.session_state[key]


def _pm(plot_key: str) -> PlotManager:
    """Create a PlotManager pre-stamped with the SCP module tag."""
    pm = PlotManager(plot_key)
    pm.module = "scp"
    return pm


def _mpl_pm(plot_key: str) -> MplPlotManager:
    """Create an MplPlotManager pre-stamped with the SCP module tag."""
    pm = MplPlotManager(plot_key)
    pm.module = "scp"
    return pm


# ─────────────────────────────────────────────────────────────────────────────
# Main class
# ─────────────────────────────────────────────────────────────────────────────

class SCPTab:
    """Encapsulates the full UI and logic for the Single-Cell Proteomics tab."""

    PLOT_KEYS = [
        "scp_qc_overview", "scp_det_hist", "scp_qc_scatter",
        "scp_stats_metric", "scp_lib_size", "scp_cov_heatmap",
        "scp_pca", "scp_elbow", "scp_umap", "scp_cluster_comp",
        "scp_expr_umap", "scp_volcano", "scp_de_heatmap",
        "scp_activity_umap", "scp_activity_violin",
    ]

    def __init__(self):
        if "scp_visualizer" not in st.session_state:
            st.session_state.scp_visualizer = None
        if "scp_de_results" not in st.session_state:
            st.session_state.scp_de_results = {}
        if "scp_score_cols" not in st.session_state:
            st.session_state.scp_score_cols = []
        if "scp_enr_results" not in st.session_state:
            st.session_state.scp_enr_results = None

    # ─────────────────────────────────────────────────────────────────────────
    # GLOBAL SETTINGS
    # ─────────────────────────────────────────────────────────────────────────

    def _render_global_settings(self) -> dict:
        with st.expander("🎨 Global Plot Settings", expanded=False):
            c1, c2 = st.columns(2)
            theme_opts = {
                "Standard White": "plotly_white",
                "Dark Mode": "plotly_dark",
                "Minimal": "simple_white",
            }
            theme = theme_opts[
                c1.selectbox("Theme", list(theme_opts.keys()), key="scp_theme")
            ]
            color_map = {}
            if _has_viz():
                groups_cols = _viz().get_available_groupby_cols()
                if groups_cols:
                    groupby_for_colors = c2.selectbox(
                        "Preview colors for column:", groups_cols, key="scp_color_col"
                    )
                    groups = _viz().adata.obs[groupby_for_colors].unique().tolist()
                    try:
                        safe_hex = [
                            "#%02x%02x%02x" % tuple(int(x) for x in pc.unlabel_rgb(c))
                            for c in px.colors.qualitative.Safe
                        ]
                    except Exception:
                        safe_hex = px.colors.qualitative.Safe
                    cols = st.columns(min(4, len(groups)))
                    for i, g in enumerate(groups):
                        with cols[i % 4]:
                            color_map[g] = st.color_picker(
                                f"{g}", value=safe_hex[i % len(safe_hex)],
                                key=f"scp_color_{g}"
                            )
            if st.button(
                "Apply to All Plots", type="primary", use_container_width=True,
                key="scp_apply_settings"
            ):
                for k in self.PLOT_KEYS:
                    st.session_state[f"{k}_fig"] = None
                st.success("Settings applied — plots will regenerate.")
                st.rerun()
        kwargs = {"template": theme}
        if color_map:
            kwargs["color_discrete_map"] = color_map
        return kwargs

    # ─────────────────────────────────────────────────────────────────────────
    # TAB 1 — UPLOAD & QC
    # ─────────────────────────────────────────────────────────────────────────

    def _tab_upload_qc(self, global_kwargs: dict):
        # ── File uploaders ──────────────────────────────────────────────────
        with st.expander(
            "📂 Upload Files", expanded=(not _has_viz())
        ):
            st.markdown(
                """
                **Required files:**
                - **PG Matrix** — DIA-NN `*_pg_matrix.tsv` (wide format: proteins × runs)
                - **Annotation** — CSV/TSV with a sample-name column (matching PG matrix columns) plus grouping metadata
                
                **Optional:**
                - **report.stats** — DIA-NN per-run QC metrics (FWHM.RT, Normalisation.Instability, etc.)
                """
            )
            c1, c2, c3 = st.columns(3)
            pg_file = c1.file_uploader(
                "PG Matrix (.tsv/.csv)", type=["tsv", "txt", "csv"], key="scp_pg_up"
            )
            stats_file = c2.file_uploader(
                "report.stats (optional)", type=["tsv", "txt", "csv"], key="scp_stats_up"
            )
            anno_file = c3.file_uploader(
                "Annotation (.csv/.tsv)", type=["tsv", "txt", "csv"], key="scp_anno_up"
            )

            st.markdown("---")
            st.markdown("**Column Configuration**")
            ca, cb, cc = st.columns(3)
            sample_col = ca.text_input(
                "Sample column (in Annotation)", value="Run", key="scp_sample_col"
            )
            protein_col = cb.text_input(
                "Protein column (in PG Matrix)", value="Protein.Group", key="scp_prot_col"
            )
            run_col = cc.text_input(
                "Run column (in report.stats)", value="File.Name", key="scp_run_col"
            )

            if st.button(
                "Load Data", type="primary", use_container_width=True, key="scp_load_btn"
            ):
                if not (pg_file and anno_file):
                    st.warning("Please upload at least the PG Matrix and Annotation files.")
                else:
                    try:
                        pg_df = pd.read_csv(pg_file, sep=None, engine="python")
                        anno_df = pd.read_csv(anno_file, sep=None, engine="python")
                        stats_df = (
                            pd.read_csv(stats_file, sep=None, engine="python")
                            if stats_file
                            else pd.DataFrame()
                        )
                        with st.spinner("Building AnnData object…"):
                            st.session_state.scp_visualizer = SCPVisualizer(
                                pg_df, stats_df, anno_df,
                                sample_col=sample_col,
                                run_col=run_col,
                                protein_col=protein_col,
                            )
                            st.session_state.scp_visualizer.compute_qc_metrics()
                        st.success(
                            f"Loaded: {_viz().adata.n_obs} cells × {_viz().adata.n_vars} proteins"
                        )
                        st.rerun()
                    except Exception as e:
                        st.error(f"Error loading data: {e}")
                        logger.error(f"SCP load error: {e}", exc_info=True)

        if not _has_viz():
            st.info("Upload files above to begin analysis.")
            return

        viz = _viz()
        st.markdown("---")

        # ── Gene Name Annotation ──────────────────────────────────────────────
        def _n_valid_genes(series):
            return int((
                series.notna()
                & ~series.astype(str).str.strip().str.lower().isin({"", "nan", "none", "na"})
            ).sum())

        _GENE_COL_CANDIDATES = ["Genes", "Gene.Names", "Gene names", "Gene", "gene_names"]
        _detected_gene_col = next(
            (c for c in _GENE_COL_CANDIDATES
             if c in viz.adata.var.columns and _n_valid_genes(viz.adata.var[c]) > 0),
            None,
        )
        if _detected_gene_col:
            _n_sym = _n_valid_genes(viz.adata.var[_detected_gene_col])
            st.caption(
                f"Gene names: `{_detected_gene_col}` — "
                f"{_n_sym}/{viz.adata.n_vars} proteins have a gene symbol."
            )
        else:
            st.warning(
                "No gene name column found in the PG matrix (or the column is empty). "
                "Gene symbols are required for pathway enrichment. "
                "Click below to fetch them from UniProt using your protein accessions."
            )
            if st.button(
                "Fetch Gene Names from UniProt",
                key="scp_fetch_uniprot",
                use_container_width=True,
            ):
                with st.spinner(
                    f"Querying UniProt for {viz.adata.n_vars} proteins… "
                    "(batched, may take ~10–30 s)"
                ):
                    try:
                        n_mapped = viz.annotate_gene_names_from_uniprot()
                        if n_mapped == 0:
                            st.error(
                                "UniProt returned 0 gene symbols. "
                                "Check that your protein IDs are UniProt accessions "
                                "(e.g. P12345). If they are a different format, "
                                "consider using a gene-name column as 'Protein column' at upload."
                            )
                        else:
                            st.success(
                                f"Annotated {n_mapped}/{viz.adata.n_vars} proteins with "
                                f"gene symbols from UniProt."
                            )
                            st.rerun()
                    except Exception as e:
                        st.error(f"UniProt lookup failed: {e}")

        # ── QC Metrics Table ─────────────────────────────────────────────────
        with st.expander("📋 Data Overview", expanded=False):
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Cells", viz.adata.n_obs)
            c2.metric("Proteins", viz.adata.n_vars)
            c3.metric("Median n_proteins",
                      int(viz.adata.obs["n_proteins"].median()))
            c4.metric("Mean detection %",
                      f"{viz.adata.obs['pct_detected'].mean():.1f}%")

            _ov_cell_tab, _ov_prot_tab = st.tabs(["Cell Metadata", "Protein Metadata"])
            with _ov_cell_tab:
                st.dataframe(viz.get_obs_df().describe(include="all").T, height=200)
            with _ov_prot_tab:
                _var_df = viz.adata.var.copy()
                _var_df.index.name = "protein_id"
                _GENE_COLS = ["Genes", "Gene.Names", "Gene names", "Gene", "gene_names"]
                _gene_col = next((c for c in _GENE_COLS if c in _var_df.columns), None)
                if _gene_col:
                    n_with_gene = _n_valid_genes(_var_df[_gene_col])
                    st.caption(
                        f"Gene name column: `{_gene_col}` — "
                        f"{n_with_gene}/{len(_var_df)} proteins have a gene symbol."
                        + (" ⚠️ Column present but all empty." if n_with_gene == 0 else "")
                    )
                else:
                    st.caption("No gene name column in adata.var. Use 'Fetch Gene Names from UniProt' above.")
                st.dataframe(_var_df.reset_index(), use_container_width=True, height=300)

        # ── Groupby selector ─────────────────────────────────────────────────
        groupby_cols = viz.get_available_groupby_cols()
        color_by = st.selectbox(
            "Color plots by:", ["None"] + groupby_cols, key="scp_qc_color"
        )
        color_by = None if color_by == "None" else color_by

        qc_overview_tab, det_hist_tab, scatter_tab, stats_tab = st.tabs([
            "QC Violin Overview", "Protein Detection", "QC Scatter", "DIA-NN Stats Metrics"
        ])

        with qc_overview_tab:
            pm = _pm("scp_qc_overview")
            pm.render_generate_button(
                viz.plot_qc_overview, color_by=color_by, **global_kwargs
            )
            pm.render_plot_and_editor()

        with det_hist_tab:
            pm = _pm("scp_det_hist")
            pm.render_generate_button(
                viz.plot_protein_detection_histogram, **global_kwargs
            )
            pm.render_plot_and_editor()

            st.markdown("**Per-Protein QC Table**")
            st.dataframe(viz.get_var_df().sort_values("detection_rate").head(50))

        with scatter_tab:
            c1, c2 = st.columns(2)
            num_cols = viz.get_available_numeric_cols()
            x_col = c1.selectbox("X axis", num_cols,
                                  index=num_cols.index("n_proteins") if "n_proteins" in num_cols else 0,
                                  key="scp_scatter_x")
            y_col = c2.selectbox("Y axis", num_cols,
                                  index=num_cols.index("total_intensity") if "total_intensity" in num_cols else min(1, len(num_cols)-1),
                                  key="scp_scatter_y")
            pm = _pm("scp_qc_scatter")
            pm.render_generate_button(
                viz.plot_n_proteins_scatter,
                x_col=x_col, y_col=y_col, color_col=color_by,
                **global_kwargs
            )
            pm.render_plot_and_editor()

        with stats_tab:
            diann_cols = viz.get_available_diann_stats()
            if not diann_cols:
                st.info("No DIA-NN stats columns available. Upload a report.stats file.")
            else:
                metric = st.selectbox(
                    "DIA-NN Stats Metric", diann_cols, key="scp_stats_metric_sel"
                )
                pm = _pm("scp_stats_metric")
                pm.render_generate_button(
                    viz.plot_stats_metric, metric=metric, groupby=color_by,
                    **global_kwargs
                )
                pm.render_plot_and_editor()

        # ── Filtering controls ────────────────────────────────────────────────
        st.markdown("---")
        st.subheader("🔧 QC Filtering")
        with st.expander("Configure & Apply Filters", expanded=False):
            fc1, fc2 = st.columns(2)
            min_prots = fc1.number_input(
                "Min proteins per cell", min_value=1, value=100, step=10, key="scp_min_prot"
            )
            max_prots = fc2.number_input(
                "Max proteins per cell (0 = no limit)", min_value=0, value=0, step=50,
                key="scp_max_prot"
            )
            min_det_pct = st.slider(
                "Min protein detection % (across cells)", 0.0, 50.0, 10.0, 1.0,
                key="scp_min_det_pct",
                help="Proteins detected in fewer than this % of cells are removed."
            )

            fc3, fc4 = st.columns(2)
            extra_metric = fc3.selectbox(
                "Additional QC metric filter", ["None"] + viz.get_available_diann_stats(),
                key="scp_extra_metric_sel"
            )
            if extra_metric != "None":
                extra_thresh = fc4.number_input(
                    f"Max {extra_metric} threshold",
                    value=float(viz.adata.obs[extra_metric].median()),
                    key="scp_extra_thresh"
                )

            if st.button("Apply Filters", type="primary", use_container_width=True,
                         key="scp_apply_filters"):
                n_cells_before = viz.adata.n_obs
                n_prot_before = viz.adata.n_vars
                with st.spinner("Filtering…"):
                    viz.filter_samples(
                        min_proteins=min_prots,
                        max_proteins=max_prots if max_prots > 0 else None
                    )
                    viz.filter_proteins(min_detection_pct=min_det_pct)
                    if extra_metric != "None":
                        viz.filter_by_qc_metric(extra_metric, extra_thresh, "below")
                st.success(
                    f"After filtering: {viz.adata.n_obs} cells (from {n_cells_before}), "
                    f"{viz.adata.n_vars} proteins (from {n_prot_before})"
                )
                # Invalidate downstream plots
                for k in self.PLOT_KEYS:
                    st.session_state[f"{k}_fig"] = None
                st.rerun()

    # ─────────────────────────────────────────────────────────────────────────
    # TAB 2 — PREPROCESSING
    # ─────────────────────────────────────────────────────────────────────────

    def _tab_preprocessing(self, global_kwargs: dict):
        if not _has_viz():
            st.info("Complete the Upload & QC step first.")
            return

        viz = _viz()

        # -- Pipeline Status --------------------------------------------------
        state = viz.pp_state
        cols = st.columns(6)
        steps = [
            ("Normalized", "normalized"), ("Log1p", "log_transformed"),
            ("Regressed", "regressed"), ("Scaled", "scaled"),
            ("PCA", "pca_computed"), ("Batch Corrected", "batch_corrected"),
        ]
        for col, (label, key) in zip(cols, steps):
            icon = "✅" if state[key] else "⬜"
            col.metric(label, icon)

        st.markdown("---")

        # -- Preprocessing controls -------------------------------------------
        with st.expander("⚙️ Preprocessing Pipeline", expanded=not state["normalized"]):
            st.markdown("Steps are applied in order. Re-running resets from the raw input layer.")

            c1, c2 = st.columns(2)
            do_norm = c1.toggle("Normalize (median total count)", value=True, key="scp_do_norm")
            do_log = c2.toggle("Log1p transform", value=True, key="scp_do_log")

            st.markdown("**Regress out technical covariates**")
            num_cols = viz.get_available_numeric_cols()
            default_covs = [
                c for c in ["log_n_proteins", "Normalisation.Instability"] if c in num_cols
            ]
            regress_covs = st.multiselect(
                "Covariates to regress out:", num_cols,
                default=default_covs, key="scp_regress_covs",
                help="Strongly correlated with PC1/PC2 = technical noise. Use the correlation heatmap below to decide."
            )
            do_scale = st.toggle(
                "Standard scale (unit variance per protein)",
                value=False, key="scp_do_scale",
                help="Gives equal weight to all proteins in PCA/clustering. "
                     "Useful when low-abundance signalling proteins carry key biology."
            )

            if st.button(
                "Run Preprocessing", type="primary", use_container_width=True,
                key="scp_run_preprocess"
            ):
                with st.spinner("Running preprocessing pipeline…"):
                    viz.preprocess(
                        normalize=do_norm,
                        log_transform=do_log,
                        regress_covariates=regress_covs if regress_covs else None,
                        scale=do_scale,
                    )
                for k in self.PLOT_KEYS:
                    st.session_state[f"{k}_fig"] = None
                st.success("✓ Preprocessing complete")
                st.rerun()

        # -- PCA controls -----------------------------------------------------
        with st.expander("🔬 PCA & Batch Correction", expanded=not state["pca_computed"]):
            pc_c1, pc_c2 = st.columns(2)
            n_comps = pc_c1.slider(
                "Number of PCs", 10, 100, 50, 5, key="scp_n_comps"
            )
            if st.button("Run PCA", type="primary", use_container_width=True,
                         key="scp_run_pca"):
                with st.spinner("Running PCA…"):
                    viz.run_pca(n_comps=n_comps)
                for k in ["scp_pca", "scp_elbow", "scp_cov_heatmap"]:
                    st.session_state[f"{k}_fig"] = None
                st.success(f"✓ PCA complete ({n_comps} PCs)")
                st.rerun()

            st.markdown("---")
            st.markdown("**Harmony Batch Correction** (optional)")
            groupby_cols = viz.get_available_groupby_cols()
            batch_key = st.selectbox(
                "Batch key", ["None"] + groupby_cols, key="scp_batch_key"
            )
            if batch_key != "None":
                if st.button("Run Harmony", use_container_width=True, key="scp_run_harmony"):
                    try:
                        with st.spinner("Running Harmony…"):
                            viz.run_harmony(batch_key)
                        st.success(f"✓ Harmony applied on '{batch_key}'")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Harmony failed: {e}")

        # -- Preprocessing visualisations ------------------------------------
        if not state["normalized"]:
            st.info("Run preprocessing to see visualisations.")
            return

        # -- Session save / load -----------------------------------------------
        with st.expander("💾 Save / Load Session", expanded=False):
            st.caption(
                "Save the full analysis state (AnnData + all preprocessing steps) "
                "to a file. Load it in a future session to resume from where you left off."
            )
            sc1, sc2 = st.columns(2)
            with sc1:
                try:
                    session_bytes = pickle.dumps(viz)
                    sc1.download_button(
                        "⬇ Download session (.pkl)",
                        data=session_bytes,
                        file_name="scp_session.pkl",
                        mime="application/octet-stream",
                        key="scp_save_session",
                        use_container_width=True,
                    )
                except Exception as e:
                    sc1.error(f"Could not serialize session: {e}")
            with sc2:
                uploaded_session = sc2.file_uploader(
                    "Load session (.pkl)", type=["pkl"], key="scp_load_session_file",
                    label_visibility="collapsed",
                )
                if uploaded_session is not None:
                    if sc2.button("Restore session", use_container_width=True, key="scp_restore_btn"):
                        try:
                            restored = pickle.loads(uploaded_session.read())
                            if not isinstance(restored, type(viz)):
                                st.error("File does not contain a valid SCP session.")
                            else:
                                st.session_state.scp_visualizer = restored
                                for k in self.PLOT_KEYS:
                                    st.session_state.pop(f"{k}_fig", None)
                                st.success("Session restored. Reload the tab to reflect the loaded state.")
                                st.rerun()
                        except Exception as e:
                            st.error(f"Failed to restore session: {e}")

        lib_tab, corr_tab = st.tabs(["📏 Library Sizes", "🔗 Covariate Correlations"])

        with lib_tab:
            pm = _pm("scp_lib_size")
            pm.render_generate_button(viz.plot_library_size_comparison, **global_kwargs)
            pm.render_plot_and_editor()

        with corr_tab:
            if not state["pca_computed"]:
                st.info("Run PCA first to see covariate correlations.")
            else:
                pm = _pm("scp_cov_heatmap")
                pm.render_generate_button(
                    viz.plot_covariate_correlation_heatmap, **global_kwargs
                )
                pm.render_plot_and_editor()
                with st.expander("Raw correlation table"):
                    st.dataframe(viz.get_covariate_pc_correlations())

    # ─────────────────────────────────────────────────────────────────────────
    # TAB 3 — EMBEDDING & CLUSTERING
    # ─────────────────────────────────────────────────────────────────────────

    def _tab_embedding(self, global_kwargs: dict):
        if not _has_viz():
            st.info("Complete Upload & Preprocessing first.")
            return

        viz = _viz()

        if not viz.pp_state["pca_computed"]:
            st.warning("Run PCA in the Preprocessing tab first.")
            return

        # -- Neighbor / UMAP / Leiden controls --------------------------------
        with st.expander("⚙️ Embedding & Clustering Controls", expanded=not viz.pp_state["umap_computed"]):
            e1, e2 = st.columns(2)
            n_neighbors = e1.slider("Neighbors (k)", 5, 50, 15, 5, key="scp_n_neighbors")
            n_pcs_for_neighbors = e2.slider("PCs for neighbors", 5, 50, 30, 5, key="scp_n_pcs_nb")
            use_harmony = st.toggle(
                "Use Harmony-corrected PCA for neighbors",
                value=viz.pp_state["batch_corrected"],
                key="scp_use_harmony_nb"
            )
            e3, e4 = st.columns(2)
            umap_min_dist = e3.slider("UMAP min_dist", 0.01, 0.9, 0.3, 0.01, key="scp_umap_md")
            umap_spread = e4.slider("UMAP spread", 0.5, 2.0, 1.0, 0.1, key="scp_umap_sp")
            leiden_res = st.slider("Leiden resolution", 0.1, 3.0, 0.5, 0.1, key="scp_leiden_res")

            if st.button(
                "Run Neighbors → UMAP → Leiden",
                type="primary", use_container_width=True, key="scp_run_embed"
            ):
                with st.spinner("Computing neighborhood graph…"):
                    viz.run_neighbors(
                        n_neighbors=n_neighbors,
                        n_pcs=n_pcs_for_neighbors,
                        use_harmony=use_harmony
                    )
                with st.spinner("Computing UMAP…"):
                    viz.run_umap(min_dist=umap_min_dist, spread=umap_spread)
                with st.spinner("Running Leiden clustering…"):
                    n_cl = viz.run_leiden(resolution=leiden_res)
                for k in ["scp_pca", "scp_elbow", "scp_umap", "scp_cluster_comp"]:
                    st.session_state[f"{k}_fig"] = None
                st.success(f"✓ UMAP + Leiden done. Found {n_cl} clusters.")
                st.rerun()

        groupby_cols = viz.get_available_groupby_cols()
        color_by = st.selectbox(
            "Colour embedding by:", groupby_cols, key="scp_embed_color"
        ) if groupby_cols else None

        elbow_tab, pca_tab, umap_tab, expr_tab, cluster_tab = st.tabs([
            "📉 Elbow Plot", "🔵 PCA", "🌐 UMAP", "🎨 Expression Overlay", "🗂 Cluster Composition"
        ])

        with elbow_tab:
            pm = _pm("scp_elbow")
            pm.render_generate_button(viz.plot_elbow, **global_kwargs)
            pm.render_plot_and_editor()

        with pca_tab:
            pm = _pm("scp_pca")
            pm.render_generate_button(viz.plot_pca, color_by=color_by, **global_kwargs)
            pm.render_plot_and_editor()

        with umap_tab:
            if not viz.pp_state["umap_computed"]:
                st.info("Run embedding above to generate UMAP.")
            else:
                pm = _pm("scp_umap")
                pm.render_generate_button(viz.plot_umap, color_by=color_by, **global_kwargs)
                pm.render_plot_and_editor()

        with expr_tab:
            if not viz.pp_state["umap_computed"]:
                st.info("Run embedding above to generate UMAP first.")
            else:
                all_proteins = viz.get_protein_names()
                sel_protein = st.selectbox(
                    "Select protein to overlay:", all_proteins, key="scp_expr_protein"
                )
                pm = _pm("scp_expr_umap")
                pm.render_generate_button(
                    viz.plot_expression_umap, protein=sel_protein, **global_kwargs
                )
                pm.render_plot_and_editor()

        with cluster_tab:
            if not viz.pp_state["clustered"]:
                st.info("Run Leiden clustering above.")
            else:
                c1, c2 = st.columns(2)
                groupby_for_comp = c1.selectbox(
                    "Primary grouping:", ["leiden"] + groupby_cols, key="scp_comp_groupby"
                )
                splitby = c2.selectbox(
                    "Split by:", ["None"] + groupby_cols, key="scp_comp_splitby"
                )
                pm = _pm("scp_cluster_comp")
                pm.render_generate_button(
                    viz.plot_cluster_composition,
                    groupby=groupby_for_comp,
                    splitby=splitby if splitby != "None" else None,
                    **global_kwargs
                )
                pm.render_plot_and_editor()

    # ─────────────────────────────────────────────────────────────────────────
    # TAB 4 — DIFFERENTIAL EXPRESSION
    # ─────────────────────────────────────────────────────────────────────────

    def _tab_de(self, global_kwargs: dict):
        if not _has_viz():
            st.info("Complete Upload & Preprocessing first.")
            return

        viz = _viz()
        groupby_cols = viz.get_available_groupby_cols()

        # -- DE controls ------------------------------------------------------
        with st.expander("⚙️ Configure DE Analysis", expanded=not viz.pp_state["de_computed"]):
            dc1, dc2 = st.columns(2)
            de_groupby = dc1.selectbox("Group by:", groupby_cols, key="scp_de_groupby")
            if de_groupby:
                all_groups = viz.adata.obs[de_groupby].astype(str).unique().tolist()
                de_reference = dc2.selectbox(
                    "Reference group:", ["rest"] + all_groups, key="scp_de_reference"
                )
                de_groups = st.multiselect(
                    "Test groups (empty = all):", all_groups, key="scp_de_groups"
                )
                de3, de4, de5 = st.columns(3)
                de_method = de3.selectbox(
                    "Method:", ["wilcoxon", "t-test"], key="scp_de_method"
                )
                de_min_cells = de4.number_input(
                    "Min cells per group", 3, 100, 3, key="scp_de_min_cells"
                )
                de_min_det = de5.slider(
                    "Min detection % (protein filter)", 0, 100, 25, 5, key="scp_de_min_det"
                )

                if st.button(
                    "Run Differential Expression",
                    type="primary", use_container_width=True, key="scp_run_de"
                ):
                    with st.spinner("Running DE analysis (Wilcoxon rank-sum)…"):
                        try:
                            de_res = viz.run_de(
                                groupby=de_groupby,
                                groups=de_groups if de_groups else None,
                                reference=de_reference,
                                method=de_method,
                                min_cells=de_min_cells,
                                min_detection_pct=de_min_det / 100,
                            )
                            st.session_state.scp_de_results = de_res
                            for k in ["scp_volcano"]:
                                st.session_state[f"{k}_fig"] = None
                            st.success(f"✓ DE complete for {len(de_res)} group(s)")
                        except Exception as e:
                            st.error(f"DE failed: {e}")
                            logger.error(e, exc_info=True)
                    st.rerun()

        if not st.session_state.scp_de_results:
            st.info("Configure and run DE analysis above.")
            return

        de_groups_avail = viz.get_de_groups()

        sel_group = st.selectbox(
            "Select group to visualize:", de_groups_avail, key="scp_de_view_group"
        )
        if not sel_group:
            return

        de_df = viz.get_de_results(sel_group)

        # Significance controls
        sc1, sc2 = st.columns(2)
        fc_thresh = sc1.slider("log₂FC threshold", 0.0, 5.0, 1.0, 0.1, key="scp_vol_fc")
        pval_thresh = sc2.select_slider(
            "adj. p-value threshold",
            options=[0.001, 0.01, 0.05, 0.1, 0.2],
            value=0.05, key="scp_vol_pval"
        )

        sig_df = de_df[
            (de_df["pval_adj"] < pval_thresh) & (de_df["log2FC"].abs() > fc_thresh)
        ]
        m1, m2, m3 = st.columns(3)
        m1.metric("Significant proteins", len(sig_df))
        m2.metric("Up-regulated", (sig_df["log2FC"] > 0).sum())
        m3.metric("Down-regulated", (sig_df["log2FC"] < 0).sum())

        volcano_tab, table_tab, heatmap_tab = st.tabs([
            "🌋 Volcano Plot", "📋 Results Table", "🔥 Heatmap"
        ])

        with volcano_tab:
            n_label = st.slider("Max protein labels", 3, 25, 10, key="scp_n_label")
            pm = _pm("scp_volcano")
            pm.render_generate_button(
                viz.plot_volcano_sc, de_df=de_df,
                title=f"Volcano — {sel_group} vs {de_df['reference'].iloc[0]}",
                fc_thresh=fc_thresh, pval_thresh=pval_thresh, n_label=n_label,
                **global_kwargs
            )
            pm.render_plot_and_editor()

        with table_tab:
            show_only_sig = st.toggle("Show only significant proteins", value=True, key="scp_sig_only")
            display_df = sig_df if show_only_sig else de_df
            st.dataframe(
                display_df.sort_values("pval_adj")[
                    ["protein", "log2FC", "pval", "pval_adj",
                     "pct_group", "pct_ref", "protein_class"]
                ],
                use_container_width=True
            )
            csv = display_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download DE Results CSV", csv,
                file_name=f"DE_{sel_group}.csv",
                mime="text/csv", key="scp_dl_de"
            )

            if "protein_class" in de_df.columns:
                st.markdown("**Protein Classification Summary**")
                cls_count = sig_df["protein_class"].value_counts().reset_index()
                cls_count.columns = ["Class", "Count"]
                st.dataframe(cls_count, use_container_width=True)

        with heatmap_tab:
            n_top_heat = st.slider(
                "Top n proteins in heatmap", 5, 100, 25, 5, key="scp_heat_n"
            )
            if groupby_cols:
                heat_groupby = st.selectbox(
                    "Color columns by:", groupby_cols, key="scp_heat_groupby"
                )
            else:
                heat_groupby = None

            pm_heat = _mpl_pm("scp_de_heatmap")
            pm_heat.render_generate_button(
                viz.plot_de_heatmap,
                de_df=de_df,
                groupby=heat_groupby or de_groupby,
                n_top=n_top_heat,
                pval_thresh=pval_thresh,
                fc_thresh=fc_thresh,
            )
            pm_heat.render_plot_and_editor()

    # ─────────────────────────────────────────────────────────────────────────
    # TAB 5 — ACTIVITY SCORING
    # ─────────────────────────────────────────────────────────────────────────

    def _tab_activity_scoring(self, global_kwargs: dict):
        if not _has_viz():
            st.info("Complete Upload & Preprocessing first.")
            return

        viz = _viz()

        st.markdown("""
        **Per-cell activity scoring** uses `scanpy.tl.score_genes` to compute a summary
        expression score for each gene/protein set. Run DE → Run Enrichment → select
        pathways below → Compute Activity Scores.
        """)

        # -- From DE Results (GSEApy Enrichment) ------------------------------
        if st.session_state.scp_de_results:
            with st.expander("🔬 From DE Results (GSEApy Enrichment)", expanded=True):
                de_groups_avail = viz.get_de_groups()
                enr_group = st.selectbox("DE group:", de_groups_avail, key="scp_enr_group")

                ec1, ec2, ec3 = st.columns(3)
                enr_pval = ec1.select_slider(
                    "adj p-val", [0.001, 0.01, 0.05, 0.1], value=0.05, key="scp_enr_pval"
                )
                enr_fc = ec2.slider("|log₂FC|", 0.0, 5.0, 1.0, 0.5, key="scp_enr_fc")
                enr_dir = ec3.selectbox("Direction", ["both", "up", "down"], key="scp_enr_dir")

                _DB_OPTIONS = [
                    "KEGG_2021_Human",
                    "Reactome_2022",
                    "GO_Biological_Process_2023",
                    "WikiPathways_2019_Human",
                    "MSigDB_Hallmark_2020",
                ]
                enr_dbs = st.multiselect(
                    "Gene set databases:", _DB_OPTIONS,
                    default=["KEGG_2021_Human", "Reactome_2022"],
                    key="scp_enr_dbs",
                )

                # Show live preview of proteins that will pass the current filters
                _de_preview = viz.get_de_results(enr_group)
                if not _de_preview.empty:
                    _sig_preview = _de_preview[
                        (_de_preview["pval_adj"] < enr_pval)
                        & (_de_preview["log2FC"].abs() > enr_fc)
                    ]
                    if enr_dir == "up":
                        _sig_preview = _sig_preview[_sig_preview["log2FC"] > 0]
                    elif enr_dir == "down":
                        _sig_preview = _sig_preview[_sig_preview["log2FC"] < 0]
                    _n_sig = len(_sig_preview)
                    _preview_names = ", ".join(_sig_preview["protein"].dropna().head(5).tolist())
                    _gene_col_candidates = ["Genes", "Gene.Names", "Gene names", "Gene", "gene_names"]
                    _detected_gene_col = next(
                        (c for c in _gene_col_candidates
                         if c in viz.adata.var.columns
                         and viz.adata.var[c].notna().any()
                         and not viz.adata.var[c].astype(str).str.strip().str.lower().isin({"", "nan", "none", "na"}).all()),
                        None,
                    )
                    _gene_col_note = (
                        f" | Gene symbols from `{_detected_gene_col}`"
                        if _detected_gene_col
                        else " | ⚠️ No gene-name column detected — splitting protein IDs on `;`"
                    )
                    st.caption(
                        f"{_n_sig} protein group(s) pass current filters"
                        + (f" — e.g. {_preview_names}" if _preview_names else "")
                        + _gene_col_note
                    )

                if st.button("Run Enrichment", type="primary",
                             use_container_width=True, key="scp_run_enr"):
                    with st.spinner("Running GSEApy Enrichr…"):
                        try:
                            enr_df = viz.run_gsea_enrichment(
                                de_group=enr_group,
                                pval_thresh=enr_pval,
                                fc_thresh=enr_fc,
                                gene_sets=enr_dbs,
                                direction=enr_dir,
                            )
                            st.session_state.scp_enr_results = enr_df
                            if enr_df.empty:
                                st.warning(
                                    "Enrichment returned no results. "
                                    "If protein names above contain semicolons or UniProt IDs "
                                    "instead of gene symbols, consider using a gene-name column "
                                    "as your Protein ID when loading data. "
                                    "Also try relaxing the p-val / FC filters or choosing different databases."
                                )
                            else:
                                st.success(f"✓ {len(enr_df)} enriched terms found.")
                        except Exception as e:
                            st.error(f"Enrichment failed: {e}")

                enr_res = st.session_state.scp_enr_results
                if enr_res is not None and not enr_res.empty:
                    available_cols = enr_res.columns.tolist()
                    display_cols = [c for c in ["Gene_set", "Term", "Overlap", "Adjusted P-value", "Genes"] if c in available_cols]
                    enr_show = enr_res[display_cols].sort_values("Adjusted P-value") if "Adjusted P-value" in display_cols else enr_res[display_cols]

                    _fc1, _fc2 = st.columns([3, 1])
                    with _fc1:
                        _term_search = st.text_input(
                            "Search terms", value="",
                            placeholder="Type to filter pathway names…",
                            key="scp_enr_term_search",
                        )
                    with _fc2:
                        _min_genes = st.number_input(
                            "Min. overlap genes", min_value=0, value=3, step=1,
                            key="scp_enr_min_genes",
                        )
                    if _term_search:
                        enr_show = enr_show[
                            enr_show["Term"].str.contains(_term_search, case=False, na=False)
                        ]
                    if _min_genes > 0 and "Overlap" in enr_show.columns:
                        try:
                            _overlap_n = enr_show["Overlap"].str.split("/").str[0].astype(int)
                            enr_show = enr_show[_overlap_n >= _min_genes]
                        except (ValueError, AttributeError):
                            pass

                    st.dataframe(enr_show, use_container_width=True)

                    selected_terms = st.multiselect(
                        "Select pathways to use as activity score inputs:",
                        enr_show["Term"].tolist(),
                        key="scp_enr_selected_terms",
                    )

                    if selected_terms:
                        # Preview how many genes from each selected pathway match the data.
                        # If the gene map only covers DEPs (<50%), note that full annotation
                        # will run automatically when the button is clicked.
                        sym_to_var = viz._build_gene_symbol_map()
                        available_prots = set(viz.adata.var_names)
                        gene_map_sparse = len(sym_to_var) < viz.adata.n_vars * 0.5
                        if gene_map_sparse:
                            st.caption(
                                f"ℹ️ Gene names only resolved for {len(sym_to_var)}/{viz.adata.n_vars} proteins so far. "
                                "Full UniProt annotation will run automatically when you click the button below."
                            )
                        for term in selected_terms:
                            row = enr_res[enr_res["Term"] == term].iloc[0]
                            genes = [g.strip() for g in str(row["Genes"]).split(";") if g.strip()]
                            matched = [g for g in genes if g in available_prots or g in sym_to_var]
                            st.caption(
                                f"**{term[:50]}**: {len(matched)}/{len(genes)} genes matched"
                                + (" (more expected after full annotation)" if gene_map_sparse and len(matched) == 0 else "")
                                + (" ⚠️ < 3 — will be skipped" if len(matched) < 3 and not gene_map_sparse else "")
                            )

                    if selected_terms and st.button(
                        "Compute Activity Scores from Selected Pathways",
                        type="primary",
                        use_container_width=True, key="scp_enr_add"
                    ):
                        # Gene symbol map may only cover DEPs from the enrichment run.
                        # If it covers <50% of proteins, fetch all gene names first so
                        # every pathway gene can be resolved to a var_name.
                        sym_to_var_check = viz._build_gene_symbol_map()
                        if len(sym_to_var_check) < viz.adata.n_vars * 0.5:
                            with st.spinner(
                                f"Fetching gene names from UniProt for all "
                                f"{viz.adata.n_vars} proteins first…"
                            ):
                                viz.annotate_gene_names_from_uniprot()

                        pathway_gene_sets = {}
                        for _, row in enr_res[enr_res["Term"].isin(selected_terms)].iterrows():
                            genes = [g.strip() for g in str(row["Genes"]).split(";") if g.strip()]
                            short_name = row["Term"][:30].replace(" ", "_").replace("/", "_")
                            pathway_gene_sets[short_name] = genes

                        with st.spinner("Computing per-cell activity scores…"):
                            computed = viz.compute_activity_scores(pathway_gene_sets)
                        st.session_state.scp_score_cols = computed
                        for k in ["scp_activity_umap", "scp_activity_violin"]:
                            st.session_state[f"{k}_fig"] = None
                        if computed:
                            st.success(f"✓ Scored: {', '.join(computed)}")
                        else:
                            st.error("No pathways had ≥3 matching proteins in your data.")
                        st.rerun()
        else:
            st.info("Run Differential Expression (Tab 4) first to enable pathway enrichment here.")

        # -- Score summary ----------------------------------------------------
        existing_scores = viz.get_available_score_cols()
        if not existing_scores:
            st.info("Compute activity scores above to proceed.")
            return

        st.markdown("---")
        st.markdown("**Score Summary**")
        score_stats = {
            s: {
                "mean": viz.adata.obs[s].mean(),
                "std": viz.adata.obs[s].std(),
                "min": viz.adata.obs[s].min(),
                "max": viz.adata.obs[s].max(),
            }
            for s in existing_scores
        }
        st.dataframe(pd.DataFrame(score_stats).T.round(3))

        # -- Per-score visualisations -----------------------------------------
        groupby_cols = viz.get_available_groupby_cols()
        score_col = st.selectbox(
            "Select score to plot:", existing_scores, key="scp_score_sel"
        )
        groupby_for_plots = st.selectbox(
            "Group by:", groupby_cols, key="scp_score_groupby"
        ) if groupby_cols else None
        splitby_for_violin = st.selectbox(
            "Split violin by (optional):", ["None"] + groupby_cols, key="scp_score_splitby"
        )

        umap_score_tab, violin_score_tab, panel_tab = st.tabs([
            "🌐 UMAP Activity", "🎻 Violin by Group", "📊 Multi-Score Panel"
        ])

        with umap_score_tab:
            if not viz.pp_state["umap_computed"]:
                st.info("Run UMAP (in the Embedding tab) first.")
            else:
                pm = _pm("scp_activity_umap")
                pm.render_generate_button(
                    viz.plot_activity_umap, score_col=score_col, **global_kwargs
                )
                pm.render_plot_and_editor()

        with violin_score_tab:
            if not groupby_for_plots:
                st.info("No grouping columns available.")
            else:
                pm = _pm("scp_activity_violin")
                pm.render_generate_button(
                    viz.plot_activity_violin,
                    score_col=score_col,
                    groupby=groupby_for_plots,
                    splitby=splitby_for_violin if splitby_for_violin != "None" else None,
                    **global_kwargs
                )
                pm.render_plot_and_editor()

        with panel_tab:
            st.markdown(
                "Static multi-panel figure: UMAP (if available) + violin "
                "for all computed scores."
            )
            selected_panel_scores = st.multiselect(
                "Scores to include:", existing_scores,
                default=existing_scores[:min(4, len(existing_scores))],
                key="scp_panel_scores"
            )
            panel_groupby = st.selectbox(
                "Group by (violin row):", groupby_cols, key="scp_panel_groupby"
            ) if groupby_cols else None

            if st.button(
                "Generate Panel Figure", use_container_width=True, key="scp_gen_panel"
            ):
                if not selected_panel_scores:
                    st.warning("Select at least one score.")
                elif not panel_groupby:
                    st.warning("Select a groupby column.")
                else:
                    try:
                        with st.spinner("Generating multi-score panel…"):
                            buf = viz.plot_activity_scores_panel(
                                selected_panel_scores, panel_groupby
                            )
                        st.image(buf, caption="Per-Cell Activity Score Panel")
                    except Exception as e:
                        st.error(f"Panel generation failed: {e}")
                        logger.error(e, exc_info=True)

    # ─────────────────────────────────────────────────────────────────────────
    # MAIN RENDER
    # ─────────────────────────────────────────────────────────────────────────

    def render(self):
        st.header("🔬 Single-Cell Proteomics Analysis")

        if _has_viz():
            summary = _viz().get_preprocessing_summary()
            st.caption(
                f"Active dataset: **{summary['n_cells']} cells × {summary['n_proteins']} proteins** | "
                f"Layers: {', '.join(summary['layers'])} | "
                f"Embeddings: {', '.join(summary['embeddings']) or 'none'}"
            )

        global_kwargs = self._render_global_settings()

        (
            tab_upload, tab_preprocess, tab_embed, tab_de, tab_score
        ) = st.tabs([
            "1️⃣  Upload & QC",
            "2️⃣  Preprocessing",
            "3️⃣  Embedding & Clustering",
            "4️⃣  Differential Expression",
            "5️⃣  Activity Scoring",
        ])

        with tab_upload:
            self._tab_upload_qc(global_kwargs)

        with tab_preprocess:
            self._tab_preprocessing(global_kwargs)

        with tab_embed:
            self._tab_embedding(global_kwargs)

        with tab_de:
            self._tab_de(global_kwargs)

        with tab_score:
            self._tab_activity_scoring(global_kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# Module entry point (matches Pro-Visualize pattern)
# ─────────────────────────────────────────────────────────────────────────────

def render():
    SCPTab().render()
