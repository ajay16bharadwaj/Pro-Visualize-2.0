import streamlit as st
import pandas as pd
import logging
import plotly.express as px
from visualizations.dilution_series import DilutionSeriesVisualizer
from utils.helpers import handle_plotting_errors
from utils.plot_manager import PlotManager

logger = logging.getLogger(__name__)


class DilutionSeriesTab:
    PLOT_KEYS = [
        "dilution_counts", "dilution_dist", "dilution_cv",
        "dilution_trends", "dilution_heatmap", "dilution_ratios",
        "dilution_pca", "dilution_completeness",
        "dilution_r2_hist", "dilution_lod_loq",
    ]

    def __init__(self):
        if 'dilution_visualizer' not in st.session_state:
            st.session_state.dilution_visualizer = None

    def _upload_data_section(self):
        with st.expander("Upload Data for Dilution Series",
                         expanded=(st.session_state.dilution_visualizer is None)):
            protein_file = st.file_uploader(
                "Upload Protein-Level Data (.txt, .csv, .tsv)", type=['txt', 'csv', 'tsv'],
                key="dilution_protein_uploader")
            metadata_file = st.file_uploader(
                "Upload Metadata File (.csv, .tsv)", type=['csv', 'tsv'],
                key="dilution_metadata_uploader")

            st.markdown("**⚙️ Column Configuration**")
            st.caption("Map your file's column names to the roles this module expects. "
                       "Leave defaults if your files already use these names.")
            with st.container():
                c1, c2 = st.columns(2)
                protein_id_col = c1.text_input("Protein ID Column", value="Protein.Group",
                                               key="dil_col_protein_id")
                gene_col = c2.text_input("Gene Symbol Column", value="Genes",
                                         key="dil_col_gene")
                c3, c4 = st.columns(2)
                sample_col = c3.text_input("Sample Column (metadata)", value="Sample",
                                           key="dil_col_sample")
                concentration_col = c4.text_input("Concentration Column (metadata)", value="Concentration",
                                                  key="dil_col_conc")
                c5, c6 = st.columns(2)
                replicate_col = c5.text_input("Replicate Column (metadata)", value="Replicate",
                                              key="dil_col_rep")
                group_col = c6.text_input("Group Column (metadata)", value="Group",
                                          key="dil_col_group")

            if st.button("Load & Process Dilution Data", use_container_width=True):
                if protein_file and metadata_file:
                    try:
                        protein_df = pd.read_csv(protein_file, sep=None, engine='python')
                        metadata_df = pd.read_csv(metadata_file, sep=None, engine='python')
                        with st.spinner("Initializing visualizer..."):
                            st.session_state.dilution_visualizer = DilutionSeriesVisualizer(
                                protein_df, metadata_df,
                                protein_id_col=protein_id_col,
                                gene_col=gene_col,
                                sample_col=sample_col,
                                concentration_col=concentration_col,
                                replicate_col=replicate_col,
                                group_col=group_col,
                            )
                        st.success("Data loaded successfully!")
                        st.rerun()
                    except Exception as e:
                        st.error(f"An error occurred: {e}")
                        logger.error(f"Dilution Series Error: {e}", exc_info=True)
                else:
                    st.warning("Please upload both data and metadata files.")

    def _render_global_settings(self):
        """Renders global plot styling and deviation threshold controls."""
        with st.expander("🎨 Global Plot Settings (Colors, Theme & Thresholds)", expanded=False):
            c1, c2 = st.columns(2)
            with c1:
                theme_options = {
                    "Standard White": "plotly_white",
                    "Dark Mode": "plotly_dark",
                    "Minimal": "simple_white",
                    "GGPlot Style": "ggplot2",
                    "Seaborn Style": "seaborn"
                }
                selected_theme_label = st.selectbox(
                    "Plot Theme", list(theme_options.keys()), key="dilution_global_theme_select")
                selected_template = theme_options[selected_theme_label]

            with c2:
                palette_mode = st.radio(
                    "Palette Type", ["Auto", "Colorblind Safe", "Custom"],
                    horizontal=True, key="dilution_global_palette_mode")

            visualizer = st.session_state.dilution_visualizer
            groups = visualizer.group_order
            color_map = {}

            if palette_mode == "Colorblind Safe":
                safe_colors = px.colors.qualitative.Safe
                for i, group in enumerate(groups):
                    color_map[group] = safe_colors[i % len(safe_colors)]
                st.caption("✅ Using high-contrast, colorblind-safe colors.")
            elif palette_mode == "Custom":
                st.markdown("**Customize Group Colors:**")
                cols = st.columns(4)
                for i, group in enumerate(groups):
                    default_c = px.colors.qualitative.Safe[i % len(px.colors.qualitative.Safe)]
                    with cols[i % 4]:
                        color_map[group] = st.color_picker(
                            f"{group}", value=default_c, key=f"dilution_color_{group}")
            else:
                color_map = None

            st.markdown("---")
            st.markdown("**Deviation Bucket Thresholds** (used in Ratio plot color coding)")
            dc1, dc2 = st.columns(2)
            good_thresh = dc1.slider(
                "Good threshold (|dev| < x → green)", 0.05, 0.50, 0.20, 0.05,
                key="dilution_dev_good",
                help="Mean deviation below this value is classified as excellent linearity.")
            warn_thresh = dc2.slider(
                "Warn threshold (|dev| < x → orange, else red)", good_thresh + 0.05, 1.0, 0.50, 0.05,
                key="dilution_dev_warn",
                help="Mean deviation below this (but above good threshold) is moderate; above is poor.")

            st.markdown("---")
            if st.button("Apply Settings to All Plots", type="primary", use_container_width=True):
                count = 0
                for key in self.PLOT_KEYS:
                    fig_key = f"{key}_fig"
                    if st.session_state.get(fig_key) is not None:
                        st.session_state[fig_key] = None
                        count += 1
                if count > 0:
                    st.success(f"Settings applied! {count} plot(s) will regenerate.")
                else:
                    st.info("Settings saved. Generate plots to see changes.")
                st.rerun()

        global_kwargs = {"template": selected_template}
        if color_map:
            global_kwargs["color_discrete_map"] = color_map

        return global_kwargs, good_thresh, warn_thresh

    @handle_plotting_errors
    def _display_plots(self):
        visualizer = st.session_state.dilution_visualizer
        st.subheader("Dilution Series Visualizations")

        global_plot_kwargs, good_thresh, warn_thresh = self._render_global_settings()

        # --- Sanity checks ---
        with st.expander("🔍 Data Quality Checks", expanded=True):
            try:
                checks = visualizer.run_sanity_checks()
                if not checks:
                    st.success("✓ All sanity checks passed.")
                for check in checks:
                    if check['level'] == 'error':
                        st.error(check['message'])
                    elif check['level'] == 'warning':
                        st.warning(check['message'])
                    else:
                        st.info(check['message'])
            except Exception as e:
                st.warning(f"Sanity checks could not complete: {e}")

        (overview_tab, analysis_tab, pca_tab, completeness_tab,
         r2_tab, lod_tab) = st.tabs([
            "📊 Data Overview", "📈 Trend Analysis", "🧬 PCA",
            "✅ Completeness", "📉 Linearity (R²)", "🔬 LOD/LOQ"
        ])

        # --- TAB 1: OVERVIEW ---
        with overview_tab:
            try:
                st.markdown("### Protein Counts per Sample")
                counts_manager = PlotManager("dilution_counts")
                counts_manager.module = "dilution"
                counts_manager.render_generate_button(
                    visualizer.plot_protein_counts_per_sample, **global_plot_kwargs)
                counts_manager.render_plot_and_editor()

                st.markdown("---")
                st.markdown("### Intensity and CV Distributions")
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown("**Intensity Distribution**")
                    with st.expander("⚙️ Configuration", expanded=True):
                        plot_type = st.selectbox("Intensity Plot Type", ["box", "violin"],
                                                 key="dist_plot_type")
                        dist_manager = PlotManager("dilution_dist")
                        dist_manager.module = "dilution"
                        dist_manager.render_generate_button(
                            visualizer.plot_intensity_distribution,
                            plot_type=plot_type, points=False, **global_plot_kwargs)
                    dist_manager.render_plot_and_editor()

                with col2:
                    st.markdown("**CV Distribution**")
                    cv_manager = PlotManager("dilution_cv")
                    cv_manager.module = "dilution"
                    cv_manager.render_generate_button(
                        visualizer.plot_cv_distribution, **global_plot_kwargs)
                    cv_manager.render_plot_and_editor()
            except Exception as e:
                st.error(f"Overview tab error: {e}")
                logger.error(f"Dilution overview tab: {e}", exc_info=True)

        # --- TAB 2: ANALYSIS ---
        with analysis_tab:
            trends_tab, heatmap_tab, ratio_tab = st.tabs(
                ["Protein Trends", "Heatmap", "Relative Abundance"])

            with trends_tab:
                try:
                    st.markdown("### Individual Protein Trends")
                    with st.expander("⚙️ Configuration", expanded=True):
                        all_proteins = visualizer.protein_df[visualizer.protein_id_col].unique().tolist()
                        selected_proteins = st.multiselect(
                            "Select specific proteins to plot (optional)", options=all_proteins,
                            help="If none selected, the most abundant proteins are shown.",
                            key="dilution_trends_protein_selector")
                        n_top = st.slider(
                            "Number of top proteins to show (if none selected above)", 1, 20, 5,
                            disabled=bool(selected_proteins), key="dilution_trends_top_n_slider")
                        trends_manager = PlotManager("dilution_trends")
                        trends_manager.module = "dilution"
                        trends_manager.render_generate_button(
                            visualizer.plot_protein_trends,
                            proteins_to_plot=selected_proteins if selected_proteins else None,
                            n_top_proteins=n_top, **global_plot_kwargs)
                    trends_manager.render_plot_and_editor()
                except Exception as e:
                    st.error(f"Protein trends error: {e}")

            with heatmap_tab:
                try:
                    st.markdown("### Protein Intensity Heatmap")
                    with st.expander("⚙️ Configuration", expanded=True):
                        apply_zscore = st.toggle("Apply Z-score normalization", value=True,
                                                 key="dilution_heat_z")
                        max_proteins = st.slider("Max proteins to display", 50, 2000, 500, 50,
                                                 key="dilution_heat_max")
                        heat_manager = PlotManager("dilution_heatmap")
                        heat_manager.module = "dilution"
                        heat_manager.render_generate_button(
                            visualizer.plot_heatmap_trends,
                            max_proteins_to_plot=max_proteins,
                            apply_zscore=apply_zscore, **global_plot_kwargs)
                    heat_manager.render_plot_and_editor()
                except Exception as e:
                    st.error(f"Heatmap error: {e}")

            with ratio_tab:
                try:
                    st.markdown("### Protein Abundance Ratios")
                    st.markdown("Shows `Log2` intensity ratio relative to lowest concentration.")
                    with st.expander("⚙️ Configuration", expanded=True):
                        rc1, rc2 = st.columns(2)
                        show_lines = rc1.toggle("Show expected ratio lines", value=True,
                                                key="dilution_ratio_lines")
                        show_deviations = rc2.toggle("Show deviation metrics", value=True,
                                                     key="dilution_ratio_deviations",
                                                     help="Mean deviation with ±SD error bars")
                        ratio_manager = PlotManager("dilution_ratios")
                        ratio_manager.module = "dilution"
                        ratio_manager.render_generate_button(
                            visualizer.plot_relative_abundance_ratios,
                            add_expected_lines=show_lines,
                            show_deviations=show_deviations,
                            good_thresh=good_thresh,
                            warn_thresh=warn_thresh,
                            points=False, **global_plot_kwargs)
                    ratio_manager.render_plot_and_editor()

                    if show_deviations:
                        with st.expander("📊 Interpreting Deviation Metrics", expanded=False):
                            st.markdown(f"""
                            **Color Coding (current thresholds):**
                            - 🟢 **Green:** Mean deviation < {good_thresh} (Excellent agreement with expected)
                            - 🟠 **Orange:** Mean deviation {good_thresh}–{warn_thresh} (Moderate deviation)
                            - 🔴 **Red:** Mean deviation ≥ {warn_thresh} (Poor linearity)
                            """)

                    # CSV export: CV-by-concentration matrix
                    st.markdown("---")
                    st.markdown("**Export: CV% by Concentration Group**")
                    if st.button("Compute CV Matrix", key="dil_cv_matrix_btn"):
                        try:
                            cv_matrix = visualizer.get_cv_by_concentration_matrix()
                            st.session_state['dilution_cv_matrix'] = cv_matrix
                        except Exception as e:
                            st.error(f"Could not compute CV matrix: {e}")
                    if 'dilution_cv_matrix' in st.session_state:
                        cv_m = st.session_state['dilution_cv_matrix']
                        st.dataframe(cv_m.head(20), use_container_width=True)
                        st.download_button(
                            "⬇ Download CV Matrix CSV",
                            data=cv_m.to_csv(),
                            file_name="dilution_cv_by_concentration.csv",
                            mime="text/csv",
                            key="dil_cv_matrix_dl")
                except Exception as e:
                    st.error(f"Ratio tab error: {e}")

        # --- TAB 3: PCA ---
        with pca_tab:
            try:
                st.markdown("### Principal Component Analysis (PCA)")
                with st.expander("⚙️ Configuration", expanded=True):
                    c1, c2 = st.columns(2)
                    color_by = c1.selectbox("Color PCA points by",
                                            options=['Group', 'Concentration', 'Replicate'],
                                            index=0, key="dilution_pca_color")
                    symbol_by = c2.selectbox("Use symbols for",
                                             options=[None, 'Replicate', 'Group'],
                                             index=1, key="dilution_pca_symbol")
                    pca_manager = PlotManager("dilution_pca")
                    pca_manager.module = "dilution"
                    pca_manager.render_generate_button(
                        visualizer.plot_pca,
                        color_by=color_by, symbol_by=symbol_by, **global_plot_kwargs)
                pca_manager.render_plot_and_editor()
            except Exception as e:
                st.error(f"PCA tab error: {e}")

        # --- TAB 4: COMPLETENESS ---
        with completeness_tab:
            try:
                st.markdown("### Detection Completeness Analysis")
                with st.expander("⚙️ Configuration", expanded=True):
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        analysis_level = st.selectbox(
                            "Select Analysis Level",
                            options=["Proteins", "Precursors", "Stripped Sequences"],
                            key="dilution_comp_level")
                        level_mapping = {
                            "Proteins": visualizer.protein_id_col,
                            "Precursors": ("Precursor.Id" if "Precursor.Id" in visualizer.protein_df.columns
                                           else visualizer.protein_id_col),
                            "Stripped Sequences": ("Stripped.Sequence" if "Stripped.Sequence" in visualizer.protein_df.columns
                                                   else visualizer.protein_id_col),
                        }
                        selected_column = level_mapping[analysis_level]
                    with col2:
                        cv_thresh = st.slider("CV Threshold (%)", 0.0, 50.0, 20.0, 1.0,
                                              key="dilution_comp_cv")
                    with col3:
                        use_log = st.toggle("Use logarithmic scale", value=False,
                                            key="dilution_comp_log")

                    comp_manager = PlotManager("dilution_completeness")
                    comp_manager.module = "dilution"
                    comp_manager.render_generate_button(
                        visualizer.plot_completeness_overview,
                        identifier_col=selected_column,
                        use_log_scale=use_log,
                        cv_threshold=cv_thresh, **global_plot_kwargs)
                comp_manager.render_plot_and_editor()

                # CSV export: completeness summary
                st.markdown("---")
                st.markdown("**Export: Detection Completeness Summary**")
                if st.button("Compute Completeness Summary", key="dil_comp_summary_btn"):
                    try:
                        comp_df = visualizer.get_completeness_summary(cv_threshold=cv_thresh)
                        st.session_state['dilution_completeness_df'] = comp_df
                    except Exception as e:
                        st.error(f"Could not compute completeness summary: {e}")
                if 'dilution_completeness_df' in st.session_state:
                    comp_df = st.session_state['dilution_completeness_df']
                    st.dataframe(comp_df, use_container_width=True, hide_index=True)
                    st.download_button(
                        "⬇ Download Completeness CSV",
                        data=comp_df.to_csv(index=False),
                        file_name="dilution_completeness_summary.csv",
                        mime="text/csv",
                        key="dil_comp_dl")
            except Exception as e:
                st.error(f"Completeness tab error: {e}")

        # --- TAB 5: LINEARITY (R²) ---
        with r2_tab:
            try:
                st.markdown("### Protein Linearity — R² Distribution")
                st.caption("Measures how well each protein's log2 intensity follows "
                           "log2 concentration (ideal dilution behavior).")

                with st.expander("⚙️ Configuration", expanded=True):
                    n_bins = st.slider("Histogram bins", 10, 50, 20, key="dil_r2_bins")
                    r2_hist_manager = PlotManager("dilution_r2_hist")
                    r2_hist_manager.module = "dilution"
                    r2_hist_manager.render_generate_button(
                        visualizer.plot_r2_histogram, n_bins=n_bins, **global_plot_kwargs)
                r2_hist_manager.render_plot_and_editor()

                st.markdown("---")
                st.markdown("### Ranked R² Table")
                if st.button("Load R² Table", key="dil_r2_table_btn"):
                    try:
                        r2_df = visualizer.get_r2_table()
                        st.session_state['dilution_r2_df'] = r2_df
                    except Exception as e:
                        st.error(f"Could not compute R² table: {e}")
                if 'dilution_r2_df' in st.session_state:
                    r2_df = st.session_state['dilution_r2_df']
                    r2_sort = r2_df.sort_values('r_squared', ascending=False)
                    st.caption(f"{len(r2_df)} proteins fitted | "
                               f"Median R² = {r2_df['r_squared'].median():.3f} | "
                               f"R² ≥ 0.95: {(r2_df['r_squared'] >= 0.95).sum()} proteins")
                    st.dataframe(r2_sort.head(200), use_container_width=True, hide_index=True)
                    st.download_button(
                        "⬇ Download Full R² Table CSV",
                        data=r2_sort.to_csv(index=False),
                        file_name="dilution_r2_table.csv",
                        mime="text/csv",
                        key="dil_r2_dl")
            except Exception as e:
                st.error(f"Linearity tab error: {e}")

        # --- TAB 6: LOD/LOQ ---
        with lod_tab:
            try:
                st.markdown("### LOD / LOQ Estimation per Protein")
                st.caption(
                    "LOD and LOQ are estimated from the log-log linear calibration fit:\n\n"
                    "**LOD** = 2^(log₂(C_min) + 3.3·σ / slope) | "
                    "**LOQ** = 2^(log₂(C_min) + 10·σ / slope)\n\n"
                    "where σ is the residual standard deviation and slope is from "
                    "log₂(Intensity) ~ log₂(Concentration). Proteins with negative slope "
                    "are excluded (intensity does not increase with concentration).")

                with st.expander("⚙️ Configuration", expanded=True):
                    top_n = st.slider("Proteins to display in plot (ranked by R²)", 10, 200, 50,
                                      key="dil_lod_top_n")
                    lod_manager = PlotManager("dilution_lod_loq")
                    lod_manager.module = "dilution"
                    lod_manager.render_generate_button(
                        visualizer.plot_lod_loq, top_n=top_n, **global_plot_kwargs)
                lod_manager.render_plot_and_editor()

                st.markdown("---")
                st.markdown("### LOD/LOQ Summary Table")
                if st.button("Load LOD/LOQ Table", key="dil_lod_table_btn"):
                    try:
                        lod_df = visualizer.get_lod_loq_table()
                        st.session_state['dilution_lod_df'] = lod_df
                    except Exception as e:
                        st.error(f"Could not compute LOD/LOQ table: {e}")
                if 'dilution_lod_df' in st.session_state:
                    lod_df = st.session_state['dilution_lod_df']
                    valid_count = lod_df['LOD (ng)'].notna().sum()
                    in_range_count = lod_df['LOD in range'].sum()
                    st.caption(f"{valid_count} of {len(lod_df)} proteins have valid LOD estimates | "
                               f"{in_range_count} have LOD within the tested concentration range "
                               f"(≤ {visualizer.metadata_df['Concentration'].max()} ng)")
                    display_cols = [c for c in [
                        visualizer.protein_id_col, visualizer.gene_col,
                        'R²', 'slope', 'residual_std', 'LOD (ng)', 'LOQ (ng)', 'LOD in range'
                    ] if c and c in lod_df.columns]
                    st.dataframe(lod_df[display_cols].head(500),
                                 use_container_width=True, hide_index=True)
                    st.download_button(
                        "⬇ Download LOD/LOQ Table CSV",
                        data=lod_df[display_cols].to_csv(index=False),
                        file_name="dilution_lod_loq.csv",
                        mime="text/csv",
                        key="dil_lod_dl")
            except Exception as e:
                st.error(f"LOD/LOQ tab error: {e}")

    def render(self):
        st.header("Dilution Series Analysis")
        self._upload_data_section()
        if st.session_state.dilution_visualizer:
            self._display_plots()


def render():
    """Renders the Dilution Series tab."""
    dilution_tab = DilutionSeriesTab()
    dilution_tab.render()
