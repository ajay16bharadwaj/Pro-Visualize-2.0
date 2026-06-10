import streamlit as st
import pandas as pd
import logging
import plotly.express as px
from plotly import colors as pc
from visualizations.comparative_visualizer import ComparativeVisualizer
from utils.plot_manager import PlotManager, MplPlotManager
from config.plot_configs import SIGNIFICANCE_PRESETS

logger = logging.getLogger(__name__)

class ComparativeTab:
    """
    Encapsulates the UI and logic for the Comparative Analysis tab.
    """
    PLOT_KEYS = [
        "comp_volcano", "comp_violin", "comp_manhattan", "comp_dotplot"
    ]

    def __init__(self):
        if 'comp_visualizer' not in st.session_state:
            st.session_state.comp_visualizer = None

    def _upload_and_config_section(self):
        with st.expander("Upload & Configure Data", expanded=(st.session_state.comp_visualizer is None)):
            st.markdown("##### 1. Upload Files")
            protein_file = st.file_uploader("Upload Protein-Level Data", type=['txt', 'csv', 'tsv'], key="comp_protein_uploader")
            annotation_file = st.file_uploader("Upload Annotation File", type=['txt', 'csv', 'tsv'], key="comp_annotation_uploader")
            comparative_file = st.file_uploader("Upload Comparative Analysis File", type=['txt', 'csv', 'tsv'], key="comp_file_uploader")

            st.markdown("---")
            st.markdown("##### 2. Configure Column Names")
            c1, c2 = st.columns(2)
            with c1:
                protein_id_col = st.text_input("Protein ID Column", value="Protein", key="comp_prot_id")
                sample_id_col = st.text_input("Sample Linking Column", value="Level3", key="comp_samp_id")
                grouping_col = st.text_input("Experimental Group Column", value="attribute_ExperimentalGroup", key="comp_group_id")
            with c2:
                comp_prot_id = st.text_input("Protein ID (Comp. File)", value="Protein", key="comp_c_prot_id")
                fc_col = st.text_input("Log2 Fold Change Column", value="log2FC", key="comp_fc_id")
                fdr_col = st.text_input("FDR / p-value Column", value="Imputed.FDR", key="comp_fdr_id")
                lbl_col = st.text_input("Comparison Label Column", value="Label2", key="comp_lbl_id")

            if st.button("Load & Process Data", use_container_width=True, type="primary", key="comp_load_data_btn"):
                if all([protein_file, annotation_file, comparative_file]):
                    try:
                        p_df = pd.read_csv(protein_file, sep=None, engine='python')
                        a_df = pd.read_csv(annotation_file, sep=None, engine='python')
                        c_df = pd.read_csv(comparative_file, sep=None, engine='python')
                        config = {
                            "protein_id": protein_id_col, "sample_id": sample_id_col,
                            "grouping": grouping_col, "comp_protein_id": comp_prot_id,
                            "fold_change": fc_col, "fdr": fdr_col, "comparison_label": lbl_col
                        }
                        with st.spinner("Initializing..."):
                            st.session_state.comp_visualizer = ComparativeVisualizer(p_df, a_df, c_df, config)
                        st.success("Data loaded!"); st.rerun()
                    except Exception as e: st.error(f"Error: {e}")
                else: st.warning("Please upload all files.")

    def _render_global_settings(self):
        with st.expander("🎨 Global Plot Settings (Colors & Theme)", expanded=False):
            c1, c2 = st.columns(2)
            with c1:
                theme_options = {"Standard White": "plotly_white", "Dark Mode": "plotly_dark", "Minimal": "simple_white"}
                selected_theme = theme_options[st.selectbox("Plot Theme", list(theme_options.keys()), key="comp_theme")]
            with c2:
                vis = st.session_state.comp_visualizer
                grp_col = vis.column_config['grouping']
                groups = vis.annotation_df[grp_col].unique().tolist() if vis else []
                color_map = {}
                if groups:
                    st.markdown("**Customize Group Colors:**")
                    cols = st.columns(4)

                    try:
                        safe_colors = ['#%02x%02x%02x' % tuple(int(x) for x in pc.unlabel_rgb(c)) for c in px.colors.qualitative.Safe]
                    except Exception:
                        safe_colors = px.colors.qualitative.Safe

                    for i, group in enumerate(groups):
                        d_c = safe_colors[i % len(safe_colors)]
                        with cols[i % 4]: color_map[group] = st.color_picker(f"{group}", value=d_c, key=f"comp_color_{group}")

            if st.button("Apply Settings to All Plots", type="primary", use_container_width=True, key="comp_apply_settings_btn"):
                for key in list(st.session_state.keys()):
                    if key.startswith("comp_") and key.endswith("_fig"):
                        st.session_state[key] = None
                st.success("Settings applied! Plots will regenerate automatically."); st.rerun()

        kwargs = {"template": selected_theme}
        if color_map: kwargs["color_discrete_map"] = color_map
        return kwargs

    def render(self):
        self._upload_and_config_section()
        if 'significant_proteins' not in st.session_state: st.session_state.significant_proteins = pd.DataFrame()

        if st.session_state.comp_visualizer:
            visualizer = st.session_state.comp_visualizer
            st.header("Comparative Analysis")
            global_plot_kwargs = self._render_global_settings()

            # Read slider / selection state before entering tabs so that a
            # failure in one tab (e.g. select) does not NameError another tab.
            fdr = st.session_state.get("comp_fdr", 0.05)
            fc = st.session_state.get("comp_fc", 1.0)
            sig_df = st.session_state.get("significant_proteins", pd.DataFrame())
            pid = visualizer.column_config['comp_protein_id']

            overview, select, volcano, heat, expr, path = st.tabs(["Overview", "Selection", "Volcano", "Heatmap", "Expression", "Pathway"])

            with overview:
                try:
                    st.markdown("#### Data Preview")
                    st.markdown("**Protein Data**")
                    st.dataframe(visualizer.get_protein_data_preview().head())
                    st.markdown("**Annotation Data**")
                    st.dataframe(visualizer.get_annotation_data_preview().head())
                    st.markdown("**Comparative Data**")
                    st.dataframe(visualizer.get_comparative_data_preview().head())
                except Exception as _e:
                    logger.error("Comparative Overview tab error", exc_info=True)
                    st.error(f"**Overview tab** failed: {_e}")

            with select:
                try:
                    comps = visualizer.get_comparison_groups()
                    if 'selected_comparison' not in st.session_state: st.session_state.selected_comparison = comps[0] if comps else None
                    sel_comp = st.selectbox("Select Comparison:", comps, key='sel_comp_main')

                    st.markdown("**Significance Thresholds**")
                    c_pre, c_app = st.columns([3, 1])
                    with c_pre:
                        preset = st.selectbox(
                            "Preset:", ["Custom"] + list(SIGNIFICANCE_PRESETS.keys()),
                            key="comp_preset",
                            help="Most stringent: FDR<0.01, |log2FC|>1.5 — Standard: FDR<0.05, |log2FC|>1.0 — Exploratory: FDR<0.10, |log2FC|>0.5",
                        )
                    with c_app:
                        st.write("")
                        if st.button("Apply", key="comp_apply_preset") and preset in SIGNIFICANCE_PRESETS:
                            _p = SIGNIFICANCE_PRESETS[preset]
                            st.session_state["comp_fdr"] = _p["fdr"]
                            st.session_state["comp_fc"] = _p["log2fc"]
                            st.rerun()

                    c1, c2 = st.columns(2)
                    fdr = c1.slider("FDR Cutoff", 0.001, 0.1, 0.05, 0.001, key="comp_fdr")
                    fc = c2.slider("Log2FC Cutoff", 0.0, 5.0, 1.0, 0.1, key="comp_fc")

                    if sel_comp:
                        st.session_state.selected_comparison = sel_comp
                        sig_df = visualizer.filter_significant_proteins(sel_comp, fdr, fc)
                        st.session_state.significant_proteins = sig_df

                        st.markdown("---")
                        st.markdown(f"**Found {len(sig_df)} significant proteins**")

                        up_count = len(sig_df[sig_df['Regulation'] == 'Up-regulated'])
                        down_count = len(sig_df[sig_df['Regulation'] == 'Down-regulated'])
                        m1, m2 = st.columns(2)
                        m1.metric("Up-regulated", f"{up_count}")
                        m2.metric("Down-regulated", f"{down_count}")

                        st.dataframe(sig_df)
                except Exception as _e:
                    logger.error("Comparative Selection tab error", exc_info=True)
                    st.error(f"**Selection tab** failed: {_e}")

            with volcano:
                try:
                    st.markdown("### Volcano Plot")
                    if not st.session_state.get("selected_comparison"): st.warning("Select comparison first.")
                    else:
                        with st.expander("Customize Volcano Plot", expanded=True):
                            auto_top10 = st.toggle("Automatically mark top 10 proteins", value=True)

                            prot_id_col = visualizer.column_config['protein_id']
                            gene_col = 'Gene Name' if 'Gene Name' in visualizer.protein_df.columns else None
                            opts = [f"{g} ({p})" for g, p in visualizer.protein_df[[gene_col, prot_id_col]].dropna().values] if gene_col else visualizer.protein_df[prot_id_col].unique()
                            sel_prots = st.multiselect("Mark specific proteins:", opts)
                            clean_sel = [x.split(" (")[0] for x in sel_prots] if gene_col else sel_prots

                            color_opt = st.selectbox("Highlight method:", ['None', 'Custom List', 'Transcription Factors', 'Keyword Search'], key="volc_col")
                            cust, keyw = None, None
                            if color_opt == 'Custom List': cust = [x.strip() for x in st.text_area("IDs/Genes (one per line):").split('\n') if x.strip()]
                            elif color_opt == 'Keyword Search': keyw = st.text_input("Enter keyword:")

                            st.markdown("---")
                            st.markdown("**Arrow Customization**")
                            ac1, ac2, ac3 = st.columns(3)
                            ah = ac1.slider("Head Size", 1, 8, 2, key="volc_ah")
                            as_ = ac2.slider("Length", 1, 5, 1, key="volc_as")
                            aw = ac3.slider("Width", 1, 5, 1, key="volc_aw")
                            arrow_config = {'arrowhead': ah, 'arrowsize': as_, 'arrowwidth': aw}

                        pm_v = PlotManager("comp_volcano")
                        pm_v.module = "comparative"
                        pm_v.render_generate_button(
                            visualizer.plot_volcano,
                            fdr_cutoff=fdr, fc_cutoff=fc, proteins_to_annotate=clean_sel,
                            color_by_option=color_opt, custom_list=cust, keyword=keyw,
                            annotate_top_10=auto_top10, arrow_config=arrow_config, **global_plot_kwargs
                        )
                        pm_v.render_plot_and_editor()
                except Exception as _e:
                    logger.error("Comparative Volcano tab error", exc_info=True)
                    st.error(f"**Volcano tab** failed: {_e}")

            with heat:
                try:
                    st.markdown("### Heatmap")
                    sel_method = st.radio("Select proteins to display:", ["All Significant", "Top 10 DE", "Up-regulated", "Down-regulated", "Custom Selection"])
                    prot_list = []
                    _sig = st.session_state.get("significant_proteins", pd.DataFrame())

                    if sel_method == "All Significant": prot_list = _sig[pid].tolist()
                    elif sel_method == "Top 10 DE": prot_list = _sig.sort_values(visualizer.column_config['fdr']).head(10)[pid].tolist()
                    elif sel_method == "Up-regulated": prot_list = _sig[_sig['Regulation']=='Up-regulated'][pid].tolist()
                    elif sel_method == "Down-regulated": prot_list = _sig[_sig['Regulation']=='Down-regulated'][pid].tolist()
                    elif sel_method == "Custom Selection":
                        protein_info_df = visualizer.protein_df[[visualizer.column_config['protein_id'], 'Gene Name']].drop_duplicates().reset_index(drop=True)
                        st.markdown("Select proteins from the table below:")
                        selection = st.dataframe(protein_info_df, use_container_width=True, hide_index=True, on_select="rerun", selection_mode="multi-row", key="heat_custom_sel")
                        if selection.selection["rows"]:
                            prot_list = protein_info_df.iloc[selection.selection["rows"]][visualizer.column_config['protein_id']].tolist()

                    if prot_list:
                        pm_heat = MplPlotManager("comp_heatmap")
                        pm_heat.module = "comparative"
                        pm_heat.render_generate_button(visualizer.plot_comparative_heatmap, protein_list=prot_list)
                        pm_heat.render_plot_and_editor()
                    elif sel_method == "Custom Selection":
                        st.info("Please select rows in the table above.")
                    else:
                        st.warning("No significant proteins to display. Run the Selection tab first.")
                except Exception as _e:
                    logger.error("Comparative Heatmap tab error", exc_info=True)
                    st.error(f"**Heatmap tab** failed: {_e}")

            with expr:
                try:
                    st.markdown("### Expression Violin Plots")
                    with st.expander("Configure Violin Plots", expanded=True):
                        sel_method_v = st.radio("Select proteins:", ["Top 10 DE", "Top 10 Up", "Top 10 Down", "Custom"], key="vio_sel")
                        p_list_v = []
                        _sig_e = st.session_state.get("significant_proteins", pd.DataFrame())

                        if sel_method_v == "Top 10 DE":
                            p_list_v = _sig_e.sort_values(visualizer.column_config['fdr']).head(10)[pid].tolist()
                        elif sel_method_v == "Top 10 Up":
                            p_list_v = _sig_e[_sig_e['Regulation']=='Up-regulated'].sort_values(visualizer.column_config['fdr']).head(10)[pid].tolist()
                        elif sel_method_v == "Top 10 Down":
                            p_list_v = _sig_e[_sig_e['Regulation']=='Down-regulated'].sort_values(visualizer.column_config['fdr']).head(10)[pid].tolist()
                        elif sel_method_v == "Custom":
                            df_info = visualizer.protein_df[[visualizer.column_config['protein_id'], 'Gene Name']].dropna()
                            st.markdown("Select proteins from the table below:")
                            sel = st.dataframe(df_info, on_select="rerun", selection_mode="multi-row", key="vio_custom_sel")
                            if sel.selection["rows"]:
                                p_list_v = df_info.iloc[sel.selection["rows"]][visualizer.column_config['protein_id']].tolist()

                    if p_list_v:
                        pm_vio = PlotManager("comp_violin")
                        pm_vio.module = "comparative"
                        pm_vio.render_generate_button(visualizer.plot_expression_violin, protein_list=p_list_v, **global_plot_kwargs)
                        pm_vio.render_plot_and_editor()
                    elif sel_method_v == "Custom":
                        st.info("Select proteins in the table above to generate plot.")
                except Exception as _e:
                    logger.error("Comparative Expression tab error", exc_info=True)
                    st.error(f"**Expression tab** failed: {_e}")

            with path:
                try:
                    st.markdown("### Pathway Enrichment")
                    with st.expander("Configure Analysis", expanded=True):
                        sel_g = st.selectbox("Select gene set:", ["All Significant", "Up-regulated", "Down-regulated", "Custom Selection"])
                        genes = []
                        _sig_p = st.session_state.get("significant_proteins", pd.DataFrame())

                        if sel_g == "All Significant":
                            genes = _sig_p.merge(visualizer.protein_df, left_on=pid, right_on=visualizer.column_config['protein_id'])['Gene Name'].dropna().unique().tolist()
                        elif sel_g == "Up-regulated":
                            genes = _sig_p[_sig_p['Regulation']=='Up-regulated'].merge(visualizer.protein_df, left_on=pid, right_on=visualizer.column_config['protein_id'])['Gene Name'].dropna().unique().tolist()
                        elif sel_g == "Down-regulated":
                            genes = _sig_p[_sig_p['Regulation']=='Down-regulated'].merge(visualizer.protein_df, left_on=pid, right_on=visualizer.column_config['protein_id'])['Gene Name'].dropna().unique().tolist()
                        elif sel_g == "Custom Selection":
                            prot_gene_df = visualizer.protein_df[[visualizer.column_config['protein_id'], 'Gene Name']].dropna().drop_duplicates()
                            st.markdown("Select genes from the table:")
                            gene_sel = st.dataframe(prot_gene_df, on_select="rerun", selection_mode="multi-row", key="path_custom_sel")
                            if gene_sel.selection["rows"]:
                                genes = prot_gene_df.iloc[gene_sel.selection["rows"]]['Gene Name'].tolist()

                        organism = st.selectbox("Select organism:", ["human", "mouse"])

                        bg_mode = st.radio(
                            "Statistical background",
                            ["Detected proteins (recommended)", "Whole genome (Enrichr default)"],
                            index=0,
                            horizontal=True,
                            key="comp_enrich_bg_mode",
                            help="Enrichment p-values depend on the background universe. Testing against "
                                 "the proteins you actually detected — rather than the whole genome — "
                                 "avoids inflating terms made of genes that were never observed in your "
                                 "experiment. Whole genome reproduces the classic Enrichr default.",
                        )

                        # Detected-protein universe = every gene quantified in the uploaded matrix.
                        detected_bg = (
                            visualizer.protein_df['Gene Name'].dropna().astype(str).unique().tolist()
                            if 'Gene Name' in visualizer.protein_df.columns else []
                        )
                        use_detected_bg = bg_mode.startswith("Detected")
                        if use_detected_bg and not detected_bg:
                            st.warning(
                                "No 'Gene Name' column found in the protein matrix — falling back to "
                                "the whole-genome background. Add gene symbols to use a custom background."
                            )
                        if use_detected_bg and detected_bg:
                            st.caption(f"Background = {len(detected_bg)} detected proteins.")

                        if st.button("Run Analysis", type="primary", key="comp_run_pathway_btn"):
                            from utils.caching import run_cached_enrichment
                            if genes:
                                background = detected_bg if (use_detected_bg and detected_bg) else None
                                try:
                                    st.session_state.enrichment_results = run_cached_enrichment(
                                        visualizer, genes, organism, background_genes=background)
                                except RuntimeError as _re:
                                    st.error(f"⚠️ {_re}")
                                except Exception as _ex:
                                    logger.error("Enrichment call failed", exc_info=True)
                                    st.error(f"Enrichment analysis failed unexpectedly: {_ex}")
                            else:
                                st.warning("No genes found/selected.")

                    if 'enrichment_results' in st.session_state and not st.session_state.enrichment_results.empty:
                        res_df = st.session_state.enrichment_results
                        st.markdown("---")

                        q_val = st.slider("Filter by q-value (adj. p):", 0.0, 1.0, 0.2, 0.01, key="comp_path_q")
                        filtered_df = res_df[res_df['adj_p_value'] <= q_val]

                        if filtered_df.empty:
                            st.warning("No terms passed the filter.")
                        else:
                            source_map = {
                                'GO:BP': 'GO Biological Process', 'GO:CC': 'GO Cellular Component',
                                'GO:MF': 'GO Molecular Function', 'KEGG': 'KEGG Pathways', 'REAC': 'Reactome'
                            }
                            sources = sorted(filtered_df['source'].unique().tolist())
                            tab_labels = ["Comprehensive Analysis"] + [source_map.get(s, s) for s in sources]

                            enr_tabs = st.tabs(tab_labels)

                            with enr_tabs[0]:
                                st.markdown("#### Global Manhattan Plot")
                                pm_man = PlotManager("comp_manhattan")
                                pm_man.module = "comparative"
                                pm_man.render_generate_button(visualizer.plot_enrichment_manhattan, enrichment_df=filtered_df, **global_plot_kwargs)
                                pm_man.render_plot_and_editor()

                            for i, source in enumerate(sources):
                                with enr_tabs[i+1]:
                                    st.markdown(f"#### {source_map.get(source, source)}")
                                    source_df = filtered_df[filtered_df['source'] == source].reset_index(drop=True)

                                    with st.expander("Select specific terms to plot", expanded=False):
                                        sel_terms_df = st.dataframe(
                                            source_df[['name', 'adj_p_value', 'intersection_size']],
                                            use_container_width=True, on_select="rerun", selection_mode="multi-row",
                                            key=f"comp_terms_sel_{source}"
                                        )
                                        selected_terms = []
                                        if sel_terms_df.selection["rows"]:
                                            selected_terms = source_df.iloc[sel_terms_df.selection["rows"]]['name'].tolist()

                                    pm_dot = PlotManager(f"comp_dotplot_{source}")
                                    pm_dot.module = "comparative"
                                    pm_dot.render_generate_button(
                                        visualizer.plot_enrichment_dotplot,
                                        enrichment_df=source_df,
                                        terms_to_plot=selected_terms if selected_terms else None,
                                        **global_plot_kwargs
                                    )
                                    pm_dot.render_plot_and_editor()
                except Exception as _e:
                    logger.error("Comparative Pathway tab error", exc_info=True)
                    st.error(f"**Pathway tab** failed: {_e}")

def render(): ComparativeTab().render()
