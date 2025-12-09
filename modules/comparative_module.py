import streamlit as st
import pandas as pd
import logging
import plotly.express as px
from plotly import colors as pc
from visualizations.comparative_visualizer import ComparativeVisualizer
from utils.helpers import handle_plotting_errors
from utils.plot_manager import PlotManager

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

            if st.button("Apply Settings to All Plots", type="primary", use_container_width=True):
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
            
            overview, select, volcano, heat, expr, path = st.tabs(["Overview", "Selection", "Volcano", "Heatmap", "Expression", "Pathway"])

            with overview:
                st.markdown("#### Data Preview")
                st.markdown("**Protein Data**")
                st.dataframe(visualizer.get_protein_data_preview().head())
                st.markdown("**Annotation Data**")
                st.dataframe(visualizer.get_annotation_data_preview().head())
                st.markdown("**Comparative Data**")
                st.dataframe(visualizer.get_comparative_data_preview().head())

            with select:
                comps = visualizer.get_comparison_groups()
                if 'selected_comparison' not in st.session_state: st.session_state.selected_comparison = comps[0] if comps else None
                sel_comp = st.selectbox("Select Comparison:", comps, key='sel_comp_main')
                c1, c2 = st.columns(2)
                fdr = c1.slider("FDR Cutoff", 0.001, 0.1, 0.05, 0.001, key="comp_fdr")
                fc = c2.slider("Log2FC Cutoff", 0.0, 5.0, 1.0, 0.1, key="comp_fc")
                
                if sel_comp:
                    st.session_state.selected_comparison = sel_comp
                    sig_df = visualizer.filter_significant_proteins(sel_comp, fdr, fc)
                    st.session_state.significant_proteins = sig_df
                    
                    st.markdown("---")
                    st.markdown(f"**Found {len(sig_df)} significant proteins**")
                    
                    # RESTORED: Metric Counts
                    up_count = len(sig_df[sig_df['Regulation'] == 'Up-regulated'])
                    down_count = len(sig_df[sig_df['Regulation'] == 'Down-regulated'])
                    m1, m2 = st.columns(2)
                    m1.metric("Up-regulated", f"{up_count}")
                    m2.metric("Down-regulated", f"{down_count}")
                    
                    st.dataframe(sig_df)

            with volcano:
                st.markdown("### Volcano Plot")
                if not st.session_state.selected_comparison: st.warning("Select comparison first.")
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
                    pm_v.render_generate_button(
                        visualizer.plot_volcano,
                        fdr_cutoff=fdr, fc_cutoff=fc, proteins_to_annotate=clean_sel,
                        color_by_option=color_opt, custom_list=cust, keyword=keyw,
                        annotate_top_10=auto_top10, arrow_config=arrow_config, **global_plot_kwargs
                    )
                    pm_v.render_plot_and_editor()

            with heat:
                st.markdown("### Static Heatmap")
                # RESTORED: Custom Selection
                sel_method = st.radio("Select proteins to display:", ["All Significant", "Top 10 DE", "Up-regulated", "Down-regulated", "Custom Selection"])
                prot_list = []
                sig_df = st.session_state.significant_proteins
                pid = visualizer.column_config['comp_protein_id']
                
                if sel_method == "All Significant": prot_list = sig_df[pid].tolist()
                elif sel_method == "Top 10 DE": prot_list = sig_df.sort_values(visualizer.column_config['fdr']).head(10)[pid].tolist()
                elif sel_method == "Up-regulated": prot_list = sig_df[sig_df['Regulation']=='Up-regulated'][pid].tolist()
                elif sel_method == "Down-regulated": prot_list = sig_df[sig_df['Regulation']=='Down-regulated'][pid].tolist()
                elif sel_method == "Custom Selection":
                    # RESTORED: Interactive Dataframe for selection
                    protein_info_df = visualizer.protein_df[[visualizer.column_config['protein_id'], 'Gene Name']].drop_duplicates().reset_index(drop=True)
                    st.markdown("Select proteins from the table below:")
                    selection = st.dataframe(protein_info_df, use_container_width=True, hide_index=True, on_select="rerun", selection_mode="multi-row", key="heat_custom_sel")
                    if selection.selection["rows"]:
                        prot_list = protein_info_df.iloc[selection.selection["rows"]][visualizer.column_config['protein_id']].tolist()

                if prot_list:
                    # Auto-generate without button for static plot
                    try: 
                        with st.spinner("Generating heatmap..."):
                            st.image(visualizer.plot_comparative_heatmap(prot_list))
                    except Exception as e: st.error(f"Failed: {e}")
                elif sel_method == "Custom Selection":
                    st.info("Please select rows in the table above.")
                else: 
                    st.warning("No proteins selected.")

            with expr:
                st.markdown("### Expression Violin Plots")
                with st.expander("Configure Violin Plots", expanded=True):
                    sel_method_v = st.radio("Select proteins:", ["Top 10 DE", "Top 10 Up", "Top 10 Down", "Custom"], key="vio_sel")
                    p_list_v = []
                    
                    if sel_method_v == "Top 10 DE": 
                        p_list_v = sig_df.sort_values(visualizer.column_config['fdr']).head(10)[pid].tolist()
                    elif sel_method_v == "Top 10 Up": 
                        p_list_v = sig_df[sig_df['Regulation']=='Up-regulated'].sort_values(visualizer.column_config['fdr']).head(10)[pid].tolist()
                    elif sel_method_v == "Top 10 Down": 
                        p_list_v = sig_df[sig_df['Regulation']=='Down-regulated'].sort_values(visualizer.column_config['fdr']).head(10)[pid].tolist()
                    elif sel_method_v == "Custom":
                        df_info = visualizer.protein_df[[visualizer.column_config['protein_id'], 'Gene Name']].dropna()
                        st.markdown("Select proteins from the table below:")
                        sel = st.dataframe(df_info, on_select="rerun", selection_mode="multi-row", key="vio_custom_sel")
                        if sel.selection["rows"]: 
                            p_list_v = df_info.iloc[sel.selection["rows"]][visualizer.column_config['protein_id']].tolist()

                if p_list_v:
                    pm_vio = PlotManager("comp_violin")
                    pm_vio.render_generate_button(visualizer.plot_expression_violin, protein_list=p_list_v, **global_plot_kwargs)
                    pm_vio.render_plot_and_editor()
                elif sel_method_v == "Custom":
                    st.info("Select proteins in the table above to generate plot.")

            with path:
                st.markdown("### Pathway Enrichment")
                with st.expander("Configure Analysis", expanded=True):
                    # RESTORED: Custom Selection
                    sel_g = st.selectbox("Select gene set:", ["All Significant", "Up-regulated", "Down-regulated", "Custom Selection"])
                    genes = []
                    
                    if sel_g == "All Significant": 
                        genes = sig_df.merge(visualizer.protein_df, left_on=pid, right_on=visualizer.column_config['protein_id'])['Gene Name'].dropna().unique().tolist()
                    elif sel_g == "Up-regulated":
                        genes = sig_df[sig_df['Regulation']=='Up-regulated'].merge(visualizer.protein_df, left_on=pid, right_on=visualizer.column_config['protein_id'])['Gene Name'].dropna().unique().tolist()
                    elif sel_g == "Down-regulated":
                        genes = sig_df[sig_df['Regulation']=='Down-regulated'].merge(visualizer.protein_df, left_on=pid, right_on=visualizer.column_config['protein_id'])['Gene Name'].dropna().unique().tolist()
                    elif sel_g == "Custom Selection":
                        # RESTORED: Interactive DataFrame for Genes
                        prot_gene_df = visualizer.protein_df[[visualizer.column_config['protein_id'], 'Gene Name']].dropna().drop_duplicates()
                        st.markdown("Select genes from the table:")
                        gene_sel = st.dataframe(prot_gene_df, on_select="rerun", selection_mode="multi-row", key="path_custom_sel")
                        if gene_sel.selection["rows"]:
                            genes = prot_gene_df.iloc[gene_sel.selection["rows"]]['Gene Name'].tolist()

                    # RESTORED: Organism Selector
                    organism = st.selectbox("Select organism:", ["human", "mouse"])

                    if st.button("Run Analysis", type="primary"):
                        from utils.caching import run_cached_enrichment
                        if genes:
                            st.session_state.enrichment_results = run_cached_enrichment(visualizer, genes, organism)
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
                                pm_dot.render_generate_button(
                                    visualizer.plot_enrichment_dotplot, 
                                    enrichment_df=source_df, 
                                    terms_to_plot=selected_terms if selected_terms else None,
                                    **global_plot_kwargs
                                )
                                pm_dot.render_plot_and_editor()

def render(): ComparativeTab().render()