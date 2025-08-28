import streamlit as st
import pandas as pd
import logging
from visualizations.comparative_visualizer import ComparativeVisualizer
from utils.helpers import handle_plotting_errors


logger = logging.getLogger(__name__)

class ComparativeTab:
    """
    Encapsulates the UI and logic for the Comparative Analysis tab.
    """
    def __init__(self):
        if 'comp_visualizer' not in st.session_state:
            st.session_state.comp_visualizer = None

    def _upload_and_config_section(self):
        """Handles file uploads and column name configuration."""
        with st.expander("Upload & Configure Data", expanded=(st.session_state.comp_visualizer is None)):
            
            st.markdown("##### 1. Upload Files")
            protein_file = st.file_uploader("Upload Protein-Level Data", type=['txt', 'csv', 'tsv'], key="comp_protein_uploader")
            annotation_file = st.file_uploader("Upload Annotation File", type=['txt', 'csv', 'tsv'], key="comp_annotation_uploader")
            comparative_file = st.file_uploader("Upload Comparative Analysis File", type=['txt', 'csv', 'tsv'], key="comp_file_uploader")

            st.markdown("---")
            st.markdown("##### 2. Configure Column Names")
            c1, c2 = st.columns(2)
            with c1:
                st.markdown("**Protein & Annotation Files**")
                protein_id_col = st.text_input("Protein ID Column", value="Protein", key="protein_id_main")
                sample_id_col = st.text_input("Sample Linking Column", value="Level3", key="sample_id")
                grouping_col = st.text_input("Experimental Group Column", value="attribute_ExperimentalGroup", key="grouping_id")
            with c2:
                st.markdown("**Comparative Analysis File**")
                comp_protein_id_col = st.text_input("Protein ID Column (in Comp. File)", value="Protein", key="protein_id_comp")
                fold_change_col = st.text_input("Log2 Fold Change Column", value="log2FC", key="fc_id")
                fdr_col = st.text_input("FDR / p-value Column", value="Imputed.FDR", key="fdr_id")
                comparison_label_col = st.text_input("Comparison Label Column", value="Label2", key="label_id")

            if st.button("Load & Process Comparative Data", use_container_width=True, type="primary"):
                if all([protein_file, annotation_file, comparative_file]):
                    try:
                        protein_df = pd.read_csv(protein_file, sep=None, engine='python')
                        annotation_df = pd.read_csv(annotation_file, sep=None, engine='python')
                        comparative_df = pd.read_csv(comparative_file, sep=None, engine='python')

                        column_config = {
                            "protein_id": protein_id_col, "sample_id": sample_id_col,
                            "grouping": grouping_col, "comp_protein_id": comp_protein_id_col,
                            "fold_change": fold_change_col, "fdr": fdr_col,
                            "comparison_label": comparison_label_col
                        }

                        with st.spinner("Initializing comparative visualizer..."):
                            st.session_state.comp_visualizer = ComparativeVisualizer(
                                protein_df, annotation_df, comparative_df, column_config
                            )
                        st.success("All data loaded successfully!")
                        st.rerun()

                    except Exception as e:
                        st.error(f"An error occurred: {e}")
                        logger.error(f"Comparative data loading error: {e}", exc_info=True)
                else:
                    st.warning("Please upload all three required files.")

    def _display_overview(self):
        """Displays previews of the loaded dataframes."""
        visualizer = st.session_state.comp_visualizer
        st.header("Data Overview")

        st.markdown("#### Protein-Level Data")
        st.dataframe(visualizer.get_protein_data_preview().head())

        st.markdown("#### Annotation Data")
        st.dataframe(visualizer.get_annotation_data_preview().head())

        st.markdown("#### Comparative Analysis Data")
        st.dataframe(visualizer.get_comparative_data_preview().head())


    def render(self):
        """Renders the entire Comparative Analysis tab."""
        st.header("Comparative Analysis")
        self._upload_and_config_section()

        if 'significant_proteins' not in st.session_state:
            st.session_state.significant_proteins = pd.DataFrame()

        if st.session_state.comp_visualizer:
            visualizer = st.session_state.comp_visualizer
            
            
            overview_tab, selection_tab, volcano_tab, heatmap_tab, expression_tab, pathway_tab = st.tabs([
                "📥 Data Overview", "🎯 Significant Protein Selection", "🌋 Volcano Plots", 
                "🔥 Heatmaps", "🎻 Expression Plots", "🌐 Pathway Analysis"
            ])

            with overview_tab:
                self._display_overview()

            with selection_tab:
                # ... (Code is unchanged)
                st.markdown("### Filter for Significant Proteins")
                st.info("Select a comparison and set thresholds. These settings will be used as the default for other plots.")
                available_comparisons = visualizer.get_comparison_groups()
                if 'selected_comparison' not in st.session_state:
                    st.session_state.selected_comparison = available_comparisons[0] if available_comparisons else None
                selected_comparison = st.selectbox("Select a comparison to analyze:", options=available_comparisons, key='selected_comparison')
                col1, col2 = st.columns(2)
                with col1:
                    fdr_cutoff = st.slider("FDR (p-value) Cutoff", 0.001, 0.1, 0.05, 0.001, format="%.3f", key="fdr_cutoff")
                with col2:
                    fc_cutoff = st.slider("Absolute Log2 Fold Change Cutoff", 0.0, 5.0, 1.0, 0.1, format="%.1f", key="fc_cutoff")
                if selected_comparison:
                    try:
                        significant_df = visualizer.filter_significant_proteins(selected_comparison, fdr_cutoff, fc_cutoff)
                        st.session_state.significant_proteins = significant_df
                        st.markdown("---"); st.markdown(f"**Found {len(significant_df)} significant proteins**")
                        up_count = len(significant_df[significant_df['Regulation'] == 'Up-regulated'])
                        down_count = len(significant_df[significant_df['Regulation'] == 'Down-regulated'])
                        c1, c2 = st.columns(2); c1.metric("Up-regulated", f"{up_count}"); c2.metric("Down-regulated", f"{down_count}")
                        st.dataframe(significant_df)
                    except Exception as e: st.error(f"An error occurred during filtering: {e}")

            with volcano_tab:
                st.markdown("### Volcano Plot")
                if not st.session_state.selected_comparison:
                    st.warning("Please select a comparison in the 'Significant Protein Selection' tab first.")
                else:
                    st.info(f"Displaying volcano plot for **{st.session_state.selected_comparison}**.")
                    with st.expander("Customize Plot"):
                        st.markdown("**1. Annotate Proteins**")
                        
                        # --- NEW: Add toggle for default annotations ---
                        annotate_top_10 = st.toggle("Automatically mark top 10 proteins", value=True)

                        protein_id_col = visualizer.column_config['protein_id']
                        gene_name_col = 'Gene Name' if 'Gene Name' in visualizer.protein_df.columns else None
                        if gene_name_col: protein_options = [f"{gene} ({pid})" for gene, pid in visualizer.protein_df[[gene_name_col, protein_id_col]].dropna().values]
                        else: protein_options = visualizer.protein_df[protein_id_col].dropna().unique().tolist()
                        proteins_to_annotate = st.multiselect("Also mark these specific proteins:", options=protein_options, key="volcano_annotate")
                        cleaned_annotations = [item.split(" (")[0] for item in proteins_to_annotate] if gene_name_col else proteins_to_annotate
                        
                        st.markdown("**2. Highlight Proteins by Category**")
                        color_by_option = st.selectbox("Choose a highlighting method:", ['None', 'Custom List', 'Transcription Factors', 'Keyword Search'], key="volcano_color")
                        custom_list, keyword = None, None
                        if color_by_option == 'Custom List':
                            custom_list_input = st.text_area("Enter a list of Protein IDs or Gene Names (one per line):", key="volcano_custom_list")
                            custom_list = [item.strip() for item in custom_list_input.split('\n') if item.strip()]
                        elif color_by_option == 'Transcription Factors':
                            tf_count = visualizer.get_transcription_factor_count()
                            if tf_count == 0:
                                st.warning("No proteins in your data matched the built-in list of transcription factors.")
                            else:
                                st.info(f"Found **{tf_count}** potential transcription factors to highlight.")
                        elif color_by_option == 'Keyword Search': keyword = st.text_input("Enter a keyword to search in 'Protein Description':", key="volcano_keyword")

                    try:
                        comp_df_subset = visualizer.comparative_df[visualizer.comparative_df[visualizer.column_config['comparison_label']] == st.session_state.selected_comparison]
                        temp_visualizer = ComparativeVisualizer(
                            visualizer.protein_df, visualizer.annotation_df, comp_df_subset, visualizer.column_config
                        )
                        with st.spinner("Generating volcano plot..."):
                            fig = temp_visualizer.plot_volcano(
                                fdr_cutoff=st.session_state.fdr_cutoff, fc_cutoff=st.session_state.fc_cutoff,
                                proteins_to_annotate=cleaned_annotations,
                                color_by_option=color_by_option,
                                custom_list=custom_list,
                                keyword=keyword,
                                annotate_top_10=annotate_top_10 # Pass toggle state
                            )
                            st.plotly_chart(fig, use_container_width=True)
                    except Exception as e:
                        st.error(f"Failed to generate volcano plot: {e}")

            with heatmap_tab:
                # ... (Code is unchanged)
                st.markdown("### Clustered Heatmap of Protein Abundance")
                protein_selection_method = st.radio(
                    "Select which proteins to display in the heatmap:",
                    options=["All Significant", "Top 10 Differentially Expressed", "Up-regulated Only", "Down-regulated Only", "Custom Selection"],
                    horizontal=True, key="heatmap_protein_select"
                )
                protein_id_col = visualizer.column_config['comp_protein_id']; fdr_col = visualizer.column_config['fdr']; protein_list = []
                significant_df = st.session_state.significant_proteins
                if protein_selection_method == "All Significant":
                    protein_list = significant_df[protein_id_col].tolist()
                    st.markdown(f"Displaying **{len(protein_list)}** significant proteins.")
                elif protein_selection_method == "Top 10 Differentially Expressed":
                    top_10_df = significant_df.sort_values(by=fdr_col).head(10)
                    protein_list = top_10_df[protein_id_col].tolist()
                    st.markdown(f"Displaying the **Top {len(protein_list)}** most significant proteins.")
                elif protein_selection_method == "Up-regulated Only":
                    protein_list = significant_df[significant_df['Regulation'] == 'Up-regulated'][protein_id_col].tolist()
                    st.markdown(f"Displaying **{len(protein_list)}** up-regulated significant proteins.")
                elif protein_selection_method == "Down-regulated Only":
                    protein_list = significant_df[significant_df['Regulation'] == 'Down-regulated'][protein_id_col].tolist()
                    st.markdown(f"Displaying **{len(protein_list)}** down-regulated significant proteins.")
                elif protein_selection_method == "Custom Selection":
                    protein_info_df = visualizer.protein_df[[visualizer.column_config['protein_id'], 'Gene Name', 'Protein Description']].drop_duplicates().reset_index(drop=True)
                    selected_proteins_df = st.dataframe(protein_info_df, use_container_width=True, hide_index=True, on_select="rerun", selection_mode="multi-row")
                    selected_indices = selected_proteins_df.selection["rows"]
                    if selected_indices: protein_list = protein_info_df.iloc[selected_indices][visualizer.column_config['protein_id']].tolist()
                    st.markdown(f"**{len(protein_list)}** proteins selected.")
                if not protein_list: st.warning("No proteins selected to display. Please adjust your filters or selection.")
                else:
                    try:
                        with st.spinner("Generating heatmap..."):
                            plot_buffer = visualizer.plot_comparative_heatmap(protein_list)
                            st.image(plot_buffer)
                    except Exception as e: st.error(f"Failed to generate heatmap: {e}")

            with expression_tab:
                st.markdown("### Protein Expression Violin Plots")
                st.info("Visualize the distribution of protein abundance for specific protein sets across experimental groups.")

                protein_id_col = visualizer.column_config['comp_protein_id']
                fdr_col = visualizer.column_config['fdr']
                fc_col = visualizer.column_config['fold_change']
                protein_list = []
                significant_df = st.session_state.significant_proteins

                protein_selection_method = st.radio(
                    "Select which proteins to display:",
                    options=["Top 10 Differentially Expressed", "Top 10 Up-regulated", "Top 10 Down-regulated", "Custom Selection"],
                    horizontal=True, key="violin_protein_select",
                    help="DE = Differentially Expressed (lowest FDR)"
                )
                
                if protein_selection_method == "Top 10 Differentially Expressed":
                    protein_list = significant_df.sort_values(by=fdr_col).head(10)[protein_id_col].tolist()
                
                elif protein_selection_method == "Top 10 Up-regulated":
                    up_df = significant_df[significant_df['Regulation'] == 'Up-regulated']
                    protein_list = up_df.sort_values(by=fdr_col).head(10)[protein_id_col].tolist()

                elif protein_selection_method == "Top 10 Down-regulated":
                    down_df = significant_df[significant_df['Regulation'] == 'Down-regulated']
                    protein_list = down_df.sort_values(by=fdr_col).head(10)[protein_id_col].tolist()

                elif protein_selection_method == "Custom Selection":
                    protein_info_df = visualizer.protein_df[[visualizer.column_config['protein_id'], 'Gene Name']].drop_duplicates()
                    protein_options = [f"{gene} ({pid})" for gene, pid in protein_info_df.values]
                    selected_options = st.multiselect("Select proteins to plot:", options=protein_options)
                    protein_list = [opt.split('(')[-1].strip(')') for opt in selected_options]
                
                if not protein_list:
                    st.warning("No proteins selected to display. Please make a selection.")
                else:
                    try:
                        with st.spinner(f"Generating violin plots for {len(protein_list)} proteins..."):
                            fig = visualizer.plot_expression_violin(protein_list)
                            st.plotly_chart(fig, use_container_width=True)
                    except Exception as e:
                        st.error(f"Failed to generate violin plot: {e}")

            with pathway_tab:
                st.markdown("### Pathway Enrichment Analysis")
                st.info("Perform pathway enrichment on a selected gene set using the Enrichr API.")

                gene_list_option = st.selectbox(
                    "Select gene set for analysis:",
                    options=["All Significant", "Up-regulated Only", "Down-regulated Only", "Custom Selection"]
                )
                
                gene_list = []
                significant_df = st.session_state.significant_proteins
                protein_info = visualizer.protein_df[[visualizer.column_config['protein_id'], 'Gene Name']].dropna()
                significant_genes_df = pd.merge(significant_df, protein_info, left_on=visualizer.column_config['comp_protein_id'], right_on=visualizer.column_config['protein_id'])

                if gene_list_option == "All Significant":
                    gene_list = significant_genes_df['Gene Name'].unique().tolist()
                elif gene_list_option == "Up-regulated Only":
                    gene_list = significant_genes_df[significant_genes_df['Regulation'] == 'Up-regulated']['Gene Name'].unique().tolist()
                elif gene_list_option == "Down-regulated Only":
                    gene_list = significant_genes_df[significant_genes_df['Regulation'] == 'Down-regulated']['Gene Name'].unique().tolist()
                
                elif gene_list_option == "Custom Selection":
                    # --- NEW: Interactive dataframe for custom selection ---
                    st.markdown("Select specific proteins from the table below to include in the analysis.")
                    protein_info_df = visualizer.protein_df[[
                        visualizer.column_config['protein_id'], 'Gene Name', 'Protein Description'
                    ]].drop_duplicates().reset_index(drop=True)

                    selection = st.dataframe(
                        protein_info_df,
                        use_container_width=True, hide_index=True, on_select="rerun",
                        selection_mode="multi-row", key="pathway_custom_select"
                    )
                    
                    selected_indices = selection.selection["rows"]
                    if selected_indices:
                        gene_list = protein_info_df.iloc[selected_indices]['Gene Name'].dropna().unique().tolist()

                st.markdown(f"**{len(gene_list)}** unique genes selected for analysis.")
                organism = st.selectbox("Select organism:", options=["human", "mouse"])

                if st.button("Run Pathway Analysis", use_container_width=True) and gene_list:
                    from utils.caching import run_cached_enrichment
                    enrichment_results = run_cached_enrichment(visualizer, gene_list, organism)
                    st.session_state.enrichment_results = enrichment_results
                    st.rerun()
                
                if 'enrichment_results' in st.session_state and not st.session_state.enrichment_results.empty:
                    st.success("Enrichment analysis complete! View results in the tabs below.")
                    
                    q_value_cutoff = st.slider(
                        "Filter results by q-value (adjusted p-value):",
                        min_value=0.0, max_value=1.0, value=0.2, step=0.01
                    )
                    
                    results_df = st.session_state.enrichment_results
                    filtered_results_df = results_df[results_df['adj_p_value'] <= q_value_cutoff]
                    
                    if filtered_results_df.empty:
                        st.warning("No enrichment terms passed the current q-value filter.")
                        return

                    source_name_map = {
                        'GO:BP': 'GO Biological Process', 'GO:CC': 'GO Cellular Component',
                        'GO:MF': 'GO Molecular Function', 'KEGG': 'KEGG Pathways', 'REAC': 'Reactome Pathways'
                    }
                    available_sources = sorted(filtered_results_df['source'].unique().tolist())
                    tab_names = ["Comprehensive Analysis"] + [source_name_map.get(s, s) for s in available_sources]
                    
                    tabs = st.tabs(tab_names)
                    
                    with tabs[0]: # Comprehensive Manhattan Plot
                        st.markdown("#### Comprehensive Manhattan Plot")
                        fig = visualizer.plot_enrichment_manhattan(filtered_results_df)
                        st.plotly_chart(fig, use_container_width=True)

                    for i, source_code in enumerate(available_sources):
                        with tabs[i+1]:
                            full_source_name = source_name_map.get(source_code, source_code)
                            st.markdown(f"#### Enriched Terms for {full_source_name}")
                            source_df = filtered_results_df[filtered_results_df['source'] == source_code].reset_index(drop=True)

                            with st.expander("Select custom terms to plot"):
                                selection = st.dataframe(
                                    source_df[['name', 'adj_p_value', 'genes', 'intersection_size']],
                                    use_container_width=True, hide_index=True, on_select="rerun",
                                    selection_mode="multi-row", key=f"table_{source_code}"
                                )
                            
                            selected_indices = selection.selection["rows"]
                            if selected_indices:
                                selected_terms = source_df.iloc[selected_indices]['name'].tolist()
                                st.markdown(f"**Displaying {len(selected_terms)} selected terms.**")
                                fig = visualizer.plot_enrichment_dotplot(source_df, terms_to_plot=selected_terms)
                            else:
                                st.markdown("**Displaying top 10 terms by q-value.**")
                                fig = visualizer.plot_enrichment_dotplot(source_df)
                            
                            st.plotly_chart(fig, use_container_width=True)

        else:
            st.info("Upload your data to begin the comparative analysis.")

def render():
    """Entry point function to render the ComparativeTab."""
    comp_tab = ComparativeTab()
    comp_tab.render()