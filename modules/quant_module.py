import streamlit as st
import pandas as pd
import logging
from visualizations.quant_visualizer import QuantificationVisualizer
from utils.helpers import handle_plotting_errors
from io import BytesIO

# Set up a logger for this module
logger = logging.getLogger(__name__)

class QuantificationTab:
    """
    Encapsulates the Streamlit UI and logic for the Quantification Analysis tab.
    """

    def __init__(self):
        """Initializes the QuantificationTab and its session state."""
        if 'quant_visualizer' not in st.session_state:
            st.session_state.quant_visualizer = None

    def _upload_data_section(self):
        """Handles the file uploader and data processing logic."""
        with st.expander("Upload & Configure Data", expanded=(st.session_state.quant_visualizer is None)):
            
            st.markdown("##### 1. Upload Your Files")
            protein_file = st.file_uploader(
                "Upload Protein-Level Data",
                type=['txt', 'csv', 'tsv'],
                key="quant_protein_uploader"
            )
            
            annotation_file = st.file_uploader(
                "Upload Annotation File (Optional)",
                type=['txt', 'csv', 'tsv'],
                key="quant_annotation_uploader"
            )
            
            st.markdown("---")
            st.markdown("##### 2. Configure Column Names")
            col1, col2, col3 = st.columns(3)
            with col1:
                protein_col = st.text_input("Protein ID Column", value="Protein")
            with col2:
                sample_col = st.text_input("Sample Linking Column (in Annotation)", value="Level3")
            with col3:
                group_col = st.text_input("Experimental Group Column (in Annotation)", value="attribute_ExperimentalGroup")


            if st.button("Load & Process Data", use_container_width=True, type="primary"):
                if protein_file:
                    try:
                        protein_df = pd.read_csv(protein_file, sep=None, engine='python')
                        annotation_df = pd.read_csv(annotation_file, sep=None, engine='python') if annotation_file else None
                        
                        with st.spinner("Initializing visualizer..."):
                            st.session_state.quant_visualizer = QuantificationVisualizer(
                                protein_df, annotation_df,
                                protein_col=protein_col,
                                sample_col=sample_col,
                                group_col=group_col
                            )
                        
                        st.success("Data loaded successfully!")
                        st.rerun()

                    except Exception as e:
                        st.error(f"An error occurred: {e}")
                        logger.error(f"Quantification data loading error: {e}", exc_info=True)
                        st.session_state.quant_visualizer = None
                else:
                    st.warning("Please upload at least a protein-level data file.")

    @handle_plotting_errors
    def _display_results(self):
        """Displays the main analysis dashboard with different tabs."""
        visualizer = st.session_state.quant_visualizer
        st.header("Quantification Analysis Dashboard")

        intersections_tab, completeness_tab, dist_corr_tab, clustering_tab = st.tabs([
            "🔎 Set Intersections", "📊 Data Completeness",
            "📈 Ranking, Distributions & Correlations", "🧬 Clustering"
        ])

        with intersections_tab:
            # ... (Code is unchanged)
            st.markdown("### Visualize Protein Set Overlaps")
            if visualizer.annotation_df is None:
                st.warning("An annotation file is required for these plots.", icon="⚠️")
            else:
                venn_inner_tab, upset_inner_tab = st.tabs(["Venn Diagram", "UpSet Plot"])
                with venn_inner_tab:
                    st.info("Select 2 to 6 experimental groups to visualize their protein set overlaps.", icon="ℹ️")
                    all_groups = visualizer.experimental_groups
                    default_selection = all_groups[:2] if len(all_groups) >= 2 else all_groups
                    selected_groups = st.multiselect("Select groups for Venn Diagram:", options=all_groups, default=default_selection, key="venn_group_select")
                    if 2 <= len(selected_groups) <= 6:
                        col1, col2, col3 = st.columns([0.2, 0.6, 0.2])
                        with col2:
                            try:
                                with st.spinner("Creating Venn Diagram..."):
                                    fig = visualizer.plot_venn_diagram(selected_groups)
                                    st.pyplot(fig, use_container_width=True)
                            except Exception as e: st.error(f"Failed to create Venn diagram: {e}")
                    elif len(selected_groups) > 0: st.error("Please select between 2 and 6 groups.")
                with upset_inner_tab:
                    if len(visualizer.experimental_groups) < 3:
                        st.warning("An UpSet plot is only generated for 3 or more experimental groups.", icon="ℹ️")
                    else:
                        st.info("The UpSet plot shows intersections for all groups, ideal for complex comparisons.")
                        try:
                            with st.spinner("Creating UpSet Plot..."):
                                col1, col2, col3 = st.columns([0.1, 0.8, 0.1])
                                with col2:
                                    fig = visualizer.plot_upset()
                                    st.pyplot(fig, use_container_width=True)
                        except Exception as e: st.error(f"Failed to create UpSet plot: {e}")

        with completeness_tab:
            # ... (Code is unchanged)
            st.markdown("### Protein Counts, Missingness, and Overlap")
            counts_inner_tab, missingness_inner_tab, overlap_inner_tab = st.tabs(["Protein Counts per Sample", "Missingness Analysis", "Overlap & Coverage"])
            with counts_inner_tab:
                try:
                    with st.spinner("Generating plot..."):
                        fig = visualizer.plot_protein_counts()
                        st.plotly_chart(fig, use_container_width=True)
                except Exception as e: st.error(f"Failed to generate protein counts plot: {e}")
            with missingness_inner_tab:
                heatmap_tab, distribution_tab = st.tabs(["Missing Value Heatmap", "Missing Value Distribution"])
                with heatmap_tab:
                    if visualizer.annotation_df is None:
                        st.warning("An annotation file is required for the missingness heatmap.", icon="⚠️")
                    else:
                        try:
                            with st.spinner("Generating missingness heatmap..."):
                                fig = visualizer.plot_missing_values_heatmap()
                                st.plotly_chart(fig, use_container_width=True)
                        except Exception as e: st.error(f"Failed to generate missingness plot: {e}")
                with distribution_tab:
                    st.info("This plot helps diagnose the nature of missing data.")
                    try:
                        with st.spinner("Generating distribution plots..."):
                            fig = visualizer.plot_missing_value_distribution()
                            st.plotly_chart(fig, use_container_width=True)
                    except Exception as e: st.error(f"Failed to generate distribution plots: {e}")
            with overlap_inner_tab:
                st.markdown("#### Protein Identification Overlap & Coverage")
                col1, col2 = st.columns(2)
                with col1:
                    try:
                        with st.spinner("Generating overlap plot..."):
                            fig_overlap = visualizer.plot_protein_overlap()
                            st.plotly_chart(fig_overlap, use_container_width=True)
                    except Exception as e: st.error(f"Failed to generate overlap plot: {e}")
                

        with dist_corr_tab:
            st.markdown("### Intensity Distributions and Sample Correlations")
            dist_inner_tab, rank_inner_tab, corr_inner_tab = st.tabs(["Intensity Distribution", "Protein Rank Order", "Correlation Matrix"])

            with dist_inner_tab:
                # ... (Code is unchanged)
                if visualizer.annotation_df is None:
                    st.warning("An annotation file is required to generate this plot.", icon="⚠️")
                else:
                    try:
                        with st.spinner("Generating intensity distribution plot..."):
                            plot_buffer = visualizer.plot_intensity_distribution()
                            st.image(plot_buffer, caption="Group-wise Intensity Distribution")
                    except Exception as e: st.error(f"Failed to generate plot: {e}")
            
            with rank_inner_tab:
                st.markdown("### Protein Rank by Average Intensity")
                st.info("This plot ranks proteins by their mean intensity. Use the options below to highlight and color proteins of interest.")

                with st.expander("Customize Plot"):
                    st.markdown("**1. Annotate Proteins**")
                    protein_id_col = visualizer.protein_col
                    gene_name_col = 'Gene Name' if 'Gene Name' in visualizer.protein_df.columns else None
                    if gene_name_col:
                        protein_options = [f"{gene} ({pid})" for gene, pid in visualizer.protein_df[[gene_name_col, protein_id_col]].dropna().values]
                    else:
                        protein_options = visualizer.protein_df[protein_id_col].dropna().unique().tolist()
                    proteins_to_annotate = st.multiselect("Select specific proteins to mark with an arrow:", options=protein_options)
                    cleaned_annotations = [item.split(" (")[0] for item in proteins_to_annotate] if gene_name_col else proteins_to_annotate

                    st.markdown("**2. Color Proteins by Category**")
                    color_by_option = st.selectbox("Choose a coloring method:", ['None', 'Custom List', 'Transcription Factors', 'Keyword Search'])

                    custom_list, keyword = None, None
                    annotate_highlighted = False # Default to False

                    if color_by_option == 'Custom List':
                        custom_list_input = st.text_area("Enter a list of Protein IDs or Gene Names (one per line):")
                        custom_list = [item.strip() for item in custom_list_input.split('\n') if item.strip()]
                        if custom_list:
                            annotate_highlighted = st.checkbox("Add arrows for highlighted proteins (Top 10)")
                    
                    elif color_by_option == 'Transcription Factors':
                        tf_count = visualizer.get_transcription_factor_count()
                        if tf_count == 0:
                            st.warning("No proteins in your data matched the built-in list of transcription factors.")
                        else:
                            st.info(f"Found **{tf_count}** potential transcription factors to highlight.")
                            annotate_highlighted = st.checkbox("Add arrows for highlighted TFs (Top 10)")

                    elif color_by_option == 'Keyword Search':
                        st.markdown("###### Preview Protein Descriptions")
                        st.dataframe(
                            visualizer.protein_df[[visualizer.protein_col, 'Gene Name', 'Protein Description']],
                            use_container_width=True, height=200
                        )
                        keyword = st.text_input("Enter a keyword to search in 'Protein Description':")
                        if keyword:
                            annotate_highlighted = st.checkbox("Add arrows for highlighted proteins (Top 10)")

                try:
                    with st.spinner("Generating rank order plot..."):
                        fig = visualizer.plot_protein_rank_order(
                            proteins_to_annotate=cleaned_annotations,
                            color_by_option=color_by_option,
                            custom_list=custom_list,
                            keyword=keyword,
                            annotate_highlighted=annotate_highlighted
                        )
                        st.plotly_chart(fig, use_container_width=True)
                except Exception as e:
                    st.error(f"Failed to generate rank order plot: {e}")

            with corr_inner_tab:
                st.info("This section is under development. 🏗️")
            
        with clustering_tab:
            st.markdown("### Sample Clustering Analysis")
            pca_anno_tab, pca_cluster_tab, dendro_tab = st.tabs([
                "PCA by Annotation", "PCA by Cluster", "Hierarchical Clustering"
            ])

            with pca_anno_tab:
                st.markdown("#### Principal Component Analysis (PCA) by Annotation")
                if visualizer.annotation_df is None:
                    st.warning("An annotation file is required to generate this plot.", icon="⚠️")
                else:
                    st.info("Color and shape the PCA plot using columns from your annotation file.")
                    
                    annotation_cols = [col for col in visualizer.annotation_df.columns if col != visualizer.sample_col]
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        color_by = st.selectbox(
                            "Color points by:", options=annotation_cols,
                            index=annotation_cols.index(visualizer.group_col) if visualizer.group_col in annotation_cols else 0
                        )
                    with col2:
                        symbol_by = st.selectbox(
                            "Shape points by:", options=[None] + annotation_cols
                        )
                    
                    # --- NEW: Add toggle for labels ---
                    show_labels_anno = st.toggle("Show sample labels", key="pca_anno_labels")
                    
                    try:
                        with st.spinner("Generating PCA plot..."):
                            fig = visualizer.plot_pca_by_annotation(
                                color_by=color_by, 
                                symbol_by=symbol_by,
                                show_labels=show_labels_anno # Pass the toggle state
                            )
                            st.plotly_chart(fig, use_container_width=True)
                    except Exception as e:
                        st.error(f"Failed to generate PCA plot: {e}")

            with pca_cluster_tab:
                st.markdown("#### PCA with Unsupervised Clustering")
                st.info("Perform k-means clustering on the samples and visualize the results in a PCA plot.")
                
                n_clusters = st.number_input(
                    "Number of clusters (k):", min_value=2, max_value=10, value=3, step=1
                )
                
                # --- NEW: Add toggle for labels ---
                show_labels_cluster = st.toggle("Show sample labels", key="pca_cluster_labels")

                try:
                    with st.spinner("Generating clustered PCA plot..."):
                        fig = visualizer.plot_pca_with_clusters(
                            n_clusters=n_clusters,
                            show_labels=show_labels_cluster # Pass the toggle state
                        )
                        st.plotly_chart(fig, use_container_width=True)
                except Exception as e:
                    st.error(f"Failed to generate plot: {e}")

            with dendro_tab:
                st.markdown("#### Hierarchical Clustering Dendrogram")
                st.info("Visualize how samples cluster together based on their protein expression profiles.")
                
                method = st.selectbox(
                    "Clustering method:",
                    options=['ward', 'complete', 'average', 'single']
                )
                
                try:
                    with st.spinner("Generating dendrogram..."):
                        plot_buffer = visualizer.plot_dendrogram(method=method)
                        st.image(plot_buffer, caption="Sample Clustering Dendrogram")
                except Exception as e:
                    st.error(f"Failed to generate dendrogram: {e}")

    def render(self):
        """Renders the entire Quantification Analysis tab."""
        st.header("Protein Quantification Analysis")
        self._upload_data_section()

        if st.session_state.quant_visualizer:
            self._display_results()

def render():
    """Entry point function to render the QuantificationTab."""
    quant_tab = QuantificationTab()
    quant_tab.render()