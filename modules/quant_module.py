import streamlit as st
import pandas as pd
import logging
from io import BytesIO

from visualizations.quant_visualizer import QuantificationVisualizer
from utils.helpers import handle_plotting_errors
from utils.plot_manager import EnhancedPlotManager
from utils.color_manager import ModuleColorManager
from config.quant_config import QUANT_PLOT_REGISTRY

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
        
        # --- 1. Initialize Global Color Manager ---
        groups = visualizer.experimental_groups if hasattr(visualizer, 'experimental_groups') else []
        color_manager = ModuleColorManager("quant", groups)
        
        # --- 2. Render Global Settings ---
        # Returns a dict with 'color_map' and 'template'
        global_settings = color_manager.render_global_settings()
        
        st.header("Quantification Analysis Dashboard")

        intersections_tab, completeness_tab, dist_corr_tab, clustering_tab = st.tabs([
            "🔎 Set Intersections", "📊 Data Completeness",
            "📈 Ranking, Distributions & Correlations", "🧬 Clustering"
        ])

        # --- TAB 1: INTERSECTIONS ---
        with intersections_tab:
            st.markdown("### Visualize Protein Set Overlaps")
            if visualizer.annotation_df is None:
                st.warning("An annotation file is required for these plots.", icon="⚠️")
            else:
                venn_inner_tab, upset_inner_tab = st.tabs(["Venn Diagram", "UpSet Plot"])
                
                with venn_inner_tab:
                    st.info("Select 2 to 6 experimental groups to visualize their protein set overlaps.", icon="ℹ️")
                    
                    all_groups = visualizer.experimental_groups
                    default_selection = all_groups[:2] if len(all_groups) >= 2 else all_groups
                    
                    # --- Form for Venn Selection ---
                    with st.form("quant_venn_form"):
                        selected_groups = st.multiselect(
                            "Select groups for Venn Diagram:", 
                            options=all_groups, 
                            default=default_selection
                        )
                        apply_venn = st.form_submit_button("Apply Selection")
                    
                    if 2 <= len(selected_groups) <= 6:
                        venn_manager = EnhancedPlotManager(
                            "quant_venn", 
                            {"title": "Venn Diagram", "default_height": 500}
                        )
                        venn_manager.render(
                            visualizer.plot_venn_diagram,
                            selected_groups=selected_groups,
                            **global_settings
                        )
                    elif len(selected_groups) > 0: 
                        st.error("Please select between 2 and 6 groups.")
                
                with upset_inner_tab:
                    if len(visualizer.experimental_groups) < 3:
                        st.warning("An UpSet plot is only generated for 3 or more experimental groups.", icon="ℹ️")
                    else:
                        st.info("The UpSet plot shows intersections for all groups.")
                        upset_manager = EnhancedPlotManager(
                            "quant_upset", 
                            {"title": "UpSet Plot", "default_height": 600}
                        )
                        upset_manager.render(
                            visualizer.plot_upset,
                            **global_settings
                        )

        # --- TAB 2: COMPLETENESS ---
        with completeness_tab:
            st.markdown("### Protein Counts, Missingness, and Overlap")
            counts_inner_tab, missingness_inner_tab, overlap_inner_tab = st.tabs(["Protein Counts per Sample", "Missingness Analysis", "Overlap & Coverage"])
            
            with counts_inner_tab:
                counts_manager = EnhancedPlotManager(
                    "quant_counts", 
                    QUANT_PLOT_REGISTRY.get("protein_counts", {"title": "Protein Counts"})
                )
                counts_manager.render(
                    visualizer.plot_protein_counts,
                    **global_settings
                )
            
            with missingness_inner_tab:
                heatmap_tab, distribution_tab = st.tabs(["Missing Value Heatmap", "Missing Value Distribution"])
                
                with heatmap_tab:
                    if visualizer.annotation_df is None:
                        st.warning("An annotation file is required for the missingness heatmap.", icon="⚠️")
                    else:
                        heatmap_manager = EnhancedPlotManager(
                            "quant_miss_heat",
                            QUANT_PLOT_REGISTRY.get("missing_heatmap", {"title": "Missing Values"})
                        )
                        heatmap_manager.render(
                            visualizer.plot_missing_values_heatmap,
                            **global_settings
                        )

                with distribution_tab:
                    st.info("This plot helps diagnose the nature of missing data.")
                    dist_manager = EnhancedPlotManager(
                        "quant_miss_dist",
                        QUANT_PLOT_REGISTRY.get("missing_dist", {"title": "Missing Distribution"})
                    )
                    dist_manager.render(
                        visualizer.plot_missing_value_distribution,
                        **global_settings
                    )
            
            with overlap_inner_tab:
                st.markdown("#### Protein Identification Overlap & Coverage")
                col1, col2 = st.columns(2)
                with col1:
                    overlap_manager = EnhancedPlotManager(
                        "quant_overlap",
                        QUANT_PLOT_REGISTRY.get("protein_overlap", {"title": "Protein Overlap"})
                    )
                    overlap_manager.render(
                        visualizer.plot_protein_overlap,
                        **global_settings
                    )
                with col2:
                    coverage_manager = EnhancedPlotManager(
                        "quant_coverage",
                        {"title": "Protein Coverage", "default_height": 500}
                    )
                    coverage_manager.render(
                        visualizer.plot_protein_coverage_chart,
                        **global_settings
                    )

        # --- TAB 3: DISTRIBUTIONS & CORRELATIONS ---
        with dist_corr_tab:
            st.markdown("### Intensity Distributions and Sample Correlations")
            dist_inner_tab, rank_inner_tab, corr_inner_tab = st.tabs(["Intensity Distribution", "Protein Rank Order", "Correlation Matrix"])

            with dist_inner_tab:
                st.markdown("### Work in progress - slow render of the plot")
                # if visualizer.annotation_df is None:
                #     st.warning("An annotation file is required to generate this plot.", icon="⚠️")
                # else:
                #     int_dist_manager = EnhancedPlotManager(
                #         "quant_int_dist",
                #         {"title": "Intensity Distribution", "default_height": 600}
                #     )
                #     int_dist_manager.render(
                #         visualizer.plot_intensity_distribution,
                #         **global_settings
                #     )
            
            with rank_inner_tab:
                st.markdown("### Protein Rank by Average Intensity")
                
                # --- Form for Rank Order Selection ---
                with st.form(key="quant_rank_data_form"):
                    st.markdown("#### 🔍 Data Selection & Highlighting")
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        protein_id_col = visualizer.protein_col
                        gene_name_col = 'Gene Name' if 'Gene Name' in visualizer.protein_df.columns else None
                        if gene_name_col:
                            protein_options = [f"{gene} ({pid})" for gene, pid in visualizer.protein_df[[gene_name_col, protein_id_col]].dropna().values]
                        else:
                            protein_options = visualizer.protein_df[protein_id_col].dropna().unique().tolist()
                        
                        proteins_to_annotate = st.multiselect("Select proteins to annotate:", options=protein_options)
                        cleaned_annotations = [item.split(" (")[0] for item in proteins_to_annotate] if gene_name_col else proteins_to_annotate

                    with col2:
                        color_by_option = st.selectbox("Highlighting Method:", ['None', 'Custom List', 'Transcription Factors', 'Keyword Search'])
                        
                        custom_list_input = ""
                        keyword = ""
                        annotate_highlighted = False
                        
                        if color_by_option == 'Custom List':
                            custom_list_input = st.text_area("Enter IDs/Genes (one per line):")
                            if custom_list_input: annotate_highlighted = st.checkbox("Annotate these?", value=True)
                        
                        elif color_by_option == 'Transcription Factors':
                            annotate_highlighted = st.checkbox("Annotate Top 10 TFs?", value=True)

                        elif color_by_option == 'Keyword Search':
                            keyword = st.text_input("Search Description:")
                            if keyword: annotate_highlighted = st.checkbox("Annotate matches?", value=True)

                    apply_rank = st.form_submit_button("Apply Data Selection")

                # Process inputs outside form
                custom_list = [item.strip() for item in custom_list_input.split('\n') if item.strip()] if custom_list_input else None

                rank_manager = EnhancedPlotManager(
                    "quant_rank_order",
                    QUANT_PLOT_REGISTRY.get("rank_order", {"title": "Rank Order", "has_markers": True})
                )
                rank_manager.render(
                    visualizer.plot_protein_rank_order,
                    proteins_to_annotate=cleaned_annotations,
                    color_by_option=color_by_option,
                    custom_list=custom_list,
                    keyword=keyword,
                    annotate_highlighted=annotate_highlighted,
                    **global_settings
                )

            with corr_inner_tab:
                st.info("This section is under development. 🏗️")
            
        # --- TAB 4: CLUSTERING ---
        with clustering_tab:
            st.markdown("### Sample Clustering Analysis")
            pca_anno_tab, pca_cluster_tab, dendro_tab = st.tabs([
                "PCA by Annotation", "PCA by Cluster", "Hierarchical Clustering"
            ])

            with pca_anno_tab:
                if visualizer.annotation_df is None:
                    st.warning("An annotation file is required.", icon="⚠️")
                else:
                    annotation_cols = [col for col in visualizer.annotation_df.columns if col != visualizer.sample_col]
                    
                    # --- Form for PCA Annotation ---
                    with st.form("quant_pca_anno_form"):
                        col1, col2 = st.columns(2)
                        with col1:
                            default_idx = annotation_cols.index(visualizer.group_col) if visualizer.group_col in annotation_cols else 0
                            color_by = st.selectbox("Color by:", options=annotation_cols, index=default_idx)
                        with col2:
                            symbol_by = st.selectbox("Shape by:", options=[None] + annotation_cols)
                        
                        apply_pca_anno = st.form_submit_button("Update PCA")

                    pca_anno_manager = EnhancedPlotManager(
                        "quant_pca_anno",
                        QUANT_PLOT_REGISTRY.get("pca_anno", {"title": "PCA", "has_markers": True})
                    )
                    pca_anno_manager.render(
                        visualizer.plot_pca_by_annotation,
                        color_by=color_by,
                        symbol_by=symbol_by,
                        **global_settings
                    )

            with pca_cluster_tab:
                # --- Form for K-Means ---
                with st.form("quant_pca_cluster_form"):
                    n_clusters = st.number_input("Number of clusters (k):", 2, 10, 3)
                    apply_pca_cluster = st.form_submit_button("Run Clustering")
                
                pca_cluster_manager = EnhancedPlotManager(
                    "quant_pca_cluster",
                    QUANT_PLOT_REGISTRY.get("pca_cluster", {"title": "Clustered PCA"})
                )
                pca_cluster_manager.render(
                    visualizer.plot_pca_with_clusters,
                    n_clusters=n_clusters,
                    **global_settings
                )

            with dendro_tab:
                # --- Form for Dendrogram ---
                with st.form("quant_dendro_form"):
                    method = st.selectbox("Clustering method:", ['ward', 'complete', 'average', 'single'])
                    apply_dendro = st.form_submit_button("Update Dendrogram")
                
                dendro_manager = EnhancedPlotManager(
                    "quant_dendro",
                    QUANT_PLOT_REGISTRY.get("dendrogram", {"title": "Dendrogram"})
                )
                dendro_manager.render(
                    visualizer.plot_dendrogram,
                    method=method,
                    **global_settings
                )

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