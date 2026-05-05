import streamlit as st
import pandas as pd
import logging
import plotly.express as px
from visualizations.quant_visualizer import QuantificationVisualizer
from utils.helpers import handle_plotting_errors, to_hex
from utils.plot_manager import PlotManager, MplPlotManager
from utils.sanity import gene_resolution_report, render_validation

# Set up a logger for this module
logger = logging.getLogger(__name__)

class QuantificationTab:
    """
    Encapsulates the Streamlit UI and logic for the Quantification Analysis tab.
    """
    
    # Keys used for PlotManager state tracking
    PLOT_KEYS = [
        "quant_counts", "quant_missing_heat", "quant_missing_dist",
        "quant_overlap", "quant_coverage", "quant_rank",
        "quant_pca_anno", "quant_pca_cluster",
        "quant_venn", "quant_upset", "quant_dendro",
        "quant_intensity_dist", "quant_corr", "quant_cv_intensity",
    ]

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


            if st.button("Load & Process Data", use_container_width=True, type="primary", key="quant_load_data_btn"):
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

    def _render_global_settings(self):
        """Renders global visual settings and returns kwargs for plot generation."""
        with st.expander("🎨 Global Plot Settings (Colors & Theme)", expanded=False):
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
                    "Plot Theme", 
                    list(theme_options.keys()), 
                    key="quant_global_theme_select"
                )
                selected_template = theme_options[selected_theme_label]

            with c2:
                visualizer = st.session_state.quant_visualizer
                groups = visualizer.experimental_groups if visualizer.experimental_groups else []
                color_map = {}
                
                if groups:
                    st.markdown("**Customize Group Colors:**")
                    cols = st.columns(4)

                    safe_colors = []
                    for c in px.colors.qualitative.Safe:
                        try:
                            safe_colors.append(to_hex(c))
                        except Exception:
                            safe_colors.append(c)

                    for i, group in enumerate(groups):
                        default_c = safe_colors[i % len(safe_colors)]
                        with cols[i % 4]: 
                            color_map[group] = st.color_picker(
                                f"{group}", 
                                value=default_c, 
                                key=f"quant_color_{group}"
                            )
                else:
                    st.caption("No experimental groups found to colorize.")

            st.markdown("---")
            
            if st.button("Apply Settings to All Plots", type="primary", use_container_width=True):
                for key in self.PLOT_KEYS:
                    for suffix in ("_fig", "_mpl_buf"):
                        sk = f"{key}{suffix}"
                        if sk in st.session_state:
                            st.session_state[sk] = None

                st.success("Settings applied! Plots will regenerate automatically.")
                st.rerun()

        global_kwargs = {"template": selected_template}
        if color_map:
            global_kwargs["color_discrete_map"] = color_map
            
        return global_kwargs

    @handle_plotting_errors
    def _display_results(self):
        """Displays the main analysis dashboard with different tabs."""
        visualizer = st.session_state.quant_visualizer
        st.header("Quantification Analysis Dashboard")
        
        # Render global settings and get kwargs
        global_plot_kwargs = self._render_global_settings()

        render_validation(
            [gene_resolution_report(visualizer.protein_df, 'Gene Name')],
            header="Gene Symbol Coverage",
        )

        intersections_tab, completeness_tab, dist_corr_tab, clustering_tab = st.tabs([
            "🔎 Set Intersections", "📊 Data Completeness",
            "📈 Ranking, Distributions & Correlations", "🧬 Clustering"
        ])

        # --- TAB 1: INTERSECTIONS (Static Plots preserved) ---
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
                    selected_groups = st.multiselect("Select groups for Venn Diagram:", options=all_groups, default=default_selection, key="venn_group_select")
                    if 2 <= len(selected_groups) <= 6:
                        mpl_venn = MplPlotManager("quant_venn")
                        mpl_venn.module = "quant"
                        mpl_venn.render_generate_button(visualizer.plot_venn_diagram, selected_groups=selected_groups)
                        mpl_venn.render_plot_and_editor()
                    elif len(selected_groups) > 0:
                        st.error("Please select between 2 and 6 groups.")
                with upset_inner_tab:
                    if len(visualizer.experimental_groups) < 3:
                        st.warning("An UpSet plot is only generated for 3 or more experimental groups.", icon="ℹ️")
                    else:
                        st.info("The UpSet plot shows intersections for all groups, ideal for complex comparisons.")
                        mpl_upset = MplPlotManager("quant_upset")
                        mpl_upset.module = "quant"
                        mpl_upset.render_generate_button(visualizer.plot_upset)
                        mpl_upset.render_plot_and_editor()

        # --- TAB 2: COMPLETENESS (Updated with PlotManager) ---
        with completeness_tab:
            st.markdown("### Protein Counts, Missingness, and Overlap")
            counts_inner_tab, missingness_inner_tab, overlap_inner_tab = st.tabs(["Protein Counts per Sample", "Missingness Analysis", "Overlap & Coverage"])
            
            with counts_inner_tab:
                pm_counts = PlotManager("quant_counts")
                pm_counts.module = "quant"
                pm_counts.render_generate_button(visualizer.plot_protein_counts, **global_plot_kwargs)
                pm_counts.render_plot_and_editor()
            
            with missingness_inner_tab:
                heatmap_tab, distribution_tab = st.tabs(["Missing Value Heatmap", "Missing Value Distribution"])
                with heatmap_tab:
                    if visualizer.annotation_df is None:
                        st.warning("An annotation file is required for the missingness heatmap.", icon="⚠️")
                    else:
                        pm_heat = PlotManager("quant_missing_heat")
                        pm_heat.module = "quant"
                        pm_heat.render_generate_button(visualizer.plot_missing_values_heatmap, **global_plot_kwargs)
                        pm_heat.render_plot_and_editor()

                with distribution_tab:
                    st.info("This plot helps diagnose the nature of missing data.")
                    pm_dist = PlotManager("quant_missing_dist")
                    pm_dist.module = "quant"
                    pm_dist.render_generate_button(visualizer.plot_missing_value_distribution, **global_plot_kwargs)
                    pm_dist.render_plot_and_editor()
            
            with overlap_inner_tab:
                st.markdown("#### Protein Identification Overlap & Coverage")
                col1, col2 = st.columns(2)
                with col1:
                    pm_overlap = PlotManager("quant_overlap")
                    pm_overlap.module = "quant"
                    pm_overlap.render_generate_button(visualizer.plot_protein_overlap, **global_plot_kwargs)
                    pm_overlap.render_plot_and_editor()
                with col2:
                    pm_cov = PlotManager("quant_coverage")
                    pm_cov.module = "quant"
                    pm_cov.render_generate_button(visualizer.plot_protein_coverage_chart, **global_plot_kwargs)
                    pm_cov.render_plot_and_editor()

        # --- TAB 3: DISTRIBUTIONS ---
        with dist_corr_tab:
            st.markdown("### Intensity Distributions and Sample Correlations")
            dist_inner_tab, rank_inner_tab, cv_inner_tab, corr_inner_tab = st.tabs([
                "Intensity Distribution", "Protein Rank Order", "CV vs Intensity", "Correlation Matrix"
            ])

            with dist_inner_tab:
                if visualizer.annotation_df is None:
                    st.warning("An annotation file is required to generate this plot.", icon="⚠️")
                else:
                    mpl_int_dist = MplPlotManager("quant_intensity_dist")
                    mpl_int_dist.module = "quant"
                    mpl_int_dist.render_generate_button(visualizer.plot_intensity_distribution)
                    mpl_int_dist.render_plot_and_editor()
            
            with rank_inner_tab:
                st.markdown("### Protein Rank by Average Intensity")
                st.info("This plot ranks proteins by their mean intensity. Use the options below to highlight and color proteins of interest.")

                # --- CUSTOMIZE PLOT EXPANDER ---
                # Changing widgets here does NOT trigger a plot update until the button is clicked
                with st.expander("Customize Plot", expanded=True):
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
                    annotate_highlighted = False 

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

                    st.markdown("---")
                    st.markdown("**3. Arrow Customization**")
                    ac1, ac2, ac3 = st.columns(3)
                    ah = ac1.slider("Head Size", 1, 8, 2, key="quant_arrow_head")
                    as_ = ac2.slider("Length", 1, 5, 1, key="quant_arrow_size")
                    aw = ac3.slider("Width", 1, 5, 1, key="quant_arrow_width")
                    arrow_config = {'arrowhead': ah, 'arrowsize': as_, 'arrowwidth': aw}

                pm_rank = PlotManager("quant_rank")
                pm_rank.module = "quant"
                pm_rank.render_generate_button(
                    visualizer.plot_protein_rank_order,
                    proteins_to_annotate=cleaned_annotations,
                    color_by_option=color_by_option,
                    custom_list=custom_list,
                    keyword=keyword,
                    annotate_highlighted=annotate_highlighted,
                    arrow_config=arrow_config,
                    **global_plot_kwargs
                )
                pm_rank.render_plot_and_editor()

            with cv_inner_tab:
                st.markdown("### Protein CV (%) vs Mean Intensity")
                st.info("High-CV proteins at low intensity are typically noise. Use this plot to choose an intensity filter before downstream analysis.")
                pm_cv = PlotManager("quant_cv_intensity")
                pm_cv.module = "quant"
                pm_cv.render_generate_button(visualizer.plot_cv_vs_intensity, **global_plot_kwargs)
                pm_cv.render_plot_and_editor()

            with corr_inner_tab:
                st.markdown("### Sample Correlation Matrix")
                st.info("Pearson/Spearman r between all sample pairs, reordered by hierarchical clustering. Upper triangle masked for clarity; values annotated for ≤30 samples.")
                corr_method = st.selectbox(
                    "Correlation method:",
                    options=['pearson', 'spearman'],
                    key="quant_corr_method",
                )
                pm_corr = PlotManager("quant_corr")
                pm_corr.module = "quant"
                pm_corr.render_generate_button(
                    visualizer.plot_correlation_matrix,
                    method=corr_method,
                    **global_plot_kwargs,
                )
                pm_corr.render_plot_and_editor()
            
        # --- TAB 4: CLUSTERING (Updated) ---
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
                    
                    show_labels_anno = st.toggle("Show sample labels", key="pca_anno_labels")
                    
                    pm_pca_anno = PlotManager("quant_pca_anno")
                    pm_pca_anno.module = "quant"
                    pm_pca_anno.render_generate_button(
                        visualizer.plot_pca_by_annotation,
                        color_by=color_by, 
                        symbol_by=symbol_by,
                        show_labels=show_labels_anno,
                        **global_plot_kwargs
                    )
                    pm_pca_anno.render_plot_and_editor()

            with pca_cluster_tab:
                st.markdown("#### PCA with Unsupervised Clustering")
                st.info("Perform k-means clustering on the samples and visualize the results in a PCA plot.")
                
                n_clusters = st.number_input(
                    "Number of clusters (k):", min_value=2, max_value=10, value=3, step=1
                )
                
                show_labels_cluster = st.toggle("Show sample labels", key="pca_cluster_labels")

                pm_pca_cluster = PlotManager("quant_pca_cluster")
                pm_pca_cluster.module = "quant"
                pm_pca_cluster.render_generate_button(
                    visualizer.plot_pca_with_clusters,
                    n_clusters=n_clusters,
                    show_labels=show_labels_cluster,
                    **global_plot_kwargs
                )
                pm_pca_cluster.render_plot_and_editor()

            with dendro_tab:
                st.markdown("#### Hierarchical Clustering Dendrogram")
                st.info("Visualize how samples cluster together based on their protein expression profiles.")

                method = st.selectbox(
                    "Clustering method:",
                    options=['ward', 'complete', 'average', 'single'],
                    key="quant_dendro_method",
                )

                mpl_dendro = MplPlotManager("quant_dendro")
                mpl_dendro.module = "quant"
                mpl_dendro.render_generate_button(visualizer.plot_dendrogram, method=method)
                mpl_dendro.render_plot_and_editor()

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