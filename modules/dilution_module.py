import streamlit as st
import pandas as pd
import logging
from visualizations.dilution_series import DilutionSeriesVisualizer
from utils.helpers import handle_plotting_errors

logger = logging.getLogger(__name__)

class DilutionSeriesTab:
    def __init__(self):
        if 'dilution_visualizer' not in st.session_state:
            st.session_state.dilution_visualizer = None

    def _upload_data_section(self):
        # (This method remains the same)
        with st.expander("Upload Data for Dilution Series", expanded=(st.session_state.dilution_visualizer is None)):
            protein_file = st.file_uploader("Upload Protein-Level Data (.txt, .csv, .tsv)", type=['txt', 'csv', 'tsv'], key="dilution_protein_uploader")
            metadata_file = st.file_uploader("Upload Metadata File (.csv, .tsv)", type=['csv', 'tsv'], key="dilution_metadata_uploader")

            if st.button("Load & Process Dilution Data", use_container_width=True):
                if protein_file and metadata_file:
                    try:
                        # Infer separator for flexibility
                        protein_df = pd.read_csv(protein_file, sep=None, engine='python')
                        metadata_df = pd.read_csv(metadata_file, sep=None, engine='python')
                        
                        with st.spinner("Initializing visualizer..."):
                            st.session_state.dilution_visualizer = DilutionSeriesVisualizer(protein_df, metadata_df)
                        st.success("Data loaded successfully!")
                        st.rerun()
                    except Exception as e:
                        st.error(f"An error occurred: {e}")
                        logger.error(f"Dilution Series Error: {e}", exc_info=True)
                else:
                    st.warning("Please upload both data and metadata files.")

        

    @handle_plotting_errors
    def _display_plots(self):
        visualizer = st.session_state.dilution_visualizer
        st.subheader("Dilution Series Visualizations")

        overview_tab, analysis_tab, pca_tab, completeness_tab = st.tabs(["📊 Data Overview", "📈 Trend Analysis", "🧬 PCA", "✅ Completeness"])

        with overview_tab:
            # (This section remains the same)
            st.markdown("### Protein Counts per Sample")
            fig_counts = visualizer.plot_protein_counts_per_sample()
            st.plotly_chart(fig_counts, use_container_width=True)
            st.markdown("### Intensity and CV Distributions")
            col1, col2 = st.columns(2)
            with col1:
                plot_type = st.selectbox("Intensity Plot Type", ["box", "violin"], key="dist_plot_type")
                fig_dist = visualizer.plot_intensity_distribution(plot_type=plot_type, points=False)
                st.plotly_chart(fig_dist, use_container_width=True)
            with col2:
                fig_cv = visualizer.plot_cv_distribution()
                st.plotly_chart(fig_cv, use_container_width=True)

        with analysis_tab:
            # (This section remains the same)
            trends_tab, heatmap_tab, ratio_tab = st.tabs(["Protein Trends", "Heatmap", "Relative Abundance"])

            with trends_tab:
                st.markdown("### Individual Protein Trends")
                all_proteins = visualizer.protein_df[visualizer.protein_id_col].unique().tolist()
                selected_proteins = st.multiselect(
                    "Select specific proteins to plot (optional)", options=all_proteins,
                    help="If none are selected, the most abundant proteins will be shown."
                )
                n_top = st.slider(
                    "Number of top proteins to show (if none selected above)", 1, 20, 5,
                    disabled=bool(selected_proteins)
                )
                fig_trends = visualizer.plot_protein_trends(
                    proteins_to_plot=selected_proteins if selected_proteins else None,
                    n_top_proteins=n_top
                )
                st.plotly_chart(fig_trends, use_container_width=True)

            with heatmap_tab:
                st.markdown("### Protein Intensity Heatmap")
                apply_zscore = st.toggle("Apply Z-score normalization", value=True)
                max_proteins = st.slider("Max proteins to display", 50, 2000, 500, 50)
                fig_heatmap = visualizer.plot_heatmap_trends(
                    max_proteins_to_plot=max_proteins, apply_zscore=apply_zscore
                )
                st.plotly_chart(fig_heatmap, use_container_width=True)

            with ratio_tab:
                st.markdown("### Protein Abundance Ratios")
                st.markdown("This plot shows the `Log2` intensity ratio of each concentration group relative to the lowest concentration. The median of the boxes should align with the dashed 'Expected' lines for accurate quantification.")
                
                show_lines = st.toggle("Show expected ratio lines", value=True)
                
                fig_ratios = visualizer.plot_relative_abundance_ratios(
                    add_expected_lines=show_lines,
                    points=False
                )
                st.plotly_chart(fig_ratios, use_container_width=True)
            
        with pca_tab:
            # (This section remains the same)
            st.markdown("### Principal Component Analysis (PCA)")
            color_by = st.selectbox("Color PCA points by", options=['Group', 'Concentration', 'Replicate'], index=0)
            symbol_by = st.selectbox("Use symbols for", options=[None, 'Replicate', 'Group'], index=1)
            fig_pca = visualizer.plot_pca(color_by=color_by, symbol_by=symbol_by)
            st.plotly_chart(fig_pca, use_container_width=True)


        with completeness_tab:
            st.markdown("### Detection Completeness Analysis")
            st.info(
                "This view shows how many proteins, peptides, or precursors are detected "
                "across your dilution series, with different quality thresholds."
            )
            
            # --- UPDATED USER CONTROLS ---
            col1, col2, col3 = st.columns(3)
            
            with col1:
                # Select the level of analysis
                analysis_level = st.selectbox(
                    "Select Analysis Level",
                    options=["Proteins", "Precursors", "Stripped Sequences"],
                    help="Choose what to count: proteins, precursors, or unique peptide sequences"
                )
                
                # Map user-friendly names to column names
                level_mapping = {
                    "Proteins": visualizer.protein_id_col,  # Usually 'Protein.Group'
                    "Precursors": "Precursor.Id" if "Precursor.Id" in visualizer.protein_df.columns else visualizer.protein_id_col,
                    "Stripped Sequences": "Stripped.Sequence" if "Stripped.Sequence" in visualizer.protein_df.columns else visualizer.protein_id_col
                }
                
                selected_column = level_mapping[analysis_level]
            
            with col2:
                # NEW: Add slider for CV threshold
                cv_thresh = st.slider(
                    "CV Threshold (%)", 0.0, 50.0, 20.0, 1.0,
                    key="completeness_cv",
                    help="The CV cutoff to define 'high-quality' data."
                )
                
            with col3:
                # Toggle for log scale
                use_log = st.toggle(
                    "Use logarithmic scale (left plot)",
                    value=False,
                    help="Enable log scale if counts vary widely"
                )
            
            # --- UPDATED PLOT CALL ---
            try:
                with st.spinner(f"Calculating {analysis_level.lower()} completeness..."):
                    fig = visualizer.plot_completeness_overview(
                        identifier_col=selected_column,
                        use_log_scale=use_log,
                        cv_threshold=cv_thresh  # Pass the new CV threshold
                    )
                    st.plotly_chart(fig, use_container_width=True)
            except ValueError as e:
                st.error(f"Could not generate plot: {e}")
                st.info(
                    f"The selected analysis level ('{analysis_level}') may not be available "
                    f"in your dataset. Try a different level or check your data columns."
                )
            except Exception as e:
                st.error(f"An unexpected error occurred: {e}")

    def render(self):
        st.header("Dilution Series Analysis")
        self._upload_data_section()
        if st.session_state.dilution_visualizer:
            self._display_plots()

def render():
    """Renders the Dilution Series tab."""
    dilution_tab = DilutionSeriesTab()
    dilution_tab.render()