import streamlit as st
import logging
from functools import wraps
from tempfile import NamedTemporaryFile
from utils.helpers import handle_plotting_errors
from visualizations.targettedQCVisualization import TargetedQcVisualizer

# Set up a logger for this module
logger = logging.getLogger(__name__)



class TargetedQcTab:
    """
    A class to encapsulate the Streamlit UI and logic for the Targeted QC tab.
    """
    
    def __init__(self):
        """Initializes the TargetedQcTab and its session state."""
        if 'targeted_qc_visualizer' not in st.session_state:
            st.session_state.targeted_qc_visualizer = None

    def _upload_and_process_data(self):
        """Handles the file uploader and data processing logic."""
        with st.expander("Upload Targeted QC Report", expanded=(st.session_state.targeted_qc_visualizer is None)):
            targeted_file = st.file_uploader(
                "Upload Targeted Proteomics Report",
                type=['csv', 'tsv'],
                key="targeted_report_uploader"
            )

            if st.button("Load Targeted Data", use_container_width=True):
                if targeted_file:
                    suffix = f".{targeted_file.name.split('.')[-1]}"
                    with NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                        tmp.write(targeted_file.getvalue())
                        temp_path = tmp.name
                    
                    try:
                        with st.spinner("Processing Targeted QC report..."):
                            st.session_state.targeted_qc_visualizer = TargetedQcVisualizer(filepath=temp_path)
                        st.success("Targeted QC report loaded successfully!")
                        st.rerun()
                    except (FileNotFoundError, ValueError, KeyError, RuntimeError) as e:
                        st.error(f"File Processing Error: {e}", icon="📄")
                        logger.error(f"Error initializing TargetedQcVisualizer: {e}", exc_info=True)
                        st.session_state.targeted_qc_visualizer = None
                else:
                    st.warning("Please upload a targeted proteomics report.")

    # --- Methods for Existing Plots (No Changes Here) ---
    @handle_plotting_errors
    def _plot_retention_time_stability(self, visualizer, top_n):
        fig = visualizer.plot_retention_time_stability(top_n_peptides=top_n, plot_config=st.session_state.get('plot_config'))
        st.plotly_chart(fig, use_container_width=True)

    @handle_plotting_errors
    def _plot_rt_distribution(self, visualizer):
        fig = visualizer.plot_rt_distribution(plot_config=st.session_state.get('plot_config'))
        st.plotly_chart(fig, use_container_width=True)

    @handle_plotting_errors
    def _plot_peak_area(self, visualizer, top_n):
        fig = visualizer.plot_peak_area(top_n_peptides=top_n, plot_config=st.session_state.get('plot_config'))
        st.plotly_chart(fig, use_container_width=True)

    @handle_plotting_errors
    def _display_peptide_stats(self, visualizer):
        stats_df = visualizer.calculate_peptide_stats()
        st.dataframe(stats_df)

    # --- NEW Methods for CV Distribution Plot and Failing Peptides Table ---
    @handle_plotting_errors
    def _plot_cv_distributions(self, visualizer, area_cv_threshold, rt_cv_threshold):
        """Generates and displays the CV distribution histograms."""
        fig = visualizer.plot_cv_distributions(
            area_cv_threshold=area_cv_threshold,
            rt_cv_threshold=rt_cv_threshold,
            plot_config=st.session_state.get('plot_config')
        )
        st.plotly_chart(fig, use_container_width=True)

    @handle_plotting_errors
    def _display_failing_peptides(self, visualizer, area_cv_threshold, rt_cv_threshold):
        """Gets and displays the table of peptides failing QC."""
        failing_df = visualizer.get_failing_peptides(
            area_cv_threshold=area_cv_threshold,
            rt_cv_threshold=rt_cv_threshold
        )
        if failing_df.empty:
            st.success("All peptides passed the specified QC thresholds! ✅")
        else:
            st.warning(f"Found {len(failing_df)} peptides failing QC thresholds.", icon="⚠️")
            st.dataframe(failing_df)

    def _display_results(self):
        """Displays the QC plots and tables within nested tabs."""
        visualizer = st.session_state.targeted_qc_visualizer
        st.markdown("---")
        st.subheader("QC Visualization Dashboard")

        rt_tab, area_tab = st.tabs(["RT Stability & Distribution", "Peak Area & Reproducibility"])

        # --- RT Tab (No Changes Here) ---
        with rt_tab:
            st.markdown("#### Retention Time Visualization")
            top_n_rt = st.slider(
                "Select Number of Peptides for RT Stability Plot", 
                min_value=5, max_value=50, value=30, key="top_n_rt"
            )
            self._plot_retention_time_stability(visualizer, top_n_rt)
            self._plot_rt_distribution(visualizer)

        # --- Peak Area & Statistics Tab (UPDATED) ---
        with area_tab:
            st.markdown("#### Peak Area Visualization")
            top_n_area = st.slider(
                "Select Number of Peptides for Peak Area Plot", 
                min_value=5, max_value=50, value=15, key="top_n_area"
            )
            self._plot_peak_area(visualizer, top_n_area)
            
            st.markdown("---") # Visual separator
            
            # --- NEW Reproducibility Analysis Section ---
            st.markdown("#### Measurement Reproducibility Analysis")
            st.markdown(
                "Analyze the Coefficient of Variation (CV) for peak area and retention time across all runs. "
                "Set the CV thresholds below to dynamically update the plots and tables."
            )
            
            # Interactive controls for CV thresholds
            col1, col2 = st.columns(2)
            with col1:
                area_cv_threshold = st.number_input(
                    "Area CV Threshold (%)", min_value=0.0, value=20.0, step=1.0,
                    help="Peptides with an Area CV above this value will be flagged."
                )
            with col2:
                rt_cv_threshold = st.number_input(
                    "RT CV Threshold (%)", min_value=0.0, value=2.0, step=0.5,
                    help="Peptides with a Retention Time CV above this value will be flagged."
                )

            # Call the new methods to display the plot and table
            self._plot_cv_distributions(visualizer, area_cv_threshold, rt_cv_threshold)
            self._display_failing_peptides(visualizer, area_cv_threshold, rt_cv_threshold)

            st.markdown("---") # Visual separator
            
            st.markdown("#### Peptide Statistics Summary")
            st.markdown("Overall summary table of calculated peptide statistics.")
            self._display_peptide_stats(visualizer)

    def render(self):
        """Renders the entire Targeted QC tab."""
        st.subheader("Targeted Quality Control")
        self._upload_and_process_data()

        if st.session_state.targeted_qc_visualizer:
            self._display_results()
        else:
            st.info("Upload a targeted QC report to begin analysis.")