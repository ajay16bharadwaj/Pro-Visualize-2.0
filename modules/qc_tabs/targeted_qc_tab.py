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
        """
        Handles the file uploader and data processing logic.
        This section's error handling ensures a bad file doesn't affect other tabs.
        """
        with st.expander("Upload Targeted QC Report", expanded=(st.session_state.targeted_qc_visualizer is None)):
            targeted_file = st.file_uploader(
                "Upload Targeted Proteomics Report",
                type=['csv', 'tsv'],
                key="targeted_report_uploader"
            )

            if st.button("Load Targeted Data", use_container_width=True):
                if targeted_file:
                    # Using a temporary file to get a path for the visualizer class
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
                        st.session_state.targeted_qc_visualizer = None # Clear on failure
                else:
                    st.warning("Please upload a targeted proteomics report.")

    @handle_plotting_errors
    def _plot_retention_time_stability(self, visualizer, top_n):
        """Generates and displays the RT stability plot."""
        fig = visualizer.plot_retention_time_stability(top_n_peptides=top_n)
        st.plotly_chart(fig, use_container_width=True)

    @handle_plotting_errors
    def _plot_rt_distribution(self, visualizer):
        """Generates and displays the RT distribution plot."""
        fig = visualizer.plot_rt_distribution()
        st.plotly_chart(fig, use_container_width=True)

    @handle_plotting_errors
    def _plot_peak_area(self, visualizer, top_n):
        """Generates and displays the peak area plot."""
        fig = visualizer.plot_peak_area(top_n_peptides=top_n)
        st.plotly_chart(fig, use_container_width=True)

    @handle_plotting_errors
    def _display_peptide_stats(self, visualizer):
        """Calculates and displays the peptide statistics table."""
        stats_df = visualizer.calculate_peptide_stats()
        st.dataframe(stats_df)

    def _display_results(self):
        """Displays the QC plots and tables within nested tabs."""
        visualizer = st.session_state.targeted_qc_visualizer
        st.markdown("---")
        st.subheader("QC Visualization Dashboard")

        # Create nested tabs for different visualization categories
        rt_tab, area_tab = st.tabs(["RT Stability & Distribution", "Peak Area & Statistics"])

        with rt_tab:
            st.markdown("#### Retention Time Visualization")
            top_n_rt = st.slider(
                "Select Number of Peptides for RT Stability Plot", 
                min_value=5, max_value=50, value=30, key="top_n_rt"
            )
            # Call the decorated plotting methods
            self._plot_retention_time_stability(visualizer, top_n_rt)
            self._plot_rt_distribution(visualizer)

        with area_tab:
            st.markdown("#### Peak Area Visualization")
            top_n_area = st.slider(
                "Select Number of Peptides for Peak Area Plot", 
                min_value=5, max_value=50, value=15, key="top_n_area"
            )
            # Call the decorated plotting methods
            self._plot_peak_area(visualizer, top_n_area)
            
            st.markdown("#### Peptide Statistics Summary")
            st.markdown("Coefficient of Variation (CV) for Retention Time (RT) and Peak Area across all runs.")
            self._display_peptide_stats(visualizer)

    def render(self):
        """Renders the entire Targeted QC tab."""
        st.subheader("Targeted Quality Control")
        
        self._upload_and_process_data()

        if st.session_state.targeted_qc_visualizer:
            self._display_results()
        else:
            st.info("Upload a targeted QC report to begin analysis.")