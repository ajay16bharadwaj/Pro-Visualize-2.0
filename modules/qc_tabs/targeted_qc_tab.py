import streamlit as st
import logging
from tempfile import NamedTemporaryFile
from utils.helpers import handle_plotting_errors
from utils.caching import load_targeted_visualizer

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
                            #st.session_state.targeted_qc_visualizer = TargetedQcVisualizer(filepath=temp_path)
                            st.session_state.targeted_qc_visualizer = load_targeted_visualizer(temp_path)
                        st.success("Targeted QC report loaded successfully!")
                        st.rerun()
                    except (FileNotFoundError, ValueError, KeyError, RuntimeError) as e:
                        st.error(f"File Processing Error: {e}", icon="📄")
                        st.session_state.targeted_qc_visualizer = None
                else:
                    st.warning("Please upload a targeted proteomics report.")

    # --- UI Wrapper Methods for each plot ---
    @handle_plotting_errors
    def _display_system_health(self, visualizer):
        st.markdown("#### System-Wide Performance Metrics")
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Total Ion Current (TIC)**")
            st.markdown("Tracks overall signal intensity. Helps spot loading issues or MS sensitivity drops.")
            # Add radio button to toggle view
            tic_view_by = st.radio("View TIC By:", ('run', 'day'), key='tic_view', horizontal=True)
            fig_tic = visualizer.plot_tic(view_by=tic_view_by)
            st.plotly_chart(fig_tic, use_container_width=True)
            
        with col2:
            st.markdown("**Daily Reproducibility Trend**")
            st.markdown("Tracks the median peptide CV per day. Helps spot gradual loss of performance.")
            fig_cv_trend = visualizer.plot_daily_cv_trends()
            st.plotly_chart(fig_cv_trend, use_container_width=True)

    @handle_plotting_errors
    def _display_peptide_stability(self, visualizer):
        st.markdown("#### Peptide Stability Across Runs and Days")
        st.markdown(
            "Use the 'View By' toggle to switch between a detailed run-by-run view and a summarized daily trend view."
        )

        col1, col2 = st.columns(2)
        with col1:
            view_by_rt = st.radio("View RT Stability By:", ('run', 'day'), key='rt_view', horizontal=True)
            top_n_rt = st.slider("Select # of Peptides", 5, 50, 30, key="top_n_rt")
            fig = visualizer.plot_retention_time_stability(top_n_peptides=top_n_rt, view_by=view_by_rt)
            st.plotly_chart(fig, use_container_width=True)
            
        with col2:
            view_by_area = st.radio("View Peak Area By:", ('run', 'day'), key='area_view', horizontal=True)
            top_n_area = st.slider("Select # of Peptides", 5, 50, 15, key="top_n_area")
            fig = visualizer.plot_peak_area(top_n_peptides=top_n_area, view_by=view_by_area)
            st.plotly_chart(fig, use_container_width=True)
        
        st.markdown("#### Peptide Elution Profile")
        fig_dist = visualizer.plot_rt_distribution()
        st.plotly_chart(fig_dist, use_container_width=True)

    @handle_plotting_errors
    def _display_reproducibility(self, visualizer):
        st.markdown("#### Measurement Reproducibility Analysis")
        st.markdown(
            "Set the CV thresholds below to see the overall distribution and identify specific peptides that fail the QC criteria."
        )
        col1, col2 = st.columns(2)
        with col1:
            area_cv_threshold = st.number_input("Area CV Threshold (%)", 0.0, 100.0, 20.0, 1.0)
        with col2:
            rt_cv_threshold = st.number_input("RT CV Threshold (%)", 0.0, 100.0, 2.0, 0.5)
        
        fig_cv = visualizer.plot_cv_distributions(area_cv_threshold, rt_cv_threshold)
        st.plotly_chart(fig_cv, use_container_width=True)
        
        failing_df = visualizer.get_failing_peptides(area_cv_threshold, rt_cv_threshold)
        if failing_df.empty:
            st.success("All peptides passed the specified QC thresholds! ✅")
        else:
            st.warning(f"Found {len(failing_df)} peptides failing QC thresholds.", icon="⚠️")
            st.dataframe(failing_df)

    @handle_plotting_errors
    def _display_investigation(self, visualizer):
        st.markdown("#### Peptide Investigation: Control Charts")
        st.markdown(
            "Select a specific peptide to generate a detailed control chart. "
            "This is useful for monitoring internal standards or key biomarkers."
        )
        peptide_list = visualizer.data['Peptide'].unique()
        selected_peptide = st.selectbox("Select a Peptide to Investigate", options=peptide_list)
        
        if selected_peptide:
            metric = st.radio("Select Metric", ('Peak Area', 'Retention Time'), horizontal=True, key='control_metric')
            view_by = st.radio("View By", ('run', 'day', 'day_distribution'), horizontal=True, key='control_view')
            
            fig = visualizer.plot_qc_peptide_chart(
                peptide_sequence=selected_peptide,
                metric=metric,
                view_by=view_by
            )
            st.plotly_chart(fig, use_container_width=True)

    def _display_results(self):
        """Displays the QC plots and tables within a logical tab structure."""
        visualizer = st.session_state.targeted_qc_visualizer
        st.markdown("---")
        st.subheader("QC Visualization Dashboard")

        # Create the three main tabs for the analysis workflow
        tab1, tab2, tab3 = st.tabs([
            "📈 System Health ", 
            "🔬 Peptide Stability", 
            "🎯 Reproducibility"
        ])

        with tab1:
            self._display_system_health(visualizer)

        with tab2:
            self._display_peptide_stability(visualizer)
            
        with tab3:
            self._display_investigation(visualizer)
            st.markdown("---")
            self._display_reproducibility(visualizer)
            st.markdown("---")
            st.markdown("#### Peptide Statistics Summary Table")
            stats_df = visualizer.calculate_peptide_stats()
            st.dataframe(stats_df)

    def render(self):
        """Renders the entire Targeted QC tab."""
        st.header("Targeted Proteomics Quality Control")
        self._upload_and_process_data()

        if st.session_state.targeted_qc_visualizer:
            self._display_results()
        else:
            st.info("Upload a targeted QC report to begin analysis.")