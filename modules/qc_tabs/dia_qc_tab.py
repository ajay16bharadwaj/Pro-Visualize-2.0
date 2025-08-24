import streamlit as st
import logging
from tempfile import NamedTemporaryFile
from visualizations.DiaQcVisualizer import DiaQcVisualizer
from utils.helpers import handle_plotting_errors
from utils.caching import load_dia_visualizer


# Set up a logger for this module
logger = logging.getLogger(__name__)

class DiaQcTab:
    """
    A class to encapsulate the Streamlit UI and logic for the DIA QC tab,
    now using a tab-based layout for workflow management.
    """
    
    def __init__(self):
        """Initializes the DiaQcTab and its session state."""
        if 'dia_qc_visualizer' not in st.session_state:
            st.session_state.dia_qc_visualizer = None
        if 'dia_metadata_confirmed' not in st.session_state:
            st.session_state.dia_metadata_confirmed = False


    def _upload_data_section(self):
        """Handles the file uploader and instantiation of the visualizer."""
        with st.expander("Upload DIA-NN Report", expanded=(st.session_state.dia_qc_visualizer is None)):
            dia_report_file = st.file_uploader(
                "Upload Precursor-level Report (.parquet)",
                type=['parquet'],
                key="dia_report_uploader"
            )

            if st.button("Load DIA Data", use_container_width=True):
                # Reset confirmation flag if a new file is loaded
                st.session_state.dia_metadata_confirmed = False
                
                if dia_report_file:
                    with NamedTemporaryFile(delete=False, suffix=".parquet") as tmp:
                        tmp.write(dia_report_file.getvalue())
                        temp_path = tmp.name
                    
                    try:
                        with st.spinner("Processing DIA-NN report..."):
                            #st.session_state.dia_qc_visualizer = DiaQcVisualizer(filepath=temp_path)
                            st.session_state.dia_qc_visualizer = load_dia_visualizer(temp_path)
                        st.success("DIA-NN report loaded. Please review the extracted metadata in the first tab.")
                        st.rerun()
                    except (FileNotFoundError, ValueError, KeyError, IOError) as e:
                        st.error(f"File Processing Error: {e}", icon="📄")
                        logger.error(f"Error initializing DiaQcVisualizer: {e}", exc_info=True)
                        st.session_state.dia_qc_visualizer = None
                else:
                    st.warning("Please upload a DIA-NN .parquet report.")

    @handle_plotting_errors
    def _metadata_editor_tab(self):
        """Displays an editable table of the metadata and its summary."""
        st.info("The metadata below was automatically extracted. Edit any incorrect values, and the summary will update instantly.", icon="✏️")
        
        visualizer = st.session_state.dia_qc_visualizer
        
        # --- The Editable Data Frame ---
        edited_df = st.data_editor(
            visualizer.get_metadata(),
            num_rows="dynamic",
            key="dia_metadata_editor",
            use_container_width=True
        )

        st.markdown("---")
        
        # --- NEW: Dynamic Metadata Summary ---
        st.subheader("Metadata Summary")
        
        # Get the summary based on the *current* state of the edited data
        summary = visualizer.summarize_metadata(edited_df)

        if "error" in summary:
            st.error(summary["error"])
        else:
            # Display summary metrics in columns
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Samples Analyzed", summary['total_samples'])
            col2.metric("Unique Acquisition Days", summary['unique_dates'])
            col3.metric("Acquisition Period", f"{summary['start_date']} to {summary['end_date']}")
            
            # Display grouped dataframes in expanders
            with st.expander("View Detailed Breakdowns"):
                sub_col1, sub_col2 = st.columns(2)
                with sub_col1:
                    st.markdown("**Samples by Project**")
                    st.dataframe(summary['project_split'], use_container_width=True)
                    st.markdown("**Samples by Amount**")
                    st.dataframe(summary['amount_split'], use_container_width=True)
                with sub_col2:
                    st.markdown("**Samples per Date**")
                    st.dataframe(summary['samples_per_date'], use_container_width=True)
            
            with st.expander("View Special Run Designations"):
                sub_col1, sub_col2, sub_col3, sub_col4 = st.columns(4)
                sub_col1.metric("Check Samples", summary['check_runs'])
                sub_col2.metric("Test Samples", summary['test_runs'])
                sub_col3.metric("Bad Injections", summary['bad_inj_runs'])
                sub_col4.metric("New Column Runs", summary['new_col_runs'])


        st.markdown("---")

        # --- The Confirmation Button ---
        if st.button("Confirm and Apply Metadata Changes", use_container_width=True, type="primary"):
            try:
                visualizer.set_metadata(edited_df)
                st.session_state.dia_metadata_confirmed = True
                st.success("Metadata updated successfully!")
                st.rerun()
            except Exception as e:
                st.error(f"Failed to apply metadata changes: {e}")
                st.session_state.dia_metadata_confirmed = False

    @handle_plotting_errors           
    def _sentinel_peptides_tab(self):
        """UI for finding and displaying sentinel peptides."""
        if not st.session_state.dia_metadata_confirmed:
            st.warning("Please confirm your metadata to find sentinel peptides.", icon="⚠️")
            return

        st.info(
            "Find a stable set of **sentinel peptides** to monitor across QC metrics. "
            "These are high-abundance precursors detected in a high percentage of runs."
        )

        # --- User Controls ---
        st.subheader("1. Set Parameters")
        col1, col2 = st.columns(2)
        with col1:
            min_detection_rate_pct = st.slider(
                "Minimum Detection Rate",
                min_value=0, max_value=100, value=90, step=1,
                format="%d%%", # Format as an integer percentage
                help="A peptide must be detected in at least this percentage of runs to be considered."
            )
        with col2:
            top_n = st.number_input(
                "Number of Sentinel Peptides to Find",
                min_value=1, max_value=100, value=15, step=1,
                help="The number of most abundant peptides to select from the stable set."
            )

        # --- Execution and Display ---
        st.subheader("2. Find and Review Peptides")
        if st.button("Find Sentinel Peptides", use_container_width=True, type="primary"):
            visualizer = st.session_state.dia_qc_visualizer
            try:
                with st.spinner("Analyzing precursor stability..."):
                    sentinel_peptides = visualizer.find_stable_peptides(
                        top_n=top_n,
                        min_detection_rate=(min_detection_rate_pct / 100.0)
                    )
                # Store the result in session state for other tabs to use
                st.session_state.sentinel_peptides = sentinel_peptides
                st.success(f"Successfully identified {len(sentinel_peptides)} sentinel peptides!")

            except (ValueError, KeyError) as e:
                st.error(f"An error occurred: {e}", icon="❗")
                st.session_state.sentinel_peptides = None # Clear previous results on error

        # Display the found peptides if they exist in the session state
        if st.session_state.get('sentinel_peptides'):
            st.markdown("**Identified Sentinel Peptides:**")
            st.dataframe(st.session_state.sentinel_peptides, use_container_width=True)
            st.caption("This list is now available for plotting in the other QC tabs.")

    def _im_qc_tab(self):
        """UI for displaying all ion mobility related QC plots."""
        if not st.session_state.get('sentinel_peptides'):
            st.info("Please select sentinel peptides in the 'Sentinel Peptides' tab first to enable these plots.", icon="🧬")
            return

        visualizer = st.session_state.dia_qc_visualizer
        q_value = st.session_state.get('q_value_cutoff', 0.01)

        # Create tabs for each plot
        control_tab, drift_tab = st.tabs([
            "📈 IM Control Chart", 
            "💨 IM Drift"
        ])

        with control_tab:
            st.markdown("Monitor the ion mobility of a single peptide over time. This helps identify instrument shifts or calibration issues.")
            
            selected_peptide = st.selectbox(
                "Select a Sentinel Peptide to Monitor:",
                options=st.session_state.sentinel_peptides,
                key="im_control_peptide_select" # Use a unique key
            )
            
            if selected_peptide:
                try:
                    fig = visualizer.plot_control_chart(
                        peptide_id=selected_peptide,
                        metric_col='IM',
                        y_axis_title='Ion Mobility (1/K0)',
                        q_value_cutoff=q_value
                    )
                    st.plotly_chart(fig, use_container_width=True)
                except (ValueError, KeyError) as e:
                    st.error(f"Could not generate IM control chart: {e}")

        with drift_tab:
            st.markdown("Visualize the ion mobility deviation for multiple sentinel peptides simultaneously.")
            
            num_peptides_to_plot = st.slider(
                "Number of Peptides to Display",
                min_value=1,
                max_value=len(st.session_state.sentinel_peptides),
                value=min(5, len(st.session_state.sentinel_peptides)),
                step=1,
                key="im_drift_slider", # Use a unique key
                help="Select how many of the top sentinel peptides to show in the plot."
            )
            
            peptides_to_plot = st.session_state.sentinel_peptides[:num_peptides_to_plot]

            try:
                fig = visualizer.plot_im_drift(peptides_to_plot, q_value)
                st.plotly_chart(fig, use_container_width=True)
            except (ValueError, KeyError) as e:
                st.error(f"Could not generate IM drift plot: {e}")

    @handle_plotting_errors
    def _rt_qc_tab(self):
        """UI for displaying all retention time related QC plots using a tabbed layout."""
        if not st.session_state.get('sentinel_peptides'):
            st.info("Please select sentinel peptides in the 'Sentinel Peptides' tab first to enable these plots.", icon="🧬")
            return

        visualizer = st.session_state.dia_qc_visualizer
        q_value = st.session_state.get('q_value_cutoff', 0.01)

        # Add a new tab for the elution plot
        control_tab, drift_tab, error_tab, peak_width_tab, elution_tab = st.tabs([
            "📈 RT Control Chart", 
            "💨 RT Drift", 
            "📉 RT Prediction Error", 
            "📊 Peak Width",
            "⏳ Peptide Elution" 
        ])

        # --- Tab 1: Control Chart ---
        with control_tab:
            st.markdown("Monitor the retention time of a single peptide over time. This Levey-Jennings plot helps identify shifts or increased variability.")
            
            selected_peptide = st.selectbox(
                "Select a Sentinel Peptide to Monitor:",
                options=st.session_state.sentinel_peptides
            )
            
            if selected_peptide:
                try:
                    fig = visualizer.plot_control_chart(
                        peptide_id=selected_peptide, metric_col='RT',
                        y_axis_title='Retention Time (min)', q_value_cutoff=q_value
                    )
                    st.plotly_chart(fig, use_container_width=True)
                except (ValueError, KeyError) as e:
                    st.error(f"Could not generate control chart: {e}")

        # --- Tab 2: RT Drift ---
        with drift_tab:
            st.markdown("Visualize the retention time deviation for multiple sentinel peptides simultaneously. This helps confirm if RT shifts are systematic.")
            
            # --- NEW SLIDER ---
            num_peptides_to_plot = st.slider(
                "Number of Peptides to Display",
                min_value=1,
                max_value=len(st.session_state.sentinel_peptides),
                value=min(5, len(st.session_state.sentinel_peptides)), # Default to 5 or max available
                step=1,
                help="Select how many of the top sentinel peptides to show in the plot."
            )
            
            # Slice the list of peptides based on the slider's value
            peptides_to_plot = st.session_state.sentinel_peptides[:num_peptides_to_plot]

            try:
                fig = visualizer.plot_rt_drift(peptides_to_plot, q_value)
                st.plotly_chart(fig, use_container_width=True)
            except (ValueError, KeyError) as e:
                st.error(f"Could not generate RT drift plot: {e}")
                
        # --- Tab 3: RT Prediction Error ---
        with error_tab:
            st.markdown("Track how well the observed retention time matches the predicted RT. A consistent, low error is desirable.")
            
            moving_avg_window = st.slider(
                "Moving Average Window", min_value=1, max_value=21, value=5, step=2,
                help="Smooths the error trend over a window of N runs."
            )
            
            try:
                fig = visualizer.plot_rt_prediction_error(q_value, moving_avg_window)
                st.plotly_chart(fig, use_container_width=True)
            except (ValueError, KeyError) as e:
                st.error(f"Could not generate RT prediction error plot: {e}")
                
        # --- Tab 4: Peak Width Distribution ---
        with peak_width_tab:
            st.markdown("Assess chromatographic performance by visualizing the distribution of Full Width at Half Maximum (FWHM) for all high-confidence peptides on each day.")
            try:
                fig = visualizer.plot_peak_width_distribution(q_value)
                st.plotly_chart(fig, use_container_width=True)
            except (ValueError, KeyError) as e:
                st.error(f"Could not generate peak width plot: {e}")
        
        with elution_tab:
            st.markdown("Visualize how evenly the unique identified peptides are distributed across the entire retention time gradient.")
            
            # Interactive widget for histogram granularity
            num_bins = st.slider(
                "Number of Bins",
                min_value=25, max_value=250, value=150, step=25,
                help="Adjust the number of bins to make the distribution plot broader or more granular."
            )
            
            try:
                fig = visualizer.plot_peptide_elution_distribution(q_value, num_bins)
                st.plotly_chart(fig, use_container_width=True)
            except (ValueError, KeyError) as e:
                st.error(f"Could not generate peptide elution plot: {e}")

    @handle_plotting_errors
    def _mass_accuracy_qc_tab(self):
        """UI for displaying all mass accuracy related QC plots."""
        if not st.session_state.get('sentinel_peptides'):
            st.info("Please select sentinel peptides in 'Sentinel Peptides' tab first.", icon="🧬")
            return

        visualizer = st.session_state.dia_qc_visualizer
        q_value = st.session_state.get('q_value_cutoff', 0.01)

        dist_tab, sentinel_tab, trend_tab = st.tabs([
            "📊 Mass Error Distribution", 
            "📈 Sentinel Mass Accuracy", 
            "📉 Mass Error Trend"
        ])

        with dist_tab:
            st.markdown("Assess the overall mass accuracy by visualizing the distribution of mass error (in ppm) for all high-confidence peptides on each day.")
            
            y_range = st.slider(
                "Set Y-axis Range (ppm)", min_value=-20.0, max_value=20.0, 
                value=(-5.0, 5.0), step=0.5,
                help="Zoom in or out to focus on the most relevant mass error range."
            )
            
            try:
                fig = visualizer.plot_mass_accuracy_distribution(q_value, list(y_range))
                st.plotly_chart(fig, use_container_width=True)
            except (ValueError, KeyError) as e:
                st.error(f"Could not generate plot: {e}")

        with sentinel_tab:
            st.markdown("Track the mass accuracy for your selected sentinel peptides across all runs to spot systematic drift or inconsistencies.")
            
            num_peptides = st.slider(
                "Number of Peptides to Display", 1, len(st.session_state.sentinel_peptides),
                min(5, len(st.session_state.sentinel_peptides)), 1,
                key="mass_acc_sentinel_slider",
                help="Select how many sentinel peptides to show."
            )
            peptides_to_plot = st.session_state.sentinel_peptides[:num_peptides]

            try:
                fig = visualizer.plot_sentinel_mass_accuracy(peptides_to_plot, q_value)
                st.plotly_chart(fig, use_container_width=True)
            except (ValueError, KeyError) as e:
                st.error(f"Could not generate plot: {e}")

        with trend_tab:
            st.markdown("Visualize the smoothed trend of the median mass error per run. This helps to see long-term instrument calibration drift more clearly.")
            
            window = st.slider(
                "Moving Average Window", 1, 21, 5, 2,
                key="mass_error_trend_slider",
                help="Smooths the error trend over N runs."
            )
            
            try:
                fig = visualizer.plot_mass_error_trend(q_value, window)
                st.plotly_chart(fig, use_container_width=True)
            except (ValueError, KeyError) as e:
                st.error(f"Could not generate plot: {e}")

    @handle_plotting_errors
    def _download_center_tab(self):
        """Provides download functionality for the processed data."""
        if not st.session_state.dia_metadata_confirmed:
            st.warning("Please confirm your metadata in the 'Metadata Review & Edit' tab to enable downloads.", icon="⚠️")
            return

        st.info("Download the full precursor-level report combined with your finalized metadata.", icon="💾")
        
        visualizer = st.session_state.dia_qc_visualizer
        processed_data = visualizer.get_processed_data()
        
        # Convert dataframe to CSV for download
        csv_data = processed_data.to_csv(index=False).encode('utf-8')
        
        st.download_button(
            label="Download Data as CSV",
            data=csv_data,
            file_name="pro_visualize_processed_data.csv",
            mime="text/csv",
            use_container_width=True
        )

    @handle_plotting_errors
    def _qc_visualizations_tab(self):
        """Displays the QC plots after metadata is confirmed."""
        if not st.session_state.dia_metadata_confirmed:
            st.warning("Please confirm your metadata in the 'Metadata Review & Edit' tab to view visualizations.", icon="⚠️")
            return

        #st.info("This is where your QC plots will appear. You can scale the application by adding new plots here.")
        # --- MASTER Q-VALUE FILTER ---
        st.subheader("Global Plot Filter")
        st.session_state.q_value_cutoff = st.number_input(
            "Q-Value Cutoff for all Plots",
            min_value=0.001, max_value=0.1, value=0.01, step=0.001, format="%.3f",
            help="This filter is applied to all plots below to ensure only high-confidence data is shown."
        )
        st.markdown("---")

        # Placeholder for future plots
        visualizer = st.session_state.dia_qc_visualizer
        processed_data = visualizer.get_processed_data()
        
        #st.markdown()
        with st.expander(f"**Data Preview:** The full dataset with **{len(processed_data)}** rows is ready."):
            st.dataframe(processed_data.head())

        # --- Create the 4 tabs for QC visualizations ---
        (
            sentinel_tab,
            rt_tab,
            im_tab,
            mass_accuracy_tab
        ) = st.tabs([
            "🧬 Sentinel Peptides",
            "⏱️ Retention Time",
            "💨 Ion Mobility",
            "🎯 Mass Accuracy"
        ])

        with sentinel_tab:
            self._sentinel_peptides_tab() # Call the new method here

        with rt_tab:
            #st.info("This tab will show Retention Time stability for the selected sentinel peptides. (Under Construction)")
            self._rt_qc_tab()

        with im_tab:
            #st.info("This tab will show Ion Mobility stability for the selected sentinel peptides. (Under Construction)")
            self._im_qc_tab()
        
        with mass_accuracy_tab:
            #st.info("This tab will show Mass Accuracy for the selected sentinel peptides. (Under in Construction)")
            self._mass_accuracy_qc_tab()

    

    def render(self):
        """Renders the entire DIA QC tab workflow."""
        # Initialize sentinel peptides in session state if it doesn't exist
        if 'sentinel_peptides' not in st.session_state:
            st.session_state.sentinel_peptides = None

        st.subheader("DIA Quality Control")
        
        self._upload_data_section()

        if st.session_state.dia_qc_visualizer:
            meta_tab, plots_tab, download_tab = st.tabs([
                "1. Metadata Review & Edit", 
                "2. QC Visualizations",
                "3. Download Center"
            ])

            with meta_tab:
                self._metadata_editor_tab()
            
            with plots_tab:
                self._qc_visualizations_tab()

            with download_tab:
                self._download_center_tab()