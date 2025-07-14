# pro-visualize/modules/qc_module.py

import streamlit as st
from utils.data_manager import DataManager
from visualizations.qc_plots import QCPlots
from config.plot_configs import QC_PLOTS_CONFIG

# Import the render functions for the sub-tabs
from .qc_tabs.dimensionality_reduction_tab import render as render_dim_reduction
from .qc_tabs.correlation_tab import render as render_correlation
from .qc_tabs.data_distribution_tab import render as render_distribution

def render():
    """Renders the entire QC module, including data uploads and analysis tabs."""
    st.header("Quality Control Analysis")

    # --- Initialize Session State for this module ---
    if 'qc_plotter' not in st.session_state:
        st.session_state.qc_plotter = None

    # --- Data Upload Section for QC ---
    with st.expander("Upload Data for QC", expanded=(st.session_state.qc_plotter is None)):
        protein_file = st.file_uploader(
            "Upload Protein Level Data",
            type=['csv', 'tsv', 'txt'],
            key="qc_protein_uploader"
        )
        
        annotation_file = st.file_uploader(
            "Upload Annotation File",
            type=['csv', 'tsv', 'txt'],
            key="qc_annotation_uploader"
        )

        if st.button("Load Data for QC", use_container_width=True):
            if protein_file and annotation_file:
                data_manager = DataManager()
                protein_df = data_manager.load_protein_data(protein_file)
                annotation_df = data_manager.load_annotation_data(annotation_file)
                
                if protein_df is not None and annotation_df is not None:
                    with st.spinner("Processing QC data..."):
                        st.session_state.qc_plotter = QCPlots(protein_df, annotation_df)
                    st.success("QC data loaded!")
                else:
                    st.error("Failed to load QC data. Please check files.")
            else:
                st.warning("Please upload both files for QC analysis.")

    # --- Analysis Section ---
    if st.session_state.qc_plotter is not None:
        # Create the sub-tabs for QC visualizations
        tab1, tab2, tab3 = st.tabs([
            "📊 Dimensionality Reduction",
            "🔗 Sample Correlation",
            "🧬 Data Distribution"
        ])

        with tab1:
            render_dim_reduction(st.session_state.qc_plotter, QC_PLOTS_CONFIG)

        with tab2:
            render_correlation(st.session_state.qc_plotter, QC_PLOTS_CONFIG.get('correlation_matrix', {}))

        with tab3:
            render_distribution(st.session_state.qc_plotter, QC_PLOTS_CONFIG.get('intensity_density', {}))
    else:
        st.info("Please upload your data above to begin the Quality Control analysis.")