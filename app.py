import streamlit as st

# --- Import Project Modules ---
# Import the main render function from your QC module
from modules.qc_module import render as render_qc

# Future modules will be imported here
# from modules.quant_module import render as render_quant
# from modules.comparative_module import render as render_comp

# --- Page Configuration ---
st.set_page_config(
    page_title="Pro-Visualize",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

st.title("Pro-Visualize: Proteomics Data Visualization")
st.markdown("A scalable application for multi-omics analysis, starting with targeted and DIA-MS proteomics.")

# --- Create Main Tabs for Different Analyses ---
tab_welcome, tab_qc, tab_quant, tab_comp = st.tabs([
    "👋 Welcome",
    "📊 QC Analysis", 
    "📈 Quantification", 
    "🆚 Comparative Analysis"
])

with tab_welcome:
    st.header("Welcome to Pro-Visualize!")
    st.markdown("""
    This application is designed to provide a comprehensive suite of tools for visualizing and analyzing proteomics data.
    
    **Current Modules:**
    - **QC Analysis**: Upload and analyze Quality Control data from both DIA and Targeted (e.g., PRM/SRM) experiments.
    
    Navigate to the **QC Analysis** tab to begin.
    """)

with tab_qc:
    # The render_qc function now builds the entire UI for the QC section,
    # including its own internal tabs for DIA and Targeted QC.
    render_qc()

with tab_quant:
    st.header("Protein Quantification")
    st.info("This section is currently under development. 🏗️")
    # Future: Call render_quant()

with tab_comp:
    st.header("Comparative Analysis")
    st.info("This section is currently under development. 🏗️")
    # Future: Call render_comp()