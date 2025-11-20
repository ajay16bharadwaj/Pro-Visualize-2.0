import streamlit as st
import warnings

warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    message="Downcasting object dtype arrays on .fillna"
)
# --- Import Project Modules ---
# Import the main render function from your QC module
from modules.qc_module import render as render_qc
from modules.dilution_module import render as render_dilution
from modules.quant_module import render as render_quant
from modules.comparative_module import render as render_comp
#from modules.rag_module import render as render_rag


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
tab_welcome, tab_qc, tab_dilution, tab_quant, tab_comp, tab_chat = st.tabs([
    "👋 Welcome",
    "📊 QC Analysis", 
    "📈 Dilution Series",
    "📈 Quantification", 
    "🆚 Comparative Analysis",
    "💬 Pro-Viz Chat"
])

with tab_welcome:
    st.header("Welcome to Pro-Visualize!")
    st.markdown("""
    This application is designed to provide a comprehensive suite of tools for visualizing and analyzing proteomics data.
    
    **Current Modules:**
    - **QC Analysis**: Upload and analyze Quality Control data from both DIA and Targeted (e.g., PRM/SRM) experiments.
    - **Dilution Series (Loading-curve)**: Upload and analyze loading curve data from DIA experiments in Protein and Precursor Level.   
    - **Quantification**: Upload protein-level data to perform quantitative analysis. 
    - **Comparative Analysis**: Visualize statistical comparisons between groups.        
    
    Navigate to the desired tab to begin.
    """)

with tab_qc:
    # The render_qc function now builds the entire UI for the QC section,
    # including its own internal tabs for DIA and Targeted QC.
    render_qc()

with tab_dilution:
    #st.header("Dilution Series Analysis")
    #st.info("This section is currently under development. 🏗️")
    # Future: 
    render_dilution()

with tab_quant:
    #st.header("Protein Quantification")
    #st.info("This section is currently under development. 🏗️")
    # Future: Call 
    render_quant()

with tab_comp:
    #st.header("Comparative Analysis")
    #st.info("This section is currently under development. 🏗️")
    # Future: Call 
    render_comp()

with tab_chat:  
    #render_rag()
    st.info("This section is currently under development. 🏗️")