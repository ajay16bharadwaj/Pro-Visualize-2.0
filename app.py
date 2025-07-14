# pro-visualize/app.py

import streamlit as st

# --- Import Project Modules ---
# Note: We no longer import DataManager or Plotters here
#from modules.qc_module import render as render_qc
# from modules.quant_module import render as render_quant
# from modules.comparative_module import render as render_comp

# --- Page Configuration ---
st.set_page_config(
    page_title="Pro-Visualize",
    page_icon="🔬",
    layout="wide"
)

st.title("Pro-Visualize: Proteomics Data Visualization")

# --- Create Main Tabs ---
Welcome, tab_qc, tab_quant, tab_comp = st.tabs([
    "Welcome",
    "QC Analysis", 
    "Quantification", 
    "Comparative Analysis"
])

with Welcome:
    st.info("Generic information on how Pro-Visualize can be used")
    st.header("Pro-Visualize: Your complete proteomics visualization suite")
with tab_qc:
    # The QC module now handles its own data loading and rendering
    st.info("This section is under development.")
    st.header("Quality Control")
    #render_qc()

with tab_quant:
    st.header("Protein Quantification")
    st.info("This section is under development.")
    # Future: Call render_quant()

with tab_comp:
    st.header("Comparative Analysis")
    st.info("This section is under development.")
    # Future: Call render_comp()