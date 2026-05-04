import streamlit as st
import warnings

warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    message="Downcasting object dtype arrays on .fillna"
)
# --- Import Project Modules ---
from modules.qc_module import render as render_qc
from modules.dilution_module import render as render_dilution
from modules.quant_module import render as render_quant
from modules.comparative_module import render as render_comp
from modules.scp_module import render as render_scp
from utils.helpers import safe_render

# --- Page Configuration ---
st.set_page_config(
    page_title="Pro-Visualize",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

st.title("Pro-Visualize: Proteomics Data Visualization")
st.markdown("A scalable application for multi-omics analysis, starting with targeted and DIA-MS proteomics.")

# --- Create Main Tabs ---
tab_welcome, tab_qc, tab_dilution, tab_quant, tab_comp, tab_scp, tab_chat = st.tabs([
    "👋 Welcome",
    "📊 QC Analysis",
    "📈 Dilution Series",
    "📈 Quantification",
    "🆚 Comparative Analysis",
    "🔬 Single-Cell Proteomics",
    "💬 Pro-Viz Chat"
])

with tab_welcome:
    st.header("Welcome to Pro-Visualize!")
    st.markdown("""
    This application is designed to provide a comprehensive suite of tools for visualizing and analyzing proteomics data.

    **Current Modules:**
    - **QC Analysis**: Upload and analyze Quality Control data from both DIA and Targeted (e.g., PRM/SRM) experiments.
    - **Dilution Series (Loading-curve)**: Upload and analyze loading curve data from DIA experiments.
    - **Quantification**: Upload protein-level data to perform quantitative analysis.
    - **Comparative Analysis**: Visualize statistical comparisons between groups.
    - **Single-Cell Proteomics**: End-to-end SCP analysis — QC filtering, normalization, regression of technical covariates,
      PCA, UMAP, Leiden clustering, differential expression (Wilcoxon), and per-cell activity scoring.

    Navigate to the desired tab to begin.
    """)

with tab_qc:
    safe_render(
        "QC Analysis",
        render_qc,
        reset_keys=[
            "dia_qc_visualizer", "targeted_qc_visualizer",
            "dia_metadata_confirmed", "q_value_cutoff", "sentinel_peptides",
        ],
    )

with tab_dilution:
    safe_render(
        "Dilution Series",
        render_dilution,
        reset_keys=["dilution_visualizer"],
    )

with tab_quant:
    safe_render(
        "Quantification",
        render_quant,
        reset_keys=["quant_visualizer"],
    )

with tab_comp:
    safe_render(
        "Comparative Analysis",
        render_comp,
        reset_keys=[
            "comp_visualizer", "enrichment_results",
            "selected_comparison", "significant_proteins",
        ],
    )

with tab_scp:
    safe_render(
        "Single-Cell Proteomics",
        render_scp,
        reset_keys=[
            "scp_visualizer", "scp_de_results",
            "scp_enr_results", "scp_score_cols",
        ],
    )

with tab_chat:
    st.info("This section is currently under development. 🏗️")