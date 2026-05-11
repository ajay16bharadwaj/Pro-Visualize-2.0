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
from utils.report_builder import ReportBuilder

# --- Page Configuration ---
st.set_page_config(
    page_title="Pro-Visualize",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

st.title("Pro-Visualize: Proteomics Data Visualization")
st.markdown("A scalable application for multi-omics analysis, starting with targeted and DIA-MS proteomics.")

# Ensure the report builder is always available before any tab renders
st.session_state.setdefault("report", ReportBuilder())

# --- Create Main Tabs ---
tab_welcome, tab_qc, tab_dilution, tab_quant, tab_comp, tab_scp, tab_report, tab_chat = st.tabs([
    "👋 Welcome",
    "📊 QC Analysis",
    "📈 Dilution Series",
    "📈 Quantification",
    "🆚 Comparative Analysis",
    "🔬 Single-Cell Proteomics",
    "📋 Report",
    "💬 Pro-Viz Chat"
])

with tab_welcome:
    st.title("Pro-Visualize 2.0")
    st.caption("Proteomics data visualization — from raw QC to single-cell DE, all in one place.")

    st.markdown("---")
    st.subheader("Quick-Start — pick your data type")

    c1, c2, c3 = st.columns(3)
    with c1:
        with st.container(border=True):
            st.markdown("### 📊 QC Analysis")
            st.markdown(
                "Have a **DIA-NN `.parquet` report**? Upload it to inspect "
                "RT/IM control charts, mass accuracy, peak width distributions, "
                "and sentinel peptide stability across runs."
            )
            st.markdown("*Also supports targeted (PRM/SRM) QC.*")
            st.info("→ Go to **QC Analysis** tab", icon="👆")

        with st.container(border=True):
            st.markdown("### 📈 Dilution Series")
            st.markdown(
                "Have a **loading-curve protein matrix** (e.g., from DIA-NN)? "
                "Visualize linearity (R²), CV across concentrations, LOD/LOQ "
                "estimates, and protein completeness."
            )
            st.info("→ Go to **Dilution Series** tab", icon="👆")

    with c2:
        with st.container(border=True):
            st.markdown("### 📈 Quantification")
            st.markdown(
                "Have a **protein-level intensity matrix** with a sample annotation "
                "file? Explore protein counts, PCA, sample correlation, missing-value "
                "patterns, Venn/UpSet overlaps, and rank-order plots."
            )
            st.info("→ Go to **Quantification** tab", icon="👆")

        with st.container(border=True):
            st.markdown("### 🆚 Comparative Analysis")
            st.markdown(
                "Have fold-change + FDR results from a statistical test? Upload your "
                "comparison table alongside the intensity matrix to generate volcano "
                "plots, expression heatmaps, violin plots, and pathway enrichment."
            )
            st.info("→ Go to **Comparative Analysis** tab", icon="👆")

    with c3:
        with st.container(border=True):
            st.markdown("### 🔬 Single-Cell Proteomics")
            st.markdown(
                "Have a **DIA-NN PG-matrix** from an SCP experiment? Run the full "
                "pipeline: QC filtering → normalization → PCA → UMAP → Leiden "
                "clustering → Wilcoxon DE → pathway enrichment → activity scoring."
            )
            st.info("→ Go to **Single-Cell Proteomics** tab", icon="👆")

        with st.container(border=True):
            st.markdown("### 📋 Report Builder")
            st.markdown(
                "Every plot in every module has an **Add to Report** button. "
                "Collect your key figures, add interpretation notes, then export "
                "an interactive HTML report or a ZIP bundle of PNG/SVG/HTML assets."
            )
            st.info("→ Go to **Report** tab", icon="👆")

    st.markdown("---")

    with st.expander("Expected file formats"):
        st.markdown("""
| Module | File 1 | File 2 (optional) |
|--------|--------|-------------------|
| QC (DIA) | DIA-NN report `.parquet` | — |
| QC (Targeted) | Skyline/Spectronaut report `.tsv/.csv` | — |
| Dilution Series | Protein matrix `.tsv/.csv` | Sample annotation `.tsv/.csv` |
| Quantification | Protein matrix `.tsv/.csv` | Sample annotation `.tsv/.csv` |
| Comparative | Protein matrix `.tsv/.csv` | Annotation + comparison table |
| SCP | DIA-NN PG-matrix `.tsv/.csv` | Sample annotation `.tsv/.csv` |

All text files are auto-detected (tab or comma separated).
Column names are **configurable** per module — you are not locked in to specific headers.
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

with tab_report:
    safe_render(
        "Report",
        st.session_state["report"].render_preview,
        reset_keys=["report"],
    )

with tab_chat:
    st.info("This section is currently under development. 🏗️")