# pro_visualize/modules/qc_module.py

import streamlit as st

# Import the new class for the Targeted QC tab
from modules.qc_tabs.targeted_qc_tab import TargetedQcTab 
from modules.qc_tabs.dia_qc_tab import DiaQcTab

# You would keep your imports and UI function for the DIA tab
# from modules.dia_qc_tab import DiaQcTab # (If you also make one for DIA)

# Dummy placeholder for the DIA tab UI function
def render_dia_qc_ui():
    #st.info("This is the placeholder for the DIA QC tab.")
    dia_qc_page = DiaQcTab()
    dia_qc_page.render()


def render():
    """Renders the main QC page with tabs for different analysis types."""
    st.header("Quality Control Analysis")

    # --- Main QC Tabs ---
    dia_tab, targeted_tab = st.tabs(["DIA QC", "Targeted QC"])

    with dia_tab:
        render_dia_qc_ui() # Call the DIA QC render function

    with targeted_tab:
        # Instantiate and render the Targeted QC tab
        targeted_qc_page = TargetedQcTab()
        targeted_qc_page.render()