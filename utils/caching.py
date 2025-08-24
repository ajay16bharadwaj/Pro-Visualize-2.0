# pro-visualize/utils/caching.py

import streamlit as st
from visualizations.targettedQCVisualization import TargetedQcVisualizer
from visualizations.DiaQcVisualizer import DiaQcVisualizer

@st.cache_data(show_spinner="Loading and processing Targeted QC data...")
def load_targeted_visualizer(filepath: str) -> TargetedQcVisualizer:
    """
    Loads, processes, and caches a TargetedQcVisualizer instance.
    Streamlit will only run this function if the filepath changes.
    """
    return TargetedQcVisualizer(filepath)

@st.cache_data(show_spinner="Loading and processing DIA-NN report...")
def load_dia_visualizer(filepath: str) -> DiaQcVisualizer:
    """
    Loads, processes, and caches a DiaQcVisualizer instance.
    Streamlit will only run this function if the filepath changes.
    """
    return DiaQcVisualizer(filepath)