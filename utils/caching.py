# pro-visualize/utils/caching.py

import streamlit as st
from visualizations.targettedQCVisualization import TargetedQcVisualizer
from visualizations.DiaQcVisualizer import DiaQcVisualizer
from visualizations.comparative_visualizer import ComparativeVisualizer

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


@st.cache_data(show_spinner="Running pathway enrichment analysis (this may take a moment)...")
def run_cached_enrichment(_visualizer: ComparativeVisualizer, gene_list: list, organism: str):
    """
    Cached wrapper for the Enrichr API call to avoid re-running it.
    The `_visualizer` parameter has an underscore to tell Streamlit's caching
    mechanism not to hash the object itself, but to treat it as a stable component.
    """
    return _visualizer.run_enrichment_analysis(gene_list, organism)