# pro-visualize/utils/data_manager.py

import pandas as pd
import streamlit as st

class DataManager:
    """Handles loading and basic processing of uploaded data."""

    @staticmethod
    def load_protein_data(uploaded_file):
        """Loads protein-level data from a tsv or csv file."""
        if uploaded_file is None:
            return None
        try:
            # Assuming tab-separated, but you can add logic to detect separator
            return pd.read_csv(uploaded_file, sep='\t')
        except Exception as e:
            st.error(f"Error loading protein data: {e}")
            return None

    @staticmethod
    def load_annotation_data(uploaded_file):
        """Loads annotation data from a tsv or csv file."""
        if uploaded_file is None:
            return None
        try:
            return pd.read_csv(uploaded_file, sep='\t')
        except Exception as e:
            st.error(f"Error loading annotation data: {e}")
            return None