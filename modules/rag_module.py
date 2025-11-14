# In modules/rag_module.py

import streamlit as st
import pandas as pd
import logging

# Import the new backend system
from visualizations.pro_viz_rag_system import ProteomicsRAGSystem
# Import existing visualizer to pass to the indexer
from visualizations.comparative_visualizer import ComparativeVisualizer

logger = logging.getLogger(__name__)

class RagTab:
    """
    Encapsulates the Streamlit UI for the Pro-Viz Chat (RAG) feature.
    (API key is now hardcoded in the backend)
    """

    def __init__(self):
        """Initializes the RagTab and its session state."""
        if 'rag_system' not in st.session_state:
            # Instantiate the main RAG system class
            st.session_state.rag_system = ProteomicsRAGSystem()
        if 'rag_chat_history' not in st.session_state:
            st.session_state.rag_chat_history = []
        
        # We no longer need the 'google_api_key_configured' state

    def _render_indexing_tab(self):
        # ... (code for st.subheader and rag_system)
        
        with st.expander("Upload & Configure Files", expanded=True):
            dataset_name = st.text_input("Unique Dataset Name", 
                                         help="e.g., 'Project_X_Phase_1_bCPA_Treated'")
            description = st.text_area("Description (Optional)",
                                       help="A brief description of this dataset.")
            
            # --- ADD THIS WIDGET ---
            storage_format_label = st.radio(
                "Select Index Storage Format",
                ["Parquet (Faster, Smaller)", "CSV (Legacy, Compatible)"],
                index=0,
                help="Parquet is recommended for new indexes. Use CSV if you need to match an older index format."
            )
            # -------------------------

            st.markdown("---")
            st.markdown("##### Upload Files")
            # ... (file uploaders) ...

            st.markdown("---")
            st.markdown("##### Configure Column Names")
            # ... (column config UI) ...

            if st.button("Create Dataset Index", use_container_width=True, type="primary"):
                if not all([dataset_name, protein_file, annotation_file, comparative_file]):
                    st.warning("Please fill in all required fields.")
                else:
                    # Map label to simple string key
                    storage_format_key = 'parquet' if "Parquet" in storage_format_label else 'csv'
                    
                    self._run_indexing(
                        dataset_name, description,
                        protein_file, annotation_file, comparative_file,
                        protein_id_col, sample_id_col, grouping_col,
                        comp_protein_id_col, fold_change_col, fdr_col, comparison_label_col,
                        storage_format_key  # <-- Pass the new key
                    )
    
    def _run_indexing(self, dataset_name, description, protein_file, annotation_file, comparative_file,
                      protein_id_col, sample_id_col, grouping_col,
                      comp_protein_id_col, fold_change_col, fdr_col, comparison_label_col,
                      storage_format: str):  # <-- ADD THIS ARGUMENT
        """Helper function to run the indexing process."""
        try:
            # ... (code to read dataframes and init visualizer) ...
            
            comp_visualizer = ComparativeVisualizer(
                protein_df, annotation_df, comparative_df, column_config
            )

            with st.spinner(f"Indexing '{dataset_name}'... This may take several minutes."):
                st.session_state.rag_system.index_new_dataset(
                    dataset_name=dataset_name,
                    protein_df=protein_df,
                    annotation_df=annotation_df,
                    comparative_df=comparative_df,
                    comp_visualizer=comp_visualizer,
                    column_config=column_config,
                    description=description,
                    storage_format=storage_format # <-- PASS IT HERE
                )
            
            st.success(f"Dataset '{dataset_name}' indexed successfully!")
            
        except Exception as e:
            st.error(f"An error occurred during indexing: {e}")
            logger.error(f"RAG indexing error: {e}", exc_info=True)


    def _render_chat_tab(self):
        """UI for selecting and chatting with an indexed dataset."""
        st.subheader("2. Chat with Your Dataset")
        rag_system = st.session_state.rag_system
        
        available_datasets = rag_system.list_available_datasets()
        dataset_names = [ds['dataset_name'] for ds in available_datasets]

        if not available_datasets:
            st.info("No datasets have been indexed yet. Please go to the 'Index New Dataset' tab first.")
            return

        selected_dataset = st.selectbox(
            "Select an indexed dataset to chat with:",
            options=dataset_names,
            key="rag_dataset_selector"
        )
        
        if st.button("Activate Dataset", use_container_width=True):
            if selected_dataset:
                with st.spinner(f"Activating '{selected_dataset}'..."):
                    rag_system.activate_dataset(selected_dataset)
                st.session_state.rag_chat_history = [] # Clear history
                st.success(f"Dataset '{selected_dataset}' is now active!")
                st.rerun()

        st.markdown("---")
        
        # Chat interface
        if rag_system.query_system:
            st.markdown(f"**Ready!** You are now chatting with `{rag_system.dataset_manager.active_dataset}`.")
            
            # --- 1. ADD THE TOGGLE ---
            show_details = st.toggle("Show Query Details ⚙️", key="rag_show_details")

            # --- 2. MODIFY HISTORY RENDERING ---
            for message in st.session_state.rag_chat_history:
                with st.chat_message(message["role"]):
                    
                    if message["role"] == "user":
                        # User message is just a string
                        st.markdown(message["content"])
                    
                    else: # Role is "assistant"
                        # Assistant message is now a dictionary
                        response_dict = message["content"]
                        
                        # A. Always display the main answer
                        st.markdown(response_dict.get('answer', 'Error: No answer found.'))
                        
                        # B. Conditionally display the details
                        if show_details:
                            # Create a clean metadata dict to show
                            details_to_show = response_dict.copy()
                            
                            # We don't need to show the 'answer' key again
                            details_to_show.pop('answer', None) 
                            
                            # Format the routing type for the title
                            route_type = details_to_show.pop('type', 'unknown').upper()
                            
                            with st.expander(f"**Query Route:** `{route_type}`", expanded=False):
                                if route_type == 'QUANTITATIVE':
                                    st.info("💡 This query was answered by the **Quantitative Handler** (the 'calculator')")
                                elif route_type == 'SEMANTIC':
                                    st.info("💡 This query was answered by the **Semantic RAG Engine** (the 'brain').")
                                elif route_type == 'LITERATURE_INTEGRATED':
                                    st.info("💡 This query was answered by the **Semantic RAG Engine** and then **cross-referenced with the Literature Index**.")
                                
                                # Show the rest of the metadata
                                st.json(details_to_show)

            # --- 3. MODIFY INPUT PROCESSING ---
            if prompt := st.chat_input("Ask a question about your data..."):
                # Add user message (this is fine, it's just a string)
                st.session_state.rag_chat_history.append({"role": "user", "content": prompt})
                with st.chat_message("user"):
                    st.markdown(prompt)
                
                with st.chat_message("assistant"):
                    with st.spinner("Thinking..."):
                        # Get the *full dictionary* from the backend
                        response_dict = rag_system.query(prompt)
                        
                        # Display just the answer part
                        response_text = response_dict.get('answer', 'Sorry, I encountered an error.')
                        st.markdown(response_text)
                
                # Store the *entire dictionary* in history, not just the text
                st.session_state.rag_chat_history.append({"role": "assistant", "content": response_dict})
        
        else:
            st.info("Please select and activate a dataset to begin chatting.")

    def render(self):
        """Renders the entire RAG tab."""
        st.header("Pro-Viz Chat 💬")
        st.info("The RAG system is configured with your API key and ready to index or chat.")

        tab1, tab2 = st.tabs(["Index New Dataset", "Chat with Dataset"])

        with tab1:
            self._render_indexing_tab()
            
        with tab2:
            self._render_chat_tab()

# ============================================
# ENTRY POINT
# ============================================
def render():
    """Entry point function to render the RagTab."""
    rag_tab = RagTab()
    rag_tab.render()