# In visualizations/pro_viz_rag_system.py

# ============================================
# IMPORTS & SETUP
# ============================================

import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import pickle
import logging
import streamlit as st # Added for caching
import shutil
import tarfile

# LlamaIndex imports
from llama_index.core import Document, VectorStoreIndex, StorageContext, load_index_from_storage, Settings
from llama_index.llms.google_genai import GoogleGenAI
from llama_index.embeddings.google_genai import GoogleGenAIEmbedding
from llama_index.core.node_parser import SentenceSplitter

# Import from our existing codebase
from visualizations.comparative_visualizer import ComparativeVisualizer

# Suppress logs
logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("google_genai").setLevel(logging.WARNING)

# A simple tokenizer that just splits by space.
# This avoids the tiktoken dependency entirely.
def simple_tokenizer(text: str) -> list[str]:
    """A simple string split tokenizer."""
    return text.split()

# ============================================
# API KEY CONFIGURATION
# ============================================

# We will set these to None initially.
# They will be configured by the UI.
Settings.llm = None
Settings.embed_model = None

os.environ["GOOGLE_API_KEY"] = ""

# --- CONFIGURE LLAMA-INDEX SETTINGS ---
try:
    # Configure LLM (using your original model)
    Settings.llm = GoogleGenAI(
        model="models/gemini-2.0-flash",
        temperature=0.1
    )

    # Configure Embeddings
    Settings.embed_model = GoogleGenAIEmbedding(
        model="models/embedding-001",
    )
    
    # Configure the node parser to use the simple tokenizer
    Settings.node_parser = SentenceSplitter(tokenizer=simple_tokenizer)
    
    logging.info("✓ Google GenAI configured successfully!")

except Exception as e:
    logging.error(f"Failed to configure Google GenAI: {e}")
    st.error(f"CRITICAL ERROR: Failed to configure Google GenAI. Check your API key. Error: {e}")
# ----------------------------------------


# Suppress logs
logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("google_genai").setLevel(logging.WARNING)

# def configure_google_genai(api_key: str):
#     """
#     Sets the Google API key and configures the LLM
#     and embedding models for LlamaIndex.
#     """
#     try:
#         os.environ["GOOGLE_API_KEY"] = api_key

#         # Configure LLM
#         Settings.llm = GoogleGenAI(
#             model="models/gemini-2.0-flash",  
#             temperature=0.1
#         )

#         # Configure Embeddings
#         Settings.embed_model = GoogleGenAIEmbedding(
#             model="models/embedding-001",
#         )

#         # Explicitly set the node parser to use our
#         # simple tokenizer instead of the tiktoken default.
#         #Settings.node_parser = SentenceSplitter(tokenizer=simple_tokenizer)
        
#         logging.info("✓ Google GenAI configured successfully!")
#         return True
#     except Exception as e:
#         logging.error(f"Failed to configure Google GenAI: {e}")
#         st.error(f"Failed to configure Google GenAI: {e}. Please check your API key.")
#         return False


# ============================================
# DATASET INDEX MANAGER (User's Code)
# ============================================
class DatasetIndexManager:
    """
    Manages multiple dataset indexes with save/load capabilities.
    (This is your class, slightly adapted for Streamlit)
    """
    
    def __init__(self, base_dir="./dataset_indexes"):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(exist_ok=True)
        self.active_dataset = None
        self.active_index = None
        self.active_metadata = None
        self.active_comp_lookup = None
        self.active_enrichment_lookup = None
    
    def create_dataset_index(
        self, 
        dataset_name: str,
        protein_df,
        annotation_df,
        comparative_df,
        comp_visualizer: ComparativeVisualizer,
        column_config: dict,
        description: str = "",
        storage_format: str = 'parquet'  # <-- ADD THIS ARGUMENT
    ):
        """
        Create and save a complete dataset index.
        (This is your method)
        """
        logging.info(f"Creating index for: {dataset_name}")
        
        dataset_dir = self.base_dir / dataset_name
        dataset_dir.mkdir(exist_ok=True)
        
        # Save metadata (first pass)
        metadata = {
            'dataset_name': dataset_name,
            'description': description,
            'created_at': datetime.now().isoformat(),
            'protein_count': len(protein_df),
            'sample_count': len(annotation_df),
            'comparison_count': comparative_df[column_config['comparison_label']].nunique(),
            'comparisons': comparative_df[column_config['comparison_label']].unique().tolist(),
            'column_config': column_config
        }
        
        with open(dataset_dir / 'metadata.json', 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logging.info("✓ Metadata saved")
        
        # --- START REPLACEMENT ---
        # Save raw data (using chosen format)
        if storage_format == 'parquet':
            protein_df.to_parquet(dataset_dir / 'protein_data.parquet')
            annotation_df.to_parquet(dataset_dir / 'annotation_data.parquet')
            comparative_df.to_parquet(dataset_dir / 'comparative_data.parquet')
            metadata['storage_format'] = 'parquet' # Record what we did
            logging.info("✓ Raw data saved in Parquet format")
        else: # Default to CSV
            protein_df.to_csv(dataset_dir / 'protein_data.csv', index=False)
            annotation_df.to_csv(dataset_dir / 'annotation_data.csv', index=False)
            comparative_df.to_csv(dataset_dir / 'comparative_data.csv', index=False)
            metadata['storage_format'] = 'csv' # Record what we did
            logging.info("✓ Raw data saved in CSV format")
        
        # Re-save metadata to include the storage_format
        with open(dataset_dir / 'metadata.json', 'w') as f:
            json.dump(metadata, f, indent=2)
        # --- END REPLACEMENT ---

        # Build RAG index
        logging.info("📊 Building RAG index...")
        
        # Run enrichment
        enrichment_results = self._run_enrichment(
            comp_visualizer, comparative_df, protein_df, column_config
        )
        
        # Build lookup tables
        comp_lookup, enrichment_lookup = self._build_lookups(
            comparative_df, enrichment_results, protein_df, column_config
        )
        
        # Save lookup tables
        with open(dataset_dir / 'comp_lookup.pkl', 'wb') as f:
            pickle.dump(comp_lookup, f)
        with open(dataset_dir / 'enrichment_lookup.pkl', 'wb') as f:
            pickle.dump(enrichment_lookup, f)
        
        logging.info("✓ Lookup tables saved")
        
        # Generate documents
        documents = self._generate_documents(
            protein_df, comp_lookup, enrichment_lookup
        )
        
        logging.info(f"🔨 Creating vector index ({len(documents)} documents)...")
        
        index = VectorStoreIndex.from_documents(documents, show_progress=True)
        
        # Save index
        index_dir = dataset_dir / "rag_index"
        index.storage_context.persist(persist_dir=str(index_dir))
        
        logging.info(f"✅ Index created and saved at {dataset_dir}")
        return dataset_dir
    
    def list_datasets(self):
        """List all available dataset indexes."""
        datasets = []
        for dataset_dir in self.base_dir.iterdir():
            if dataset_dir.is_dir():
                metadata_file = dataset_dir / 'metadata.json'
                if metadata_file.exists():
                    with open(metadata_file) as f:
                        metadata = json.load(f)
                    datasets.append(metadata)
        return datasets
    
    # Use Streamlit's caching for loading data
    @st.cache_resource(show_spinner="Loading dataset index...")
    def load_dataset(_self, dataset_name: str):
        """
        Load a saved dataset index.
        Wrapped with st.cache_resource for performance.
        """
        dataset_dir = _self.base_dir / dataset_name
        
        if not dataset_dir.exists():
            raise ValueError(f"Dataset '{dataset_name}' not found")
        
        logging.info(f"📂 Loading dataset: {dataset_name}")
        
        with open(dataset_dir / 'metadata.json') as f:
            metadata = json.load(f)
        
        index_dir = dataset_dir / "rag_index"
        storage_context = StorageContext.from_defaults(persist_dir=str(index_dir))
        index = load_index_from_storage(storage_context)
        
        with open(dataset_dir / 'comp_lookup.pkl', 'rb') as f:
            comp_lookup = pickle.load(f)
        with open(dataset_dir / 'enrichment_lookup.pkl', 'rb') as f:
            enrichment_lookup = pickle.load(f)
        
        logging.info(f"✅ Loaded dataset: {dataset_name}")
        
        # Return loaded components
        return index, metadata, comp_lookup, enrichment_lookup
    
    def set_active_dataset(self, dataset_name: str):
        """Loads and activates a dataset for querying."""
        try:
            index, metadata, comp_lookup, enrichment_lookup = self.load_dataset(dataset_name)
            self.active_dataset = dataset_name
            self.active_index = index
            self.active_metadata = metadata
            self.active_comp_lookup = comp_lookup
            self.active_enrichment_lookup = enrichment_lookup
            logging.info(f"✅ Activated dataset: {dataset_name}")
            return True
        except Exception as e:
            st.error(f"Failed to activate dataset {dataset_name}: {e}")
            logging.error(f"Failed to activate dataset: {e}", exc_info=True)
            return False

    def export_dataset(self, dataset_name: str, export_path: str):
        """Export a dataset index as a portable archive."""
        dataset_dir = self.base_dir / dataset_name
        if not dataset_dir.exists():
            raise ValueError(f"Dataset '{dataset_name}' not found")
        logging.info(f"📦 Exporting dataset: {dataset_name}")
        shutil.make_archive(
            export_path.replace('.tar.gz', ''), 'gztar', dataset_dir
        )
        logging.info(f"✅ Exported to: {export_path}")
    
    def import_dataset(self, archive_path: str):
        """Import a dataset from an archive."""
        logging.info(f"📥 Importing dataset from: {archive_path}")
        with tarfile.open(archive_path, 'r:gz') as tar:
            tar.extractall(self.base_dir)
        logging.info(f"✅ Dataset imported successfully")
    
    # Helper methods
    def _run_enrichment(self, visualizer, comparative_df, protein_df, column_config):
        """Run enrichment analysis."""
        enrichment_results = {}
        comp_label = column_config['comparison_label']
        prot_id_col = column_config['protein_id']

        for comparison in comparative_df[comp_label].unique():
            logging.info(f"  Processing enrichment: {comparison}")
            
            sig_df = visualizer.filter_significant_proteins(
                comparison=comparison, fdr_cutoff=0.05, fc_cutoff=1.0
            )
            
            if len(sig_df) > 0:
                gene_list = sig_df.merge(
                    protein_df[[prot_id_col, 'Gene Name']], 
                    left_on=column_config['comp_protein_id'], 
                    right_on=prot_id_col
                )['Gene Name'].dropna().unique().tolist()
                
                logging.info(f"    Running enrichment on {len(gene_list)} genes...")
                
                try:
                    # Assuming 'mouse' for now, this could be a UI option
                    enrichment_df = visualizer.run_enrichment_analysis(
                        gene_list=gene_list, organism="mouse" 
                    )
                    enrichment_results[comparison] = enrichment_df
                except Exception as e:
                    logging.warning(f"    ⚠ Enrichment failed: {e}")
                    enrichment_results[comparison] = pd.DataFrame()
            else:
                enrichment_results[comparison] = pd.DataFrame()
        return enrichment_results
    
    def _build_lookups(self, comparative_df, enrichment_results, protein_df, column_config):
        """Build lookup dictionaries."""
        comp_label = column_config['comparison_label']
        comp_prot_id = column_config['comp_protein_id']
        fc_col = column_config['fold_change']
        fdr_col = column_config['fdr']

        # comp_lookup
        comp_lookup = {}
        for comparison in comparative_df[comp_label].unique():
            comp_data = comparative_df[comparative_df[comp_label] == comparison]
            for _, row in comp_data.iterrows():
                protein = row[comp_prot_id]
                if protein not in comp_lookup:
                    comp_lookup[protein] = {}
                comp_lookup[protein][comparison] = {
                    'log2FC': row[fc_col],
                    'FDR': row[fdr_col],
                    'significant': (row[fdr_col] <= 0.05) and (abs(row[fc_col]) >= 1.0),
                    'regulation': 'up' if row[fc_col] > 0 else 'down'
                }
        
        # enrichment_lookup (unchanged)
        enrichment_lookup = {}
        for comparison, enrich_df in enrichment_results.items():
            if not enrich_df.empty:
                for _, row in enrich_df.iterrows():
                    genes = row['genes'].split(';')
                    pathway_info = f"{row['source']}: {row['name']} (p={row['p_value']:.2e})"
                    for gene in genes:
                        if gene not in enrichment_lookup:
                            enrichment_lookup[gene] = []
                        enrichment_lookup[gene].append(pathway_info)
        return comp_lookup, enrichment_lookup
    
    def _generate_documents(self, protein_df, comp_lookup, enrichment_lookup):
        """Generate document objects."""
        documents = []
        logging.info("  Generating documents...")
        
        for idx, row in protein_df.iterrows():
            protein_id = row['Protein']
            gene_name = row.get('Gene Name', 'Unknown')
            comp_data = comp_lookup.get(protein_id, {})
            
            comp_info_lines = []
            if comp_data:
                for comparison, data in comp_data.items():
                    sig_marker = "**SIGNIFICANT**" if data['significant'] else ""
                    comp_info_lines.append(
                        f"{comparison}: {data['regulation']}-regulated "
                        f"(log2FC={data['log2FC']:.2f}, FDR={data['FDR']:.2e}) {sig_marker}"
                    )
                comp_info_text = "\n".join(comp_info_lines)
            else:
                comp_info_text = "No differential expression data available"
            
            enrich_info_text = "No pathway enrichment data available"
            if gene_name in enrichment_lookup:
                pathways = enrichment_lookup[gene_name][:5]
                enrich_info_text = "\n".join(pathways)
            
            doc_text = f"""
PROTEIN: {protein_id}
GENE: {gene_name}
DESCRIPTION:
{row.get('Protein Description', 'No description available')}
FUNCTIONAL ANNOTATIONS:
Biological Process: {row.get('Gene ontology (biological process)', 'N/A')}
Cellular Component: {row.get('Gene ontology (cellular component)', 'N/A')}
Molecular Function: {row.get('Gene ontology (molecular function)', 'N/A')}
LOCALIZATION:
Subcellular Location: {row.get('Subcellular Location[CC]', 'N/A')}
Tissue Specificity: {row.get('Tissue Specificity', 'N/A')}
EXPRESSION DATA:
{comp_info_text}
PATHWAY ASSOCIATIONS:
{enrich_info_text}
"""
            metadata = {
                'protein_id': protein_id, 'gene_name': gene_name,
                'document_type': 'protein',
            }
            if comp_data:
                for comp_label, comp_info in comp_data.items():
                    metadata[f"{comp_label}_log2fc"] = comp_info['log2FC']
                    metadata[f"{comp_label}_fdr"] = comp_info['FDR']
                    metadata[f"{comp_label}_significant"] = comp_info['significant']
            
            documents.append(Document(text=doc_text, metadata=metadata))
        
        logging.info(f"  ✓ Generated {len(documents)} documents")
        return documents

# ============================================
# LITERATURE INDEX MANAGER (User's Code)
# ============================================
class LiteratureIndexManager:
    """
    Manages pre-indexed scientific literature.
    (This is your class, with minor adaptations)
    """
    def __init__(self, literature_dir="./literature_index"):
        self.literature_dir = Path(literature_dir)
        self.literature_dir.mkdir(exist_ok=True)
        self.literature_index = None
    
    def build_literature_index(self, papers_directory: str):
        """Build literature index with LlamaIndex."""
        from llama_index.core import SimpleDirectoryReader
        
        logging.info("📚 Building literature index with LlamaIndex...")
        
        reader = SimpleDirectoryReader(
            input_dir=papers_directory,
            required_exts=[".pdf"],
            filename_as_id=True
        )
        documents = reader.load_data()
        
        logging.info(f"✓ Loaded {len(documents)} document chunks from PDFs")
        
        for doc in documents:
            doc.metadata['source'] = 'lab_publication'
            doc.metadata['document_type'] = 'research_paper'
        
        index = VectorStoreIndex.from_documents(documents, show_progress=True)
        
        index_dir = self.literature_dir / "llamaindex_index"
        index.storage_context.persist(persist_dir=str(index_dir))
        
        metadata = {
            'created_at': datetime.now().isoformat(),
            'paper_count': len(set([d.metadata.get('file_name', '') for d in documents])),
            'chunk_count': len(documents), 'index_type': 'llamaindex'
        }
        with open(self.literature_dir / 'metadata.json', 'w') as f:
            json.dump(metadata, f, indent=2)
        
        self.literature_index = index
        logging.info("✅ Literature index created with LlamaIndex!")
        return index
    
    @st.cache_resource(show_spinner="Loading literature index...")
    def load_literature_index(_self):
        """Load pre-built literature index."""
        logging.info("📚 Loading literature index...")
        metadata_file = _self.literature_dir / 'metadata.json'
        if not metadata_file.exists():
            raise FileNotFoundError(
                f"No literature index found at {_self.literature_dir}. "
            )
        
        with open(metadata_file) as f:
            metadata = json.load(f)
        
        index_dir = _self.literature_dir / "llamaindex_index"
        storage_context = StorageContext.from_defaults(persist_dir=str(index_dir))
        _self.literature_index = load_index_from_storage(storage_context)
        
        logging.info("✅ Loaded literature index (llamaindex)")
        return _self.literature_index
    
    def search_literature(self, query: str, top_k: int = 5):
        """Search literature with citations."""
        if self.literature_index is None:
            raise ValueError("Literature index not loaded")
        
        retriever = self.literature_index.as_retriever(similarity_top_k=top_k)
        nodes = retriever.retrieve(query)
        
        formatted_results = []
        for node in nodes:
            formatted_results.append({
                'text': node.node.text,
                'source': node.node.metadata.get('file_name', 'Unknown'),
                'page': node.node.metadata.get('page_label', 'N/A'),
                'score': node.score,
                'citation': f"{node.node.metadata.get('file_name', 'Unknown')}, p.{node.node.metadata.get('page_label', 'N/A')}"
            })
        return formatted_results

# ============================================
# ENHANCED QUANTITATIVE HANDLER (User's Code)
# ============================================
class QuantitativeHandler:
    """
    Handles quantitative queries using direct data operations.
    (This is your class)
    """
    def __init__(self, dataset_manager: DatasetIndexManager):
        self.dataset_manager = dataset_manager
    
    def handle(self, query: str, comparison: str = None):
        """Process quantitative queries."""
        query_lower = query.lower()
        comp_lookup = self.dataset_manager.active_comp_lookup
        protein_df = self._load_protein_df() # Helper to load
        
        available_comparisons = self.dataset_manager.active_metadata['comparisons']

        if not comparison:
            for comp in available_comparisons:
                if comp.lower() in query_lower:
                    comparison = comp
                    break
        
        if self._is_exclusive_query(query_lower):
            return self._handle_exclusive_comparison(query, available_comparisons, comp_lookup, protein_df)

        import re
        if 'top' in query_lower or 'most' in query_lower or 'highest' in query_lower or 'lowest' in query_lower:
            top_match = re.search(r'top\s+(\d+)', query_lower)
            most_match = re.search(r'(\d+)\s+most', query_lower)
            n = int(top_match.group(1)) if top_match else (int(most_match.group(1)) if most_match else 10)
            return self._handle_top_proteins(query, comparison, n, comp_lookup, protein_df)

        if 'how many' in query_lower or 'count' in query_lower or 'number of' in query_lower:
            return self._handle_count(query, comparison, comp_lookup, protein_df)

        elif 'percentage' in query_lower or '%' in query_lower:
            return self._handle_percentage(query, comparison, comp_lookup, protein_df)

        elif 'list' in query_lower or 'show me' in query_lower or 'what are' in query_lower:
            return self._handle_list(query, comparison, comp_lookup, protein_df)
        else:
            return None
    
    def _is_exclusive_query(self, query_lower):
        exclusive_patterns = [
            'in one condition and not', 'in one comparison and not', 'only in',
            'unique to', 'exclusively in', 'specific to', 'not in the other',
            'but not in'
        ]
        return any(pattern in query_lower for pattern in exclusive_patterns)
    
    def _load_protein_df(self):
        """
        Load protein dataframe from active dataset.
        Smart-loads either parquet or csv.
        """
        dataset_dir = self.dataset_manager.base_dir / self.dataset_manager.active_dataset
        
        parquet_file = dataset_dir / 'protein_data.parquet'
        csv_file = dataset_dir / 'protein_data.csv'
        
        if parquet_file.exists():
            logging.info("Loading protein data from Parquet file.")
            return pd.read_parquet(parquet_file)
        elif csv_file.exists():
            logging.info("Loading protein data from CSV file.")
            return pd.read_csv(csv_file)
        else:
            raise FileNotFoundError(f"No protein_data.parquet or protein_data.csv found in {dataset_dir}")
    
    def _get_gene_name(self, protein_id, protein_df):
        """Get gene name for a protein ID."""
        match = protein_df[protein_df['Protein'] == protein_id]
        if len(match) > 0:
            return match.iloc[0].get('Gene Name', 'Unknown')
        return 'Unknown'
    
    # (All other _handle_... methods from your script go here, unchanged)
    # ...
    def _handle_exclusive_comparison(self, query, available_comparisons, comp_lookup, protein_df):
        """
        Handle queries like:
        - "What proteins are upregulated in Control/Ethanol but not in Control/NES?"
        - "Proteins unique to Ethanol/NES"
        - "Show me proteins only in Control/Ethanol"
        """
        query_lower = query.lower()
        
        # Parse which comparisons are involved
        mentioned_comparisons = []
        for comp in available_comparisons:
            if comp.lower() in query_lower:
                mentioned_comparisons.append(comp)
        
        if len(mentioned_comparisons) < 2:
            # Try to find target "only in X"
            only_in_match = re.search(r'(only in|unique to|exclusively in|specific to)\s+([\w\/\-\_]+)', query_lower)
            if only_in_match and len(mentioned_comparisons) == 1:
                # This is a valid "only in X" query
                target_comparison = mentioned_comparisons[0]
                exclusion_comparisons = [c for c in available_comparisons if c != target_comparison]
            else:
                return (
                    f"For exclusive comparison queries, please mention at least 2 comparisons "
                    f"or use a format like 'only in {available_comparisons[0]}'. "
                    f"Available: {', '.join(available_comparisons)}"
                )
        else:
            # Determine which is the "target" and which are "exclusions"
            # Pattern: "in X but not in Y" or "only in X"
            target_comparison = mentioned_comparisons[0]
            exclusion_comparisons = mentioned_comparisons[1:]
        
        # Detect regulation
        regulation = None
        if 'upregulated' in query_lower or 'up-regulated' in query_lower:
            regulation = 'up'
        elif 'downregulated' in query_lower or 'down-regulated' in query_lower:
            regulation = 'down'
        
        significant_only = 'significant' in query_lower
        
        # Find exclusive proteins
        exclusive_proteins = []
        
        for protein_id, comparisons_data in comp_lookup.items():
            # Check if protein is in target comparison
            if target_comparison not in comparisons_data:
                continue
            
            target_data = comparisons_data[target_comparison]
            
            # Check significance
            if significant_only and not target_data['significant']:
                continue
            
            # Check regulation
            if regulation:
                if regulation == 'up' and target_data['log2FC'] <= 0:
                    continue
                if regulation == 'down' and target_data['log2FC'] >= 0:
                    continue
            
            # Check if NOT significant in exclusion comparisons
            is_exclusive = True
            for excl_comp in exclusion_comparisons:
                if excl_comp in comparisons_data:
                    excl_data = comparisons_data[excl_comp]
                    # If significant in exclusion comparison, not exclusive
                    if excl_data['significant']:
                        is_exclusive = False
                        break
            
            if is_exclusive:
                gene_name = self._get_gene_name(protein_id, protein_df)
                exclusive_proteins.append({
                    'protein_id': protein_id,
                    'gene_name': gene_name,
                    'log2fc': target_data['log2FC'],
                    'fdr': target_data['FDR']
                })
        
        # Format response
        desc = "significant " if significant_only else ""
        reg_desc = f"{regulation}regulated " if regulation else ""
        
        response = f"**Proteins that are {desc}{reg_desc}in {target_comparison} "
        response += f"but NOT significant in {', '.join(exclusion_comparisons)}:**\n\n"
        response += f"Found **{len(exclusive_proteins)}** exclusive proteins:\n\n"
        
        if len(exclusive_proteins) > 0:
            # Sort by fold change
            exclusive_proteins.sort(key=lambda x: abs(x['log2fc']), reverse=True)
            
            # Show top 20
            for i, p in enumerate(exclusive_proteins[:20], 1):
                response += f"{i}. **{p['gene_name']}** ({p['protein_id']}): "
                response += f"log2FC={p['log2fc']:.2f}, FDR={p['fdr']:.2e}\n"
            
            if len(exclusive_proteins) > 20:
                response += f"\n... and {len(exclusive_proteins) - 20} more."
        else:
            response += "*No proteins found matching these criteria.*"
        
        return response
    
    def _handle_count(self, query, comparison, comp_lookup, protein_df):
        """Count proteins matching criteria - with gene names."""
        if not comparison:
            available = self.dataset_manager.active_metadata['comparisons']
            return f"Please specify a comparison in your query. Available: {', '.join(available)}"
        
        query_lower = query.lower()
        significant_only = 'significant' in query_lower
        regulation = None
        if 'upregulated' in query_lower or 'up-regulated' in query_lower:
            regulation = 'up'
        elif 'downregulated' in query_lower or 'down-regulated' in query_lower:
            regulation = 'down'
        
        count = 0
        matching_proteins = []
        
        for protein_id, data in comp_lookup.items():
            if comparison not in data:
                continue
            
            comp_data = data[comparison]
            
            if significant_only and not comp_data['significant']:
                continue
            
            if regulation:
                if regulation == 'up' and comp_data['log2FC'] <= 0:
                    continue
                if regulation == 'down' and comp_data['log2FC'] >= 0:
                    continue
            
            count += 1
            gene_name = self._get_gene_name(protein_id, protein_df)
            matching_proteins.append({
                'protein_id': protein_id,
                'gene_name': gene_name,
                'log2fc': comp_data['log2FC'],
                'fdr': comp_data['FDR']
            })
        
        desc = "significant " if significant_only else ""
        reg_desc = f"{regulation}regulated " if regulation else ""
        
        response = f"There are **{count}** {desc}{reg_desc}proteins in the {comparison} comparison."
        
        if count > 0 and count <= 5:
            response += "\n\nProteins:\n"
            for p in matching_proteins:
                response += f"- **{p['gene_name']}** ({p['protein_id']}): log2FC={p['log2fc']:.2f}\n"
        elif count > 5:
            response += f"\n\nTop 5 examples:\n"
            sorted_proteins = sorted(matching_proteins, key=lambda x: abs(x['log2fc']), reverse=True)
            for p in sorted_proteins[:5]:
                response += f"- **{p['gene_name']}** ({p['protein_id']}): log2FC={p['log2fc']:.2f}\n"
        
        return response
    
    def _handle_list(self, query, comparison, comp_lookup, protein_df):
        """
        List proteins matching criteria - with gene names.
        Handles queries like "List all upregulated proteins" or "Show me downregulated proteins"
        """
        if not comparison:
            available = self.dataset_manager.active_metadata['comparisons']
            return f"Please specify a comparison in your query. Available: {', '.join(available)}"
        
        query_lower = query.lower()
        significant_only = 'significant' in query_lower
        regulation = None
        if 'upregulated' in query_lower or 'up-regulated' in query_lower:
            regulation = 'up'
        elif 'downregulated' in query_lower or 'down-regulated' in query_lower:
            regulation = 'down'
        
        matching_proteins = []
        
        for protein_id, data in comp_lookup.items():
            if comparison not in data:
                continue
            
            comp_data = data[comparison]
            
            if significant_only and not comp_data['significant']:
                continue
            
            if regulation:
                if regulation == 'up' and comp_data['log2FC'] <= 0:
                    continue
                if regulation == 'down' and comp_data['log2FC'] >= 0:
                    continue
            
            gene_name = self._get_gene_name(protein_id, protein_df)
            matching_proteins.append({
                'protein_id': protein_id,
                'gene_name': gene_name,
                'log2fc': comp_data['log2FC'],
                'fdr': comp_data['FDR'],
                'significant': comp_data['significant']
            })
        
        matching_proteins.sort(key=lambda x: abs(x['log2fc']), reverse=True)
        
        desc = "significant " if significant_only else ""
        reg_desc = f"{regulation}regulated " if regulation else ""
        
        response = f"**{desc.capitalize()}{reg_desc}proteins in {comparison}:**\n\n"
        response += f"Found **{len(matching_proteins)}** proteins.\n\n"
        
        max_show = min(50, len(matching_proteins))
        
        for i, p in enumerate(matching_proteins[:max_show], 1):
            sig_marker = "✓" if p['significant'] else ""
            response += f"{i}. **{p['gene_name']}** ({p['protein_id']}): "
            response += f"log2FC={p['log2fc']:.2f}, FDR={p['fdr']:.2e} {sig_marker}\n"
        
        if len(matching_proteins) > max_show:
            response += f"\n... and {len(matching_proteins) - max_show} more proteins."
        
        return response
    
    def _handle_percentage(self, query, comparison, comp_lookup, protein_df):
        """Calculate percentage - with examples."""
        if not comparison:
            available = self.dataset_manager.active_metadata['comparisons']
            return f"Please specify a comparison in your query. Available: {', '.join(available)}"
        
        query_lower = query.lower()
        significant_only = 'significant' in query_lower
        regulation = None
        if 'upregulated' in query_lower or 'up-regulated' in query_lower:
            regulation = 'up'
        elif 'downregulated' in query_lower or 'down-regulated' in query_lower:
            regulation = 'down'
        
        total = len([p for p, d in comp_lookup.items() if comparison in d])
        count = 0
        examples = []
        
        for protein_id, data in comp_lookup.items():
            if comparison not in data:
                continue
            
            comp_data = data[comparison]
            
            if significant_only and not comp_data['significant']:
                continue
            
            if regulation:
                if regulation == 'up' and comp_data['log2FC'] <= 0:
                    continue
                if regulation == 'down' and comp_data['log2FC'] >= 0:
                    continue
            
            count += 1
            if len(examples) < 3:
                gene_name = self._get_gene_name(protein_id, protein_df)
                examples.append({
                    'gene_name': gene_name,
                    'protein_id': protein_id,
                    'log2fc': comp_data['log2FC']
                })
        
        percentage = (count / total * 100) if total > 0 else 0
        
        desc = "significant " if significant_only else ""
        reg_desc = f"{regulation}regulated " if regulation else ""
        
        response = f"**{percentage:.1f}%** ({count}/{total}) of proteins are {desc}{reg_desc}in the {comparison} comparison."
        
        if examples:
            response += "\n\nExamples:\n"
            for ex in examples:
                response += f"- **{ex['gene_name']}** ({ex['protein_id']}): log2FC={ex['log2fc']:.2f}\n"
        
        return response
    
    def _handle_top_proteins(self, query, comparison, n, comp_lookup, protein_df):
        """Get top N proteins - with gene names."""
        if not comparison:
            available = self.dataset_manager.active_metadata['comparisons']
            return f"Please specify a comparison in your query. Available: {', '.join(available)}"
        
        query_lower = query.lower()
        regulation = None
        if 'upregulated' in query_lower or 'up-regulated' in query_lower:
            regulation = 'up'
        elif 'downregulated' in query_lower or 'down-regulated' in query_lower:
            regulation = 'down'
        
        protein_data = []
        for protein_id, data in comp_lookup.items():
            if comparison not in data:
                continue
            
            comp_data = data[comparison]
            
            if regulation:
                if regulation == 'up' and comp_data['log2FC'] <= 0:
                    continue
                if regulation == 'down' and comp_data['log2FC'] >= 0:
                    continue
            
            gene_name = self._get_gene_name(protein_id, protein_df)
            protein_data.append({
                'protein_id': protein_id,
                'gene_name': gene_name,
                'log2FC': comp_data['log2FC'],
                'FDR': comp_data['FDR'],
                'significant': comp_data['significant']
            })
        
        protein_data.sort(key=lambda x: abs(x['log2FC']), reverse=True)
        top_proteins = protein_data[:n]
        
        reg_desc = f"{regulation}regulated " if regulation else ""
        response = f"**Top {n} {reg_desc}proteins in {comparison}:**\n\n"
        
        for i, p in enumerate(top_proteins, 1):
            sig = "✓" if p['significant'] else ""
            response += f"{i}. **{p['gene_name']}** ({p['protein_id']}): "
            response += f"log2FC={p['log2FC']:.2f}, FDR={p['FDR']:.2e} {sig}\n"
        
        return response

# ============================================
# UNIFIED QUERY SYSTEM (User's Code)
# ============================================
class UnifiedQuerySystem:
    """
    Query system that combines active dataset with literature knowledge.
    (This is your class)
    """
    def __init__(
        self,
        dataset_manager: DatasetIndexManager,
        literature_manager: LiteratureIndexManager
    ):
        self.dataset_manager = dataset_manager
        self.literature_manager = literature_manager
        self.quant_handler = QuantitativeHandler(dataset_manager)
    
    def query(self, question: str, use_literature: bool = True):
        """
        Unified query that routes to appropriate handler.
        """
        logging.info(f"Query: {question}")
        
        # Check if quantitative query
        query_lower = question.lower()
        
        # Use quant handler
        quant_response = self.quant_handler.handle(question)
        if quant_response:
            logging.info("Detected: Quantitative query")
            return {'answer': quant_response, 'type': 'quantitative'}
        
        # Semantic query with RAG
        logging.info("🔍 Searching dataset...")
        
        query_engine = self.dataset_manager.active_index.as_query_engine(
            similarity_top_k=10,
            response_mode="tree_summarize"
        )
        
        dataset_response = query_engine.query(question)
        
        dataset_proteins = self._extract_proteins_from_response(dataset_response)
        logging.info(f"   Found {len(dataset_proteins)} relevant proteins")
        
        if not use_literature or self.literature_manager.literature_index is None:
            return {
                'answer': dataset_response.response,
                'type': 'semantic',
                'proteins': dataset_proteins
            }
        
        # Search literature
        logging.info("📚 Searching literature...")
        literature_results = []
        for protein in dataset_proteins[:10]: # Limit to top 10 proteins
            query_str = f"{protein['gene_name']} {protein['protein_id']}"
            try:
                lit_results = self.literature_manager.search_literature(query_str, top_k=3)
                if lit_results:
                    literature_results.append({
                        'protein': protein, 'papers': lit_results
                    })
            except Exception as e:
                logging.warning(f"   ⚠ Error searching for {protein['gene_name']}: {e}")
        
        logging.info(f"   Found {len(literature_results)} proteins with literature mentions")
        
        combined_answer = self._synthesize_answer(
            dataset_response.response, literature_results
        )
        
        return {
            'answer': combined_answer, 'type': 'literature_integrated',
            'proteins': dataset_proteins, 'literature': literature_results
        }
    
    def cross_reference_proteins(self, comparison: str, regulation: str = None):
        """Find proteins in your data that appear in literature."""
        logging.info(f"Cross-referencing: {comparison}")
        comp_lookup = self.dataset_manager.active_comp_lookup
        
        proteins = []
        for protein_id, data in comp_lookup.items():
            if comparison not in data: continue
            comp_data = data[comparison]
            if not comp_data['significant']: continue
            
            if regulation:
                if regulation == 'up' and comp_data['log2FC'] <= 0: continue
                if regulation == 'down' and comp_data['log2FC'] >= 0: continue
            
            proteins.append({
                'protein_id': protein_id, 'log2fc': comp_data['log2FC'],
                'fdr': comp_data['FDR']
            })
        
        logging.info(f"🔍 Found {len(proteins)} proteins matching criteria")
        
        if self.literature_manager.literature_index is None:
            logging.warning("⚠ No literature index available")
            return pd.DataFrame(proteins)
        
        logging.info("📚 Searching literature...")
        results = []
        for protein in proteins:
            try:
                lit_results = self.literature_manager.search_literature(
                    f"{protein['protein_id']}", top_k=5
                )
                if lit_results:
                    results.append({
                        **protein, 'literature_mentions': len(lit_results),
                        'papers': [r['source'] for r in lit_results],
                        'citations': [r['citation'] for r in lit_results]
                    })
            except Exception as e:
                continue
        
        results_df = pd.DataFrame(results)
        logging.info(f"✅ Found {len(results_df)} proteins with literature support")
        return results_df
    
    def _extract_proteins_from_response(self, response):
        """Extract protein information from response."""
        proteins = []
        if hasattr(response, 'source_nodes'):
            for node in response.source_nodes:
                metadata = node.node.metadata
                if metadata.get('document_type') == 'protein':
                    proteins.append({
                        'protein_id': metadata.get('protein_id'),
                        'gene_name': metadata.get('gene_name')
                    })
        return proteins
    
    def _synthesize_answer(self, dataset_response, literature_results):
        """Synthesize final answer."""
        answer = f"{dataset_response}\n\n"
        if literature_results:
            answer += "**Literature Support:**\n\n"
            for result in literature_results:
                protein = result['protein']
                papers = result['papers']
                answer += f"• **{protein['gene_name']}** ({protein['protein_id']}): "
                answer += f"Mentioned in {len(papers)} paper(s)\n"
                for paper in papers[:2]:
                    answer += f"  - {paper['citation']}\n"
                answer += "\n"
        else:
            answer += "\n*No literature mentions found for these proteins.*"
        return answer

# ============================================
# COMPLETE PROTEOMICS RAG SYSTEM (User's Code)
# ============================================
class ProteomicsRAGSystem:
    """
    Complete system integrating datasets and literature.
    (This is your class, as the main entry point)
    """
    def __init__(self):
        self.dataset_manager = DatasetIndexManager()
        self.literature_manager = LiteratureIndexManager()
        self.query_system = None
    
    def setup_literature(self, papers_directory: str):
        """One-time setup: Index all lab papers."""
        logging.info("LITERATURE SETUP (One-time)")
        self.literature_manager.build_literature_index(papers_directory)
        logging.info("✅ Literature index ready for all future queries!")
    
    def index_new_dataset(
        self,
        dataset_name: str,
        protein_df,
        annotation_df,
        comparative_df,
        comp_visualizer,
        column_config: dict,
        description: str = "",
        storage_format: str = 'parquet' 
    ):
        """Index a new experiment dataset."""
        logging.info(f"INDEXING NEW DATASET: {dataset_name}")
        self.dataset_manager.create_dataset_index(
            dataset_name=dataset_name, protein_df=protein_df,
            annotation_df=annotation_df, comparative_df=comparative_df,
            comp_visualizer=comp_visualizer, column_config=column_config,
            description=description,
            storage_format=storage_format 
        )
        logging.info(f"✅ Dataset '{dataset_name}' indexed and saved!")
    
    def list_available_datasets(self):
        """Show all indexed datasets."""
        return self.dataset_manager.list_datasets()
    
    def activate_dataset(self, dataset_name: str):
        """Load a dataset for querying."""
        logging.info(f"ACTIVATING DATASET: {dataset_name}")
        
        success = self.dataset_manager.set_active_dataset(dataset_name)
        if not success:
            return
        
        if self.literature_manager.literature_index is None:
            try:
                self.literature_manager.load_literature_index()
            except FileNotFoundError:
                logging.warning("⚠ No literature index found. Literature features disabled.")
        
        self.query_system = UnifiedQuerySystem(
            self.dataset_manager, self.literature_manager
        )
        logging.info(f"✅ Dataset '{dataset_name}' is now active!")
    
    def query(self, question: str, use_literature: bool = True):
        """Ask questions combining dataset + literature."""
        if self.query_system is None:
            logging.warning("⚠ No active dataset. Activate one first!")
            st.warning("Please activate a dataset first.")
            return {'answer': 'Error: No active dataset.', 'type': 'error'}
        
        return self.query_system.query(question, use_literature=use_literature)
    
    def find_proteins_in_literature(self, comparison: str, regulation: str = None):
        """Cross-reference your DE proteins with literature."""
        if self.query_system is None:
            logging.warning("⚠ No active dataset. Activate one first!")
            return None
        
        return self.query_system.cross_reference_proteins(
            comparison=comparison, regulation=regulation
        )