import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
from io import BytesIO
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import requests
import time
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from config.plot_configs import HUMAN_TRANSCRIPTION_FACTORS
from utils.helpers import handle_plotting_errors

logger = logging.getLogger(__name__)

class ComparativeVisualizer:
    """
    Handles data loading, validation, and visualization for comparative analysis.
    """

    def __init__(self, protein_df: pd.DataFrame, annotation_df: pd.DataFrame, 
                 comparative_df: pd.DataFrame, column_config: dict):
        """
        Initializes the visualizer with all necessary data and configuration.
        """
        logger.info("Initializing ComparativeVisualizer...")
        self.column_config = column_config
        
        # Validate and store the dataframes
        self._validate_data(protein_df, annotation_df, comparative_df)
        self.protein_df = protein_df.copy()
        self.annotation_df = annotation_df.copy()
        self.comparative_df = comparative_df.copy()

        sample_id_col = self.column_config['sample_id']
        self.sample_cols = self.annotation_df[sample_id_col].unique().tolist()
        
        logger.info("ComparativeVisualizer initialized successfully.")

    def _validate_data(self, protein_df, annotation_df, comparative_df):
        """Performs initial validation of all input dataframes."""
        protein_col = self.column_config['protein_id']
        sample_col = self.column_config['sample_id']
        group_col = self.column_config['grouping']
        comp_protein_col = self.column_config['comp_protein_id']
        fc_col = self.column_config['fold_change']
        fdr_col = self.column_config['fdr']
        label_col = self.column_config['comparison_label']

        if protein_col not in protein_df.columns:
            raise ValueError(f"Protein data is missing the specified protein column: '{protein_col}'")
        if sample_col not in annotation_df.columns:
            raise ValueError(f"Annotation data is missing the sample column: '{sample_col}'")
        if group_col not in annotation_df.columns:
            raise ValueError(f"Annotation data is missing the grouping column: '{group_col}'")
        if any(c not in comparative_df.columns for c in [comp_protein_col, fc_col, fdr_col, label_col]):
            raise ValueError("Comparative data is missing one or more required columns.")

    def get_protein_data_preview(self):
        return self.protein_df

    def get_annotation_data_preview(self):
        return self.annotation_df

    def get_comparative_data_preview(self):
        return self.comparative_df
    
    def get_transcription_factor_count(self) -> int:
        """
        Counts how many proteins in the current comparative dataset match the built-in TF list.
        """
        protein_id_col = self.column_config['protein_id']
        comp_protein_col = self.column_config['comp_protein_id']
        
        if 'Gene Name' not in self.protein_df.columns:
            return 0
        
        merged_df = pd.merge(
            self.comparative_df[[comp_protein_col]],
            self.protein_df[[protein_id_col, 'Gene Name']],
            left_on=comp_protein_col,
            right_on=protein_id_col,
            how='left'
        )
        
        if 'Gene Name' not in merged_df.columns:
            return 0
            
        tf_found = merged_df['Gene Name'].isin(HUMAN_TRANSCRIPTION_FACTORS)
        return tf_found.sum()
    
    def get_comparison_groups(self) -> list:
        """
        Returns a list of unique comparison labels from the comparative data.
        """
        label_col = self.column_config['comparison_label']
        return self.comparative_df[label_col].unique().tolist()

    def filter_significant_proteins(self, comparison: str, fdr_cutoff: float, fc_cutoff: float):
        """
        Filters the comparative data for a specific comparison.
        """
        logger.info(f"Filtering for significant proteins in '{comparison}'...")
        label_col = self.column_config['comparison_label']
        fdr_col = self.column_config['fdr']
        fc_col = self.column_config['fold_change']
        
        comparison_df = self.comparative_df[self.comparative_df[label_col] == comparison].copy()
        significant_mask = (comparison_df[fdr_col] <= fdr_cutoff) & (comparison_df[fc_col].abs() >= fc_cutoff)
        significant_df = comparison_df[significant_mask]

        significant_df['Regulation'] = 'Down-regulated'
        significant_df.loc[significant_df[fc_col] > 0, 'Regulation'] = 'Up-regulated'
        
        return significant_df

    @handle_plotting_errors
    def plot_volcano(self, fdr_cutoff: float, fc_cutoff: float,
                     proteins_to_annotate=None, color_by_option='None',
                     custom_list=None, keyword=None, annotate_top_10=False,
                     arrow_config=None, **kwargs):
        """
        Generates a volcano plot with customizable arrows and themes.
        """
        logger.info("Generating volcano plot...")
        
        # Default arrow config
        if arrow_config is None:
            arrow_config = {'arrowhead': 2, 'arrowsize': 1, 'arrowwidth': 1}
        
        comp_protein_col = self.column_config['comp_protein_id']
        fc_col = self.column_config['fold_change']
        fdr_col = self.column_config['fdr']
        label_col = self.column_config['comparison_label']
        
        comparison = self.comparative_df[label_col].iloc[0]
        plot_df = self.comparative_df.copy()
        plot_df['-log10(FDR)'] = -np.log10(plot_df[fdr_col])

        protein_info_cols = ['Gene Name', 'Protein Description']
        protein_id_col = self.column_config['protein_id']
        for col in protein_info_cols:
            if col not in self.protein_df.columns: self.protein_df[col] = '' 
        plot_df = pd.merge(
            plot_df,
            self.protein_df[[protein_id_col] + protein_info_cols].drop_duplicates(subset=[protein_id_col]),
            left_on=comp_protein_col, right_on=protein_id_col, how='left'
        )

        # Color logic
        plot_df['Color'] = 'Non-significant'
        plot_df.loc[plot_df[fdr_col] <= fdr_cutoff, 'Color'] = 'Significant'
        plot_df.loc[(plot_df[fdr_col] <= fdr_cutoff) & (plot_df[fc_col] > fc_cutoff), 'Color'] = 'Up-regulated'
        plot_df.loc[(plot_df[fdr_col] <= fdr_cutoff) & (plot_df[fc_col] < -fc_cutoff), 'Color'] = 'Down-regulated'
        color_map = {'Up-regulated': '#D55E00', 'Down-regulated': '#0072B2', 'Significant': '#F0E442', 'Non-significant': '#999999'}
        
        if color_by_option == 'Custom List' and custom_list:
            mask = plot_df[comp_protein_col].isin(custom_list) | plot_df['Gene Name'].isin(custom_list)
            plot_df.loc[mask, 'Color'] = 'In Custom List'; color_map['In Custom List'] = '#009E73'
        elif color_by_option == 'Transcription Factors':
            mask = plot_df['Gene Name'].isin(HUMAN_TRANSCRIPTION_FACTORS)
            plot_df.loc[mask, 'Color'] = 'Transcription Factor'; color_map['Transcription Factor'] = '#CC79A7'
        elif color_by_option == 'Keyword Search' and keyword:
            mask = plot_df['Protein Description'].str.contains(keyword, case=False, na=False)
            plot_df.loc[mask, 'Color'] = f"Contains '{keyword}'"; color_map[f"Contains '{keyword}'"] = '#56B4E9'

        # Extract kwargs but remove global color map as Volcano uses specific significance colors
        plot_kwargs = kwargs.copy()
        plot_kwargs.pop('color_discrete_map', None)

        fig = px.scatter(
            plot_df, x=fc_col, y='-log10(FDR)', color='Color',
            color_discrete_map=color_map, hover_data=[comp_protein_col, 'Gene Name', 'Protein Description'],
            title=f"Volcano Plot for {comparison}",
            **plot_kwargs
        )
        fig.add_vline(x=fc_cutoff, line_dash="dash", line_color="grey")
        fig.add_vline(x=-fc_cutoff, line_dash="dash", line_color="grey")
        fig.add_hline(y=-np.log10(fdr_cutoff), line_dash="dash", line_color="grey")

        # Annotations
        final_annotations = set(proteins_to_annotate if proteins_to_annotate else [])
        if annotate_top_10:
            significant_df = plot_df[(plot_df[fdr_col] <= fdr_cutoff) & (plot_df[fc_col].abs() >= fc_cutoff)]
            top_10_df = significant_df.nsmallest(10, fdr_col)
            for protein in top_10_df['Gene Name'].fillna(top_10_df[comp_protein_col]).tolist():
                final_annotations.add(protein)
        
        if final_annotations:
            for protein_id in final_annotations:
                protein_data = plot_df[(plot_df[comp_protein_col] == protein_id) | (plot_df['Gene Name'] == protein_id)]
                if not protein_data.empty:
                    row = protein_data.iloc[0]
                    fig.add_annotation(
                        x=row[fc_col], y=row['-log10(FDR)'],
                        text=row.get('Gene Name') or row.get(comp_protein_col),
                        showarrow=True, 
                        arrowhead=arrow_config.get('arrowhead', 2),
                        arrowsize=arrow_config.get('arrowsize', 1),
                        arrowwidth=arrow_config.get('arrowwidth', 1),
                        ax=20, ay=-40, font=dict(color="black", size=12)
                    )
        
        fig.update_layout(height=700)
        return fig

    @handle_plotting_errors
    def plot_comparative_heatmap(self, protein_list: list, title: str = None, figsize: tuple = None, dpi: int = 150, **_kwargs):
        """
        Generates a clustered heatmap using Seaborn (Static).
        Accepts title/figsize/dpi kwargs so MplPlotManager's edit panel can re-render with user changes.
        """
        logger.info(f"Generating Seaborn heatmap for {len(protein_list)} proteins.")
        
        protein_id_col = self.column_config['protein_id']
        sample_id_col = self.column_config['sample_id']
        group_col = self.column_config['grouping']

        if not protein_list:
            raise ValueError("The list of proteins to plot cannot be empty.")
        heatmap_df = self.protein_df[self.protein_df[protein_id_col].isin(protein_list)]
        
        if 'Gene Name' in heatmap_df.columns and not heatmap_df['Gene Name'].isnull().all():
            row_labels = heatmap_df['Gene Name'].fillna(heatmap_df[protein_id_col])
        else:
            row_labels = heatmap_df[protein_id_col]
        
        intensity_matrix = heatmap_df.set_index(row_labels)[self.sample_cols]

        imputer = SimpleImputer(strategy='mean')
        imputed_matrix = imputer.fit_transform(intensity_matrix)
        scaler = StandardScaler()
        scaled_matrix = scaler.fit_transform(imputed_matrix.T).T
        scaled_df = pd.DataFrame(scaled_matrix, index=intensity_matrix.index, columns=intensity_matrix.columns)

        groups = self.annotation_df.set_index(sample_id_col)[group_col]
        unique_groups = groups.unique()
        palette = sns.color_palette("viridis", len(unique_groups))
        color_map = dict(zip(unique_groups, palette))
        col_colors = groups.map(color_map)

        show_yticklabels = len(scaled_df) <= 50 
        
        _figsize = figsize or (12, max(10, len(scaled_df) * 0.08))
        g = sns.clustermap(
            scaled_df,
            method='ward',
            cmap='RdBu_r',
            col_colors=col_colors.to_frame(),
            yticklabels=show_yticklabels,
            figsize=_figsize,
        )

        handles = [plt.Rectangle((0,0),1,1, color=color_map[group]) for group in unique_groups]
        plt.legend(handles, unique_groups, title='Group', bbox_to_anchor=(1, 1),
                   bbox_transform=plt.gcf().transFigure, frameon=False)
        g.fig.suptitle(title or "Heatmap of Protein Abundance (Z-Score)", y=1.02)

        buf = BytesIO()
        plt.savefig(buf, format="png", dpi=dpi, bbox_inches='tight')
        buf.seek(0)
        plt.close(g.fig)
        return buf

    @handle_plotting_errors
    def plot_expression_violin(self, protein_list: list, **kwargs):
        """
        Generates clean, faceted violin plots for a given list of proteins.
        Accepts kwargs for global styling.
        """
        logger.info(f"Generating clean violin plots for {len(protein_list)} proteins.")
        
        protein_id_col = self.column_config['protein_id']
        sample_id_col = self.column_config['sample_id']
        group_col = self.column_config['grouping']

        if not protein_list:
            raise ValueError("The list of proteins to plot cannot be empty.")
        
        plot_df = self.protein_df[self.protein_df[protein_id_col].isin(protein_list)]
        melted_df = plot_df.melt(
            id_vars=[protein_id_col, 'Gene Name'],
            value_vars=self.sample_cols,
            var_name=sample_id_col,
            value_name='Intensity'
        ).dropna(subset=['Intensity'])
        
        plot_data = pd.merge(melted_df, self.annotation_df, on=sample_id_col)
        plot_data = plot_data[plot_data['Intensity'] > 0].copy()
        plot_data['log2(Intensity)'] = np.log2(plot_data['Intensity'])

        valid_proteins_to_plot = plot_data[protein_id_col].unique()
        if len(valid_proteins_to_plot) == 0:
            raise ValueError("None of the selected proteins have valid intensity data to plot.")

        num_proteins = len(valid_proteins_to_plot)
        cols = 2
        rows = (num_proteins + 1) // cols
        
        subplot_titles = []
        for protein in valid_proteins_to_plot:
            gene_name = plot_data[plot_data[protein_id_col] == protein]['Gene Name'].iloc[0]
            subplot_titles.append(gene_name or protein)

        fig = make_subplots(rows=rows, cols=cols, subplot_titles=subplot_titles)

        # Colors
        default_colors = px.colors.qualitative.Plotly
        unique_groups = plot_data[group_col].unique()
        legend_added = []
        
        # Extract user custom colors if present
        custom_colors = kwargs.get('color_discrete_map', {})

        for i, protein in enumerate(valid_proteins_to_plot):
            row = (i // cols) + 1
            col = (i % cols) + 1
            protein_data = plot_data[plot_data[protein_id_col] == protein]
            
            for j, group in enumerate(unique_groups):
                group_data = protein_data[protein_data[group_col] == group]
                show_legend = group not in legend_added
                if show_legend: legend_added.append(group)
                
                # Determine color
                fill_color = custom_colors.get(group, default_colors[j % len(default_colors)])
                
                fig.add_trace(go.Violin(
                    x=group_data[group_col], y=group_data['log2(Intensity)'], name=group,
                    box_visible=True, meanline_visible=True, points='all', pointpos=0, jitter=0.05,
                    fillcolor=fill_color, line_color='black',
                    showlegend=show_legend, legendgroup=group
                ), row=row, col=col)

        # Apply template from kwargs
        layout_args = dict(
            title_text="Protein Expression Distribution", height=400 * rows,
            violingap=0, violingroupgap=0, violinmode='overlay'
        )
        
        if 'template' in kwargs:
            layout_args['template'] = kwargs['template']

        fig.update_layout(**layout_args)
        fig.update_traces(marker=dict(size=3))
        fig.update_yaxes(title_text="log2(Intensity)")

        return fig

    def run_enrichment_analysis(self, gene_list: list, organism: str = "human"):
        """
        Perform comprehensive enrichment analysis for a list of genes using Enrichr.
        """
        logger.info(f"Running Enrichr analysis for {len(gene_list)} genes on organism '{organism}'...")
        if not gene_list:
            raise ValueError("Gene list for enrichment analysis cannot be empty.")

        organism_libraries = {
            "human": ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human', 'Reactome_2022'],
            "mouse": ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2019_Mouse', 'Reactome_2022'],
        }
        source_mapping = {
            'GO_Biological_Process_2023': 'GO:BP', 'GO_Cellular_Component_2023': 'GO:CC',
            'GO_Molecular_Function_2023': 'GO:MF', 'KEGG_2021_Human': 'KEGG',
            'KEGG_2019_Mouse': 'KEGG', 'Reactome_2022': 'REAC'
        }

        if organism.lower() not in organism_libraries:
            raise ValueError(f"Organism '{organism}' is not supported by this function.")

        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr'
        genes_str = '\n'.join(gene_list)
        payload = {'list': (None, genes_str), 'description': (None, 'Pro-Visualize Analysis')}

        _retry = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
        _session = requests.Session()
        _session.mount("https://", HTTPAdapter(max_retries=_retry))
        _session.mount("http://", HTTPAdapter(max_retries=_retry))

        try:
            response = _session.post(f'{ENRICHR_URL}/addList', files=payload, timeout=30)
            response.raise_for_status()
            user_list_id = response.json()['userListId']
        except requests.exceptions.RequestException as e:
            raise RuntimeError(
                "Enrichr API may be down — try again or use the offline gene-list export. "
                f"(Network error: {e})"
            )

        all_results = []
        for library in organism_libraries[organism.lower()]:
            time.sleep(0.5)
            try:
                query_url = f'{ENRICHR_URL}/enrich?userListId={user_list_id}&backgroundType={library}'
                response = _session.get(query_url, timeout=30)
                response.raise_for_status()
                data = response.json().get(library, [])
                
                for result in data:
                    all_results.append({
                        'source': source_mapping.get(library, library),
                        'name': result[1],
                        'p_value': float(result[2]),
                        'adj_p_value': float(result[6]),
                        'genes': ";".join(result[5]),
                        'intersection_size': len(result[5])
                    })
            except requests.exceptions.RequestException as e:
                logger.warning(f"Could not retrieve results from Enrichr library {library}: {e}")
        
        if not all_results:
            return pd.DataFrame()

        results_df = pd.DataFrame(all_results).sort_values('p_value').reset_index(drop=True)
        results_df['-log10(p)'] = -np.log10(results_df['p_value'])
        return results_df

    @handle_plotting_errors
    def plot_enrichment_manhattan(self, enrichment_df: pd.DataFrame, **kwargs):
        """
        Creates an enhanced Manhattan-like plot with colored blocks for each source.
        Accepts kwargs for global theme.
        """
        if enrichment_df.empty:
            raise ValueError("Enrichment data is empty, cannot generate Manhattan plot.")

        source_order = sorted(enrichment_df['source'].unique())
        enrichment_df['source'] = pd.Categorical(enrichment_df['source'], categories=source_order, ordered=True)
        enrichment_df = enrichment_df.sort_values('source')
        enrichment_df['index'] = range(len(enrichment_df))

        fig = px.scatter(
            enrichment_df,
            x='index',
            y='-log10(p)',
            color='source',
            size='intersection_size',
            hover_name='name',
            hover_data=['p_value', 'genes'],
            title="Pathway Enrichment Manhattan Plot",
            labels={'index': 'Enrichment Terms by Source', '-log10(p)': '-log10(p-value)'},
            **kwargs
        )

        start_idx = 0
        for i, source in enumerate(source_order):
            source_data = enrichment_df[enrichment_df['source'] == source]
            if not source_data.empty:
                count = len(source_data)
                fig.add_shape(
                    type='rect',
                    x0=start_idx - 0.5, x1=start_idx + count - 0.5,
                    y0=0, y1=enrichment_df['-log10(p)'].max() * 1.05,
                    fillcolor=px.colors.qualitative.Plotly[i % len(px.colors.qualitative.Plotly)],
                    opacity=0.1, line_width=0, layer='below'
                )
                start_idx += count

        fig.update_xaxes(showticklabels=False)
        fig.update_layout(height=600)
        return fig

    @handle_plotting_errors
    def plot_enrichment_dotplot(self, enrichment_df: pd.DataFrame, terms_to_plot: list = None, **kwargs):
        """
        Generates a dot plot for enriched terms.
        Accepts kwargs for global theme.
        """
        if enrichment_df.empty:
            raise ValueError("Enrichment data is empty, cannot generate dot plot.")
        
        enrichment_df['adj_p_value'] = pd.to_numeric(enrichment_df['adj_p_value'], errors='coerce')
        enrichment_df = enrichment_df.dropna(subset=['adj_p_value'])
        enrichment_df['-log10(q)'] = -np.log10(enrichment_df['adj_p_value'] + 1e-10)
        
        if terms_to_plot:
            plot_data = enrichment_df[enrichment_df['name'].isin(terms_to_plot)]
        else:
            plot_data = enrichment_df.sort_values(by='adj_p_value').head(10)
        
        plot_data = plot_data.sort_values(by='-log10(q)', ascending=True)

        fig = px.scatter(
            plot_data,
            x='-log10(q)',
            y='name',
            size='intersection_size',
            color='intersection_size',
            color_continuous_scale='viridis',
            hover_name='name',
            hover_data=['source', 'adj_p_value', 'genes'],
            title="Top Enriched Pathways",
            **kwargs
        )
        
        min_count = plot_data['intersection_size'].min()
        max_count = plot_data['intersection_size'].max()
        tick_values = np.linspace(min_count, max_count, num=min(5, max_count - min_count + 1), dtype=int)

        fig.update_layout(
            height=max(500, len(plot_data) * 25),
            yaxis={'categoryorder':'total ascending'},
            xaxis_title="-log10(q-value)",
            yaxis_title="Pathway Term",
            coloraxis_colorbar_title_text='Gene Count',
            coloraxis_colorbar_tickvals=tick_values
        )
        return fig