import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
from io import BytesIO
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
# Import third-party libraries for plotting
from venn import venn
from upsetplot import UpSet, from_memberships
import plotly.figure_factory as ff
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.spatial import ConvexHull
import scipy.cluster.hierarchy as sch
import seaborn as sns
import matplotlib.pyplot as plt
from io import BytesIO
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from plotly.subplots import make_subplots
import plotly.graph_objects as go

logger = logging.getLogger(__name__)

class ComparativeVisualizer:
    """
    Handles data loading, validation, and visualization for comparative analysis.
    """

    def __init__(self, protein_df: pd.DataFrame, annotation_df: pd.DataFrame, 
                 comparative_df: pd.DataFrame, column_config: dict):
        """
        Initializes the visualizer with all necessary data and configuration.

        Args:
            protein_df (pd.DataFrame): DataFrame with protein abundance data.
            annotation_df (pd.DataFrame): DataFrame with sample metadata.
            comparative_df (pd.DataFrame): DataFrame with statistical comparison results.
            column_config (dict): A dictionary mapping standard roles to actual column names.
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
        
        # --- Expected columns based on the config ---
        protein_col = self.column_config['protein_id']
        sample_col = self.column_config['sample_id']
        group_col = self.column_config['grouping']
        comp_protein_col = self.column_config['comp_protein_id']
        fc_col = self.column_config['fold_change']
        fdr_col = self.column_config['fdr']
        label_col = self.column_config['comparison_label']

        # --- Validation checks ---
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
    
    def get_comparison_groups(self) -> list:
        """
        Returns a list of unique comparison labels from the comparative data.
        """
        label_col = self.column_config['comparison_label']
        return self.comparative_df[label_col].unique().tolist()

    def filter_significant_proteins(self, comparison: str, fdr_cutoff: float, fc_cutoff: float):
        """
        Filters the comparative data for a specific comparison and returns
        significantly up- and down-regulated proteins.

        Args:
            comparison (str): The comparison group to filter (e.g., 'Control/Ethanol').
            fdr_cutoff (float): The FDR threshold (e.g., 0.05).
            fc_cutoff (float): The absolute Log2 Fold Change threshold (e.g., 1.0).

        Returns:
            A pandas DataFrame containing only the significant proteins for that comparison.
        """
        logger.info(f"Filtering for significant proteins in '{comparison}' with FDR <= {fdr_cutoff} and |log2FC| >= {fc_cutoff}")
        
        # Column names from config
        label_col = self.column_config['comparison_label']
        fdr_col = self.column_config['fdr']
        fc_col = self.column_config['fold_change']
        
        # Filter for the selected comparison
        comparison_df = self.comparative_df[self.comparative_df[label_col] == comparison].copy()
        
        # Apply FDR and FC filters
        significant_mask = (comparison_df[fdr_col] <= fdr_cutoff) & (comparison_df[fc_col].abs() >= fc_cutoff)
        significant_df = comparison_df[significant_mask]

        # Add a status column for clarity
        significant_df['Regulation'] = 'Down-regulated'
        significant_df.loc[significant_df[fc_col] > 0, 'Regulation'] = 'Up-regulated'
        
        return significant_df
    

    # ... (inside the ComparativeVisualizer class)

    def plot_volcano(self, fdr_cutoff: float, fc_cutoff: float,
                     proteins_to_annotate=None, color_by_option='None',
                     custom_list=None, keyword=None):
        """
        Generates a volcano plot with extensive customization options.
        """
        logger.info(f"Generating volcano plot with options: {color_by_option}")
        
        # 1. Get column names from config
        comp_protein_col = self.column_config['comp_protein_id']
        fc_col = self.column_config['fold_change']
        fdr_col = self.column_config['fdr']
        label_col = self.column_config['comparison_label']
        
        # 2. Prepare the data for the selected comparison
        comparison = self.comparative_df[label_col].iloc[0] # Assume df is pre-filtered
        plot_df = self.comparative_df.copy()
        plot_df['-log10(FDR)'] = -np.log10(plot_df[fdr_col])

        # Merge with protein info for tooltips and searching
        protein_info_cols = ['Gene Name', 'Protein Description']
        protein_id_col = self.column_config['protein_id']
        
        # Ensure columns exist in the main protein_df before merging
        for col in protein_info_cols:
            if col not in self.protein_df.columns:
                self.protein_df[col] = '' 
                
        plot_df = pd.merge(
            plot_df,
            self.protein_df[[protein_id_col] + protein_info_cols].drop_duplicates(subset=[protein_id_col]),
            left_on=comp_protein_col, right_on=protein_id_col, how='left'
        )

        # 3. Define significance categories
        plot_df['Color'] = 'Non-significant'
        plot_df.loc[plot_df[fdr_col] <= fdr_cutoff, 'Color'] = 'Significant'
        plot_df.loc[(plot_df[fdr_col] <= fdr_cutoff) & (plot_df[fc_col] > fc_cutoff), 'Color'] = 'Up-regulated'
        plot_df.loc[(plot_df[fdr_col] <= fdr_cutoff) & (plot_df[fc_col] < -fc_cutoff), 'Color'] = 'Down-regulated'

        # 4. Apply custom coloring logic (overrides significance)
        color_map = {
            'Up-regulated': '#D55E00',
            'Down-regulated': '#0072B2',
            'Significant': '#F0E442',
            'Non-significant': '#999999'
        }
        
        if color_by_option == 'Custom List' and custom_list:
            mask = plot_df[comp_protein_col].isin(custom_list) | plot_df['Gene Name'].isin(custom_list)
            plot_df.loc[mask, 'Color'] = 'In Custom List'
            color_map['In Custom List'] = '#009E73' # Teal
        elif color_by_option == 'Transcription Factors':
            mask = plot_df['Gene Name'].isin(self.HUMAN_TRANSCRIPTION_FACTORS)
            plot_df.loc[mask, 'Color'] = 'Transcription Factor'
            color_map['Transcription Factor'] = '#CC79A7' # Pink/Purple
        elif color_by_option == 'Keyword Search' and keyword:
            mask = plot_df['Protein Description'].str.contains(keyword, case=False, na=False)
            plot_df.loc[mask, 'Color'] = f"Contains '{keyword}'"
            color_map[f"Contains '{keyword}'"] = '#56B4E9' # Sky Blue

        # 5. Create the scatter plot
        fig = px.scatter(
            plot_df,
            x=fc_col,
            y='-log10(FDR)',
            color='Color',
            color_discrete_map=color_map,
            hover_data=[comp_protein_col, 'Gene Name', 'Protein Description'],
            title=f"Volcano Plot for {comparison}"
        )

        # 6. Add threshold lines
        fig.add_vline(x=fc_cutoff, line_dash="dash", line_color="grey")
        fig.add_vline(x=-fc_cutoff, line_dash="dash", line_color="grey")
        fig.add_hline(y=-np.log10(fdr_cutoff), line_dash="dash", line_color="grey")

        # 7. Add annotations
        if proteins_to_annotate:
            for protein_id in proteins_to_annotate:
                protein_data = plot_df[(plot_df[comp_protein_col] == protein_id) | (plot_df['Gene Name'] == protein_id)]
                if not protein_data.empty:
                    row = protein_data.iloc[0]
                    fig.add_annotation(
                        x=row[fc_col], y=row['-log10(FDR)'],
                        text=row.get('Gene Name') or row.get(comp_protein_col),
                        showarrow=True, arrowhead=2, arrowsize=1,
                        ax=20, ay=-40, font=dict(color="black", size=12)
                    )
        
        fig.update_layout(height=700)
        return fig
    

    def plot_comparative_heatmap(self, protein_list: list):
        """
        Generates a clustered heatmap for a given list of proteins using Seaborn,
        with a color bar annotating the sample groups.

        Args:
            protein_list (list): A list of Protein IDs to include in the heatmap.
        """
        logger.info(f"Generating Seaborn heatmap for {len(protein_list)} proteins.")
        
        # 1. Get column names from config
        protein_id_col = self.column_config['protein_id']
        sample_id_col = self.column_config['sample_id']
        group_col = self.column_config['grouping']

        # 2. Filter protein data for the selected list
        if not protein_list:
            raise ValueError("The list of proteins to plot cannot be empty.")
        heatmap_df = self.protein_df[self.protein_df[protein_id_col].isin(protein_list)]
        
        if 'Gene Name' in heatmap_df.columns and not heatmap_df['Gene Name'].isnull().all():
            row_labels = heatmap_df['Gene Name'].fillna(heatmap_df[protein_id_col])
        else:
            row_labels = heatmap_df[protein_id_col]
        
        intensity_matrix = heatmap_df.set_index(row_labels)[self.sample_cols]

        # 3. Impute and Z-score scale the data
        imputer = SimpleImputer(strategy='mean')
        imputed_matrix = imputer.fit_transform(intensity_matrix)
        scaler = StandardScaler()
        scaled_matrix = scaler.fit_transform(imputed_matrix.T).T
        scaled_df = pd.DataFrame(scaled_matrix, index=intensity_matrix.index, columns=intensity_matrix.columns)

        # 4. Create the color mapping for the sample groups
        groups = self.annotation_df.set_index(sample_id_col)[group_col]
        unique_groups = groups.unique()
        palette = sns.color_palette("viridis", len(unique_groups))
        color_map = dict(zip(unique_groups, palette))
        col_colors = groups.map(color_map)

        # 5. Generate the clustermap
        # Dynamically hide row labels if there are too many proteins to avoid clutter
        show_yticklabels = len(scaled_df) <= 50 
        
        g = sns.clustermap(
            scaled_df,
            method='ward',       # Hierarchical clustering method
            cmap='RdBu_r',       # Diverging color map for Z-scores
            col_colors=col_colors.to_frame(), # Add the group color bar
            yticklabels=show_yticklabels,
            figsize=(12, max(10, len(scaled_df) * 0.08)) # Dynamic height
        )
        
        # Add a proper legend for the group colors
        handles = [plt.Rectangle((0,0),1,1, color=color_map[group]) for group in unique_groups]
        plt.legend(handles, unique_groups, title='Group', bbox_to_anchor=(1, 1), 
                   bbox_transform=plt.gcf().transFigure, frameon=False)
        g.fig.suptitle("Heatmap of Protein Abundance (Z-Score)", y=1.02)
        
        # 6. Save plot to an in-memory buffer to pass to Streamlit
        buf = BytesIO()
        plt.savefig(buf, format="png", bbox_inches='tight')
        buf.seek(0)
        plt.close(g.fig)
        return buf
    

    def plot_expression_violin(self, protein_list: list):
        """
        Generates clean, faceted violin plots for a given list of proteins using
        the logic from the original ProteinVisualization class for a cleaner look.
        """
        logger.info(f"Generating clean violin plots for {len(protein_list)} proteins.")
        
        # 1. Get column names and prepare data (same as before)
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

        # 2. Set up subplots
        num_proteins = len(protein_list)
        cols = 2  # Arrange plots in 2 columns
        rows = (num_proteins + 1) // cols
        
        # Use Gene Name for subplot titles if available
        subplot_titles = []
        for protein in protein_list:
            gene_name = plot_data[plot_data[protein_id_col] == protein]['Gene Name'].iloc[0]
            subplot_titles.append(gene_name or protein)

        fig = make_subplots(rows=rows, cols=cols, subplot_titles=subplot_titles)

        # 3. Create the plot using a loop, similar to your original code
        colors = px.colors.qualitative.Plotly
        unique_groups = plot_data[group_col].unique()
        legend_added = []

        for i, protein in enumerate(protein_list):
            row = (i // cols) + 1
            col = (i % cols) + 1
            
            protein_data = plot_data[plot_data[protein_id_col] == protein]
            
            for j, group in enumerate(unique_groups):
                group_data = protein_data[protein_data[group_col] == group]
                show_legend = group not in legend_added
                if show_legend:
                    legend_added.append(group)
                
                fig.add_trace(go.Violin(
                    x=group_data[group_col],
                    y=group_data['log2(Intensity)'],
                    name=group,
                    box_visible=True,
                    meanline_visible=True,
                    points='all',
                    pointpos=0,
                    jitter=0.05,
                    fillcolor=colors[j % len(colors)],
                    line_color='black',
                    showlegend=show_legend,
                    legendgroup=group
                ), row=row, col=col)

        # 4. Final layout adjustments
        fig.update_layout(
            title_text="Protein Expression Distribution",
            height=400 * rows, # Dynamic height
            violingap=0, 
            violingroupgap=0, 
            violinmode='overlay'
        )
        fig.update_traces(marker=dict(size=3)) # Smaller points
        fig.update_yaxes(title_text="log2(Intensity)")

        return fig