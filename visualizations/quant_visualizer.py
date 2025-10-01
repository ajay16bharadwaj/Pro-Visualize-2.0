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




# Set up a logger for this module
logger = logging.getLogger(__name__)

class QuantificationVisualizer:
    """
    Handles data processing and visualization for quantification analysis.
    """

    def __init__(self, protein_df: pd.DataFrame, annotation_df: pd.DataFrame = None,
                 protein_col: str = 'Protein', sample_col: str = 'Level3',
                 group_col: str = 'attribute_ExperimentalGroup'):
        """
        Initializes the visualizer with data and configurable column names.

        Args:
            protein_df (pd.DataFrame): DataFrame with protein abundance data.
            annotation_df (pd.DataFrame, optional): DataFrame with sample metadata.
            protein_col (str): The name of the column containing unique protein identifiers.
            sample_col (str): The name of the column in the annotation file used for matching sample names.
            group_col (str): The name of the column in the annotation file defining experimental groups.
        """
        logger.info("Initializing QuantificationVisualizer with custom columns...")
        
        # Store column names
        self.protein_col = protein_col
        self.sample_col = sample_col
        self.group_col = group_col

        self._validate_input_data(protein_df, annotation_df)
        self.protein_df = protein_df.copy()
        self.annotation_df = annotation_df.copy() if annotation_df is not None else None
        
        if self.annotation_df is not None:
            self.sample_cols = self.annotation_df[self.sample_col].unique().tolist()
            self.experimental_groups = self.annotation_df[self.group_col].unique().tolist()
        else:
            known_non_sample_cols = [
                'ProteinIds', 'Gene Name', 'Protein Description', 'Tissue Specificity',
                'Gene ontology (biological process)', 'Gene ontology (cellular component)',
                'Gene ontology (molecular function)', 'Subcellular Location[CC]', self.protein_col
            ]
            self.sample_cols = [col for col in self.protein_df.columns if col not in known_non_sample_cols]
            self.experimental_groups = []
        
        logger.info(f"Identified {len(self.sample_cols)} samples for analysis.")

    def _validate_input_data(self, protein_df, annotation_df):
        """Performs initial validation using the specified column names."""
        if self.protein_col not in protein_df.columns:
            raise ValueError(f"Protein data must contain the specified protein column: '{self.protein_col}'.")
        
        if annotation_df is not None:
            required_meta_cols = [self.sample_col, self.group_col]
            missing_cols = [col for col in required_meta_cols if col not in annotation_df.columns]
            if missing_cols:
                raise ValueError(f"Annotation data is missing required columns: {', '.join(missing_cols)}")

    def _get_protein_sets_by_group(self) -> dict:
        """Groups proteins into sets based on their experimental group."""
        if self.annotation_df is None:
            raise ValueError("Annotation data is required to group proteins.")
            
        protein_sets = {}
        grouped_samples = self.annotation_df.groupby(self.group_col)[self.sample_col].apply(list)
        
        for group, samples in grouped_samples.items():
            group_df = self.protein_df[[self.protein_col] + [s for s in samples if s in self.protein_df.columns]]
            detected_proteins = group_df.dropna(subset=samples, how='all')[self.protein_col].unique()
            protein_sets[group] = set(detected_proteins)
            
        return protein_sets

    def plot_venn_diagram(self, selected_groups: list):
        """Generates a Venn diagram for the intersection of proteins."""
        if not (2 <= len(selected_groups) <= 6):
            raise ValueError("Venn diagrams are supported for 2 to 6 groups.")
            
        all_protein_sets = self._get_protein_sets_by_group()
        sets_to_plot = {group: all_protein_sets.get(group, set()) for group in selected_groups}

        fig, ax = plt.subplots(figsize=(4, 4))
        venn(sets_to_plot, ax=ax)
        plt.tight_layout()
        
        return fig

    def plot_upset(self):
        """
        Generates an UpSet plot to visualize intersections of proteins across all experimental groups.
        This implementation is optimized for performance.
        """
        protein_sets = self._get_protein_sets_by_group()

        # Get the universe of all proteins present
        all_proteins = set.union(*protein_sets.values())
        
        # Create the membership list required by from_memberships
        memberships = []
        for protein in all_proteins:
            protein_membership = [
                group for group, proteins in protein_sets.items() if protein in proteins
            ]
            memberships.append(protein_membership)

        # Convert the membership list to the format required by the UpSet library
        upset_data = from_memberships(memberships)

        # Create the plot object, showing all intersections by default
        upset_plot = UpSet(upset_data, subset_size='count', show_counts='%d', sort_by='cardinality')
        
        # Render the plot to a matplotlib Figure
        fig = plt.figure(figsize=(4, 8))
        upset_plot.plot(fig=fig)

        fig.suptitle("Protein Intersections Across Experimental Groups", fontsize=12)
        return fig
    
    def plot_protein_counts(self):
        """
        Generates a bar plot of the number of quantified proteins per sample.

        The plot is color-coded by experimental group if annotation data is available.
        Otherwise, a single color is used for all bars.

        Returns:
            A Plotly figure object.
        """
        logger.info("Generating protein counts per sample plot...")
        
        # Calculate the number of non-NA values for each sample column
        protein_counts = self.protein_df[self.sample_cols].notna().sum().reset_index()
        protein_counts.columns = ['Sample', 'ProteinCount']

        plot_args = {
            "x": 'Sample',
            "y": 'ProteinCount',
            "title": "Quantified Proteins per Sample",
            "labels": {'Sample': 'Sample', 'ProteinCount': 'Number of Proteins'}
        }

        # Check if annotation data is available for coloring
        if self.annotation_df is not None:
            # Merge with annotation to get group information
            plot_data = pd.merge(protein_counts, self.annotation_df, 
                                 left_on='Sample', right_on=self.sample_col, how='left')
            
            # Sort bars for better visualization, grouped by experimental group
            plot_data = plot_data.sort_values(by=[self.group_col, 'Sample'])
            
            plot_args["color"] = self.group_col
            plot_args["category_orders"] = {'Sample': plot_data['Sample'].tolist()} # Enforce sorting
        else:
            plot_data = protein_counts.sort_values(by='ProteinCount', ascending=False)
            plot_args["category_orders"] = {'Sample': plot_data['Sample'].tolist()}


        fig = px.bar(plot_data, **plot_args)
        
        fig.update_layout(
            xaxis_tickangle=-45,
            height=500,
            width=800,
            margin=dict(b=120) # Add bottom margin for sample names
        )
        return fig
    

    def plot_missing_values_heatmap(self):
        """
        Generates an enhanced heatmap to visualize the pattern of missing proteins.
        This version ensures unique samples are plotted to prevent duplicate labels.
        """
        logger.info("Generating enhanced missing values heatmap...")
        if self.annotation_df is None:
            raise ValueError("Annotation data is required to group samples for the missing values heatmap.")

        # --- FIX: A more robust method to get a unique, sorted list of samples ---
        # 1. Start with the unique sample identifiers from the annotation file.
        meta_unique_samples = self.annotation_df[[self.sample_col, self.group_col]].drop_duplicates(subset=[self.sample_col])

        # 2. Find the intersection of samples that exist in BOTH the protein data and the unique metadata.
        valid_samples = list(set(self.sample_cols) & set(meta_unique_samples[self.sample_col]))
        
        # 3. Filter the unique metadata to this final valid list and sort it by group.
        sorted_meta = meta_unique_samples[meta_unique_samples[self.sample_col].isin(valid_samples)].sort_values(by=self.group_col)
        sorted_samples = sorted_meta[self.sample_col].tolist()

        # 4. Filter for proteins with at least one missing value across the valid samples
        protein_subset = self.protein_df[sorted_samples]
        missing_mask = protein_subset.isnull().any(axis=1)
        df_missing = self.protein_df.loc[missing_mask]

        if df_missing.empty:
            fig = go.Figure()
            fig.update_layout(title="No Missing Proteins Found", xaxis_showticklabels=False, yaxis_showticklabels=False)
            return fig

        binary_matrix = df_missing[sorted_samples].notna().astype(int)

        # 5. Create plots (logic is unchanged)
        fig = make_subplots(rows=2, cols=1, row_heights=[0.9, 0.1], vertical_spacing=0.02)

        fig.add_trace(go.Heatmap(
            z=binary_matrix.values, x=binary_matrix.columns, y=[f"Protein {i}" for i in range(len(binary_matrix))],
            colorscale=[[0, 'white'], [1, 'black']], showscale=False, name=""
        ), row=1, col=1)

        group_names = sorted_meta[self.group_col].unique()
        group_map = {name: i for i, name in enumerate(group_names)}
        group_colors_numeric = sorted_meta[self.group_col].map(group_map)
        color_scale = px.colors.qualitative.Plotly
        
        fig.add_trace(go.Heatmap(
            z=[group_colors_numeric.values], x=sorted_samples, y=['Group'],
            colorscale=[[i / (len(group_names) - 1), color_scale[i % len(color_scale)]] for i in range(len(group_names))] if len(group_names) > 1 else [[0, color_scale[0]], [1, color_scale[0]]],
            showscale=False, name=""
        ), row=2, col=1)

        # 6. Add legend (logic is unchanged)
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name='Observed', marker=dict(size=10, color='black', symbol='square')))
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name='Missing', marker=dict(size=10, color='white', symbol='square', line=dict(width=1, color='black'))))
        for i, group in enumerate(group_names):
            fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=group, marker=dict(size=10, color=color_scale[i % len(color_scale)], symbol='square')))

        # 7. Final layout (logic is unchanged)
        fig.update_layout(
            title_text='Missing Proteins Pattern', height=600, xaxis_tickangle=-90,
            xaxis_tickfont=dict(size=8), yaxis_title="Proteins missing<br>from at least one sample",
            yaxis_showticklabels=False, yaxis_ticks="", yaxis2_showticklabels=False, yaxis2_ticks="",
            margin=dict(b=150, l=80), legend_title_text='Legend', legend_tracegroupgap=20
        )
        return fig
    
    def plot_missing_value_distribution(self):
        """
        Generates density and cumulative fraction plots for mean protein intensities,
        comparing proteins with and without missing values.
        """
        logger.info("Generating missing value distribution plots...")

        # 1. Prepare data (unchanged)
        df_long = self.protein_df.melt(
            id_vars=[self.protein_col],
            value_vars=self.sample_cols,
            var_name='Sample',
            value_name='Intensity'
        )
        df_long['Intensity'] = df_long['Intensity'].replace(0, np.nan)
        protein_stats = df_long.groupby(self.protein_col)['Intensity'].agg(
            mean_intensity=lambda x: x.mean(skipna=True),
            has_missing_values=lambda x: x.isnull().any()
        ).dropna(subset=['mean_intensity'])
        protein_stats['log2_intensity'] = np.log2(protein_stats['mean_intensity'])
        
        dist_data_missing = protein_stats[protein_stats['has_missing_values'] == True]['log2_intensity']
        dist_data_observed = protein_stats[protein_stats['has_missing_values'] == False]['log2_intensity']

        # --- NEW: Check for empty data groups before plotting ---
        if dist_data_missing.empty or dist_data_observed.empty:
            message = "No missing values found." if dist_data_missing.empty else "All proteins have missing values."
            fig = go.Figure()
            fig.update_layout(
                title="Cannot Generate Plot",
                xaxis={'visible': False},
                yaxis={'visible': False},
                annotations=[{"text": f"{message}<br>Comparison requires both groups.", "showarrow": False}]
            )
            return fig

        # 2. Prepare data for plots (unchanged)
        def get_ecdf_data(series):
            sorted_series = series.sort_values()
            cumulative_fraction = np.arange(1, len(sorted_series) + 1) / len(sorted_series)
            return pd.DataFrame({'log2_intensity': sorted_series, 'cumulative_fraction': cumulative_fraction})
        ecdf_missing = get_ecdf_data(dist_data_missing)
        ecdf_observed = get_ecdf_data(dist_data_observed)

        # 3. Create plots (unchanged)
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=("Density of missing values", "Cumulative fraction of missing values")
        )
        dist_fig = ff.create_distplot(
            [dist_data_observed, dist_data_missing],
            ['FALSE', 'TRUE'],
            show_hist=False, show_rug=False,
            colors=['#e15759', '#4e79a7']
        )
        fig.add_trace(dist_fig['data'][0], row=1, col=1)
        fig.add_trace(dist_fig['data'][1], row=1, col=1)
        fig.add_trace(go.Scatter(
            x=ecdf_observed['log2_intensity'], y=ecdf_observed['cumulative_fraction'],
            mode='lines', name='FALSE', line=dict(color='#e15759')
        ), row=1, col=2)
        fig.add_trace(go.Scatter(
            x=ecdf_missing['log2_intensity'], y=ecdf_missing['cumulative_fraction'],
            mode='lines', name='TRUE', line=dict(color='#4e79a7')
        ), row=1, col=2)

        # 4. Final layout (unchanged)
        fig.update_layout(
            height=500, showlegend=True, legend_title_text='Missing values',
            xaxis_title_text='log2 Intensity', yaxis_title_text='Density',
            xaxis2_title_text='log2 Intensity', yaxis2_title_text='Cumulative fraction'
        )
        fig.data[0].showlegend = False
        fig.data[1].showlegend = False

        return fig
    
    def plot_protein_overlap(self):
        """
        Generates a bar plot showing the number of proteins identified in a given
        number of samples.
        """
        logger.info("Generating protein overlap plot...")
        
        # Count in how many samples each protein appears (is not null)
        protein_presence = self.protein_df[self.sample_cols].notna().sum(axis=1)
        
        # Count how many proteins fall into each presence category
        overlap_counts = protein_presence.value_counts().sort_index()
        
        plot_data = overlap_counts.reset_index()
        plot_data.columns = ['Num_Samples', 'Protein_Count']

        fig = px.bar(
            plot_data,
            x='Num_Samples',
            y='Protein_Count',
            title="Protein Overlap Across All Samples",
            labels={
                'Num_Samples': 'Identified in X Samples',
                'Protein_Count': 'Number of Shared Proteins'
            }
        )
        fig.update_layout(
            xaxis_type='category',
            height=500,
            width=600,
            margin=dict(b=100)
        )
        return fig

    def plot_protein_coverage_chart(self):
        """
        Generates a bar chart that visualizes the total number of unique proteins,
        with each protein represented as a segment in a single bar.
        """
        logger.info("Generating protein coverage plot...")

        # Create a dataframe where each protein is a row, for stacking in the bar chart
        plot_data = pd.DataFrame({
            'protein_id': self.protein_df[self.protein_col].unique(),
            'value': 1,
            'category': 'all'
        })

        total_proteins = len(plot_data)

        # Create a stacked bar chart. Each protein is a tiny segment of the bar.
        fig = px.bar(
            plot_data,
            x='category',
            y='value',
            title="Protein Coverage",
            labels={'category': ''} # Hide x-axis label
        )

        # Style the plot to match the monochrome, segmented example
        fig.update_layout(
            showlegend=False,
            plot_bgcolor='white',
            yaxis_title="Number of proteins",
            xaxis_showticklabels=True,
            uniformtext_minsize=8,
            uniformtext_mode='hide',
            height=500
        )
        
        # This loop sets all segments to the same dark color and adds a thin white line
        # between them to create the segmented effect.
        fig.for_each_trace(lambda t: t.update(
            marker_color='#333333',
            marker_line_color='white',
            marker_line_width=0.5
        ))
        
        # Ensure the y-axis is scaled correctly and has a grid
        fig.update_yaxes(
            range=[0, total_proteins * 1.05], 
            gridcolor='#e5e5e5'
        )
        fig.update_xaxes(linecolor='black', showline=True)

        return fig
    
    def plot_intensity_distribution(self):
        """
        Generates a faceted density plot of protein intensities, grouped by
        experimental condition and colored by sample.
        """
        logger.info("Generating intensity distribution plot...")
        if self.annotation_df is None:
            raise ValueError("Annotation data is required to generate the intensity distribution plot.")

        # 1. Melt the protein data to a long format
        melted_data = self.protein_df.melt(
            id_vars=[self.protein_col],
            value_vars=self.sample_cols,
            var_name='Sample',
            value_name='Intensity'
        ).dropna(subset=['Intensity'])

        # 2. Merge with annotation data to get group information
        merged_data = melted_data.merge(self.annotation_df, left_on='Sample', right_on=self.sample_col)
        
        # 3. Calculate log10 intensity, handling zero or negative values
        merged_data = merged_data[merged_data['Intensity'] > 0].copy()
        merged_data['log10(Intensity)'] = np.log10(merged_data['Intensity'])

        # 4. Create the FacetGrid plot using seaborn
        g = sns.FacetGrid(
            merged_data,
            col=self.group_col,
            hue='Sample',
            col_wrap=3, # Adjust number of columns in the plot grid
            height=4,
            aspect=1.2,
            sharex=True,
            sharey=True
        )
        g.map(sns.kdeplot, "log10(Intensity)", fill=True, alpha=0.2)
        g.map(sns.kdeplot, "log10(Intensity)")

        # 5. Final styling
        g.set_axis_labels("log10(Intensity)", "Density")
        g.set_titles(col_template="{col_name}")
        g.fig.suptitle("Protein Intensity Density Distribution", y=1.03, fontsize=16)
        
        # 6. Save the plot to an in-memory buffer
        buf = BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight")
        buf.seek(0)
        plt.close(g.fig) # Close the figure to free up memory

        return buf


    HUMAN_TRANSCRIPTION_FACTORS = set([
    'ATF1', 'ATF2', 'ATF3', 'ATF4', 'ATF5', 'ATF6', 'ATF7', 'CREB1', 'CREB3',
    'CREB5', 'FOS', 'FOSB', 'JUN', 'JUNB', 'JUND', 'MYC', 'MYCN', 'MAX', 'MAD',
    'CEBPA', 'CEBPB', 'CEBPD', 'EGR1', 'EGR2', 'EGR3', 'EGR4', 'SP1', 'SP2',
    'SP3', 'SP4', 'KLF1', 'KLF2', 'KLF4', 'KLF5', 'STAT1', 'STAT2', 'STAT3',
    'STAT4', 'STAT5A', 'STAT5B', 'STAT6', 'NFKB1', 'NFKB2', 'RELA', 'RELB',
    'TP53', 'TP63', 'TP73', 'SOX2', 'SOX9', 'POU5F1', 'NANOG', 'GATA1', 'GATA2',
    'GATA3', 'GATA4', 'GATA6', 'RUNX1', 'RUNX2', 'RUNX3', 'TCF3', 'TCF4', 'LEF1'
]) 

    def get_transcription_factor_count(self) -> int:
        """Counts how many proteins in the data match the built-in TF list."""
        if 'Gene Name' not in self.protein_df.columns:
            return 0
        
        tf_found = self.protein_df['Gene Name'].isin(self.HUMAN_TRANSCRIPTION_FACTORS)
        return tf_found.sum()   


    def plot_protein_rank_order(self, proteins_to_annotate=None, color_by_option='None',
                                custom_list=None, keyword=None, annotate_highlighted=False):
        """
        Generates a protein rank order plot with extensive customization options.

        Args:
            proteins_to_annotate (list, optional): List of proteins to explicitly mark with arrows.
            color_by_option (str): The coloring strategy.
            custom_list (list, optional): A user-provided list of proteins for coloring.
            keyword (str, optional): A keyword to search for in protein descriptions.
            annotate_highlighted (bool): If True, automatically add arrows for colored proteins.
        """
        logger.info(f"Generating protein rank order plot with options: {color_by_option}")

        # 1. Calculate average intensity and rank (code is unchanged)
        avg_intensity = self.protein_df[self.sample_cols].mean(axis=1)
        avg_intensity = avg_intensity[avg_intensity > 0]
        rank_df = pd.DataFrame({
            self.protein_col: self.protein_df.loc[avg_intensity.index, self.protein_col],
            'log10_intensity': np.log10(avg_intensity)
        }).sort_values('log10_intensity', ascending=False).reset_index(drop=True)
        rank_df['Rank'] = rank_df.index + 1
        protein_info_cols = ['Gene Name', 'Protein Description']
        for col in protein_info_cols:
            if col not in self.protein_df.columns: self.protein_df[col] = ''
        rank_df = pd.merge(
            rank_df,
            self.protein_df[[self.protein_col] + protein_info_cols].drop_duplicates(),
            on=self.protein_col, how='left'
        )

        # 2. Apply coloring logic (code is unchanged)
        rank_df['Color'] = 'Other Proteins'
        color_map = {'Other Proteins': '#4e79a7'}
        highlighted_proteins = pd.Series(dtype=str)

        if color_by_option == 'Custom List' and custom_list:
            mask = rank_df[self.protein_col].isin(custom_list) | rank_df['Gene Name'].isin(custom_list)
            rank_df.loc[mask, 'Color'] = 'In Custom List'
            color_map['In Custom List'] = '#f28e2b'
            highlighted_proteins = rank_df.loc[mask, 'Gene Name'].fillna(rank_df.loc[mask, self.protein_col])
        elif color_by_option == 'Transcription Factors':
            mask = rank_df['Gene Name'].isin(self.HUMAN_TRANSCRIPTION_FACTORS)
            rank_df.loc[mask, 'Color'] = 'Transcription Factor'
            color_map['Transcription Factor'] = '#e15759'
            highlighted_proteins = rank_df.loc[mask, 'Gene Name']
        elif color_by_option == 'Keyword Search' and keyword:
            mask = rank_df['Protein Description'].str.contains(keyword, case=False, na=False)
            rank_df.loc[mask, 'Color'] = f"Contains '{keyword}'"
            color_map[f"Contains '{keyword}'"] = '#59a14f'
            highlighted_proteins = rank_df.loc[mask, 'Gene Name'].fillna(rank_df.loc[mask, self.protein_col])
            
        # 3. Create the scatter plot (code is unchanged)
        fig = px.scatter(
            rank_df, x='Rank', y='log10_intensity', color='Color',
            color_discrete_map=color_map,
            hover_data=[self.protein_col, 'Gene Name', 'Protein Description'],
            title="Protein Rank Order by Average Intensity"
        )

        # 4. Add annotations
        # Combine explicitly selected and automatically highlighted proteins for annotation
        final_annotations = set(proteins_to_annotate if proteins_to_annotate else [])
        if annotate_highlighted and not highlighted_proteins.empty:
            # Add top 10 highlighted proteins to the annotation set to prevent clutter
            for protein in highlighted_proteins.head(10).tolist():
                final_annotations.add(protein)

        if final_annotations:
            for protein_id in final_annotations:
                protein_data = rank_df[(rank_df[self.protein_col] == protein_id) | (rank_df['Gene Name'] == protein_id)]
                if not protein_data.empty:
                    row = protein_data.iloc[0]
                    fig.add_annotation(
                        x=row['Rank'], y=row['log10_intensity'],
                        text=row.get('Gene Name') or row.get(self.protein_col),
                        showarrow=True, arrowhead=2, arrowsize=1,
                        ax=20, ay=-40, font=dict(color="black", size=12)
                    )

        fig.update_layout(
            xaxis_title="Protein Rank", yaxis_title="log10(Average Intensity)",
            legend_title_text="Category", height=600
        )
        return fig

    def _prepare_data_for_clustering(self):
        """
        Prepares the protein data for PCA and other clustering methods by imputing
        missing values and scaling the data. This version correctly handles
        samples that may be entirely empty.
        """
        logger.info("Preparing data for clustering...")
        
        numeric_df = self.protein_df[self.sample_cols].copy()
        numeric_df.replace(0, np.nan, inplace=True)

        # --- FIX: Drop columns that are entirely NaN before imputation ---
        # This prevents a shape mismatch if some samples have no data.
        numeric_df.dropna(axis=1, how='all', inplace=True)
        
        # Get the list of remaining valid samples
        valid_samples = numeric_df.columns.tolist()

        # Impute missing values (e.g., with the mean of each protein)
        imputer = SimpleImputer(strategy='mean')
        imputed_data = imputer.fit_transform(numeric_df)
        df_imputed = pd.DataFrame(imputed_data, columns=valid_samples, index=numeric_df.index)

        # Scale the data (important for PCA)
        scaler = StandardScaler()
        # We transpose (.T) because PCA works on samples (rows) vs features (columns)
        scaled_data = scaler.fit_transform(df_imputed.T)
        
        # Return scaled data and the CORRECT list of sample names
        return scaled_data, df_imputed.columns

    def plot_pca_by_annotation(self, color_by: str, symbol_by: str = None, 
                               pc_x: int = 1, pc_y: int = 2, show_labels: bool = False):
        """
        Generates a PCA plot where samples are colored and shaped by metadata columns.
        """
        if self.annotation_df is None:
            raise ValueError("Annotation data is required for this PCA plot.")
            
        scaled_data, samples = self._prepare_data_for_clustering()
        
        pca = PCA(n_components=max(pc_x, pc_y))
        principal_components = pca.fit_transform(scaled_data)
        
        pc_x_label = f'PC{pc_x} ({pca.explained_variance_ratio_[pc_x-1]*100:.1f}%)'
        pc_y_label = f'PC{pc_y} ({pca.explained_variance_ratio_[pc_y-1]*100:.1f}%)'
        
        pca_df = pd.DataFrame(
            principal_components[:, [pc_x-1, pc_y-1]],
            columns=[pc_x_label, pc_y_label],
            index=samples
        ).reset_index().rename(columns={'index': self.sample_col})

        plot_df = pd.merge(pca_df, self.annotation_df, on=self.sample_col)

        fig = px.scatter(
            plot_df,
            x=pc_x_label,
            y=pc_y_label,
            color=color_by,
            symbol=symbol_by,
            title=f"PCA Colored by '{color_by}'",
            hover_data=[self.sample_col],
            # --- NEW: Conditionally add text labels ---
            text=plot_df[self.sample_col] if show_labels else None
        )
        fig.update_traces(
            marker=dict(size=10),
            textposition='top center'
        )
        fig.update_layout(height=600)
        return fig

    def plot_pca_with_clusters(self, n_clusters: int = 3, show_labels: bool = False):
        """
        Generates a PCA plot with unsupervised k-means clustering and convex hulls.
        """
        scaled_data, samples = self._prepare_data_for_clustering()

        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(scaled_data)

        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        clusters = kmeans.fit_predict(principal_components)

        pca_df = pd.DataFrame(principal_components, columns=['PC1', 'PC2'])
        pca_df['Cluster'] = clusters
        pca_df['Sample'] = samples

        fig = go.Figure()
        colors = px.colors.qualitative.Plotly

        for cluster_id in range(n_clusters):
            cluster_points = pca_df[pca_df['Cluster'] == cluster_id]
            color = colors[cluster_id % len(colors)]
            
            fig.add_trace(go.Scatter(
                x=cluster_points['PC1'], y=cluster_points['PC2'],
                # --- NEW: Conditionally change mode to include text ---
                mode='markers+text' if show_labels else 'markers',
                marker=dict(color=color, size=10),
                name=f'Cluster {cluster_id}',
                text=cluster_points['Sample'],
                hoverinfo='text',
                textposition='top center'
            ))

            if len(cluster_points) > 2:
                hull = ConvexHull(cluster_points[['PC1', 'PC2']].values)
                hull_points = np.append(hull.vertices, hull.vertices[0])
                fig.add_trace(go.Scatter(
                    x=cluster_points['PC1'].iloc[hull_points],
                    y=cluster_points['PC2'].iloc[hull_points],
                    fill="toself",
                    fillcolor=color,
                    opacity=0.2,
                    line=dict(color=color, width=1),
                    showlegend=False
                ))
        
        fig.update_layout(
            title=f'PCA with {n_clusters} Clusters',
            xaxis_title=f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)',
            yaxis_title=f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)',
            height=600
        )
        return fig

    def plot_dendrogram(self, method='ward'):
        """
        Generates a hierarchical clustering dendrogram of the samples.
        """
        scaled_data, samples = self._prepare_data_for_clustering()
        
        # Generate linkage matrix
        Z = sch.linkage(scaled_data, method=method)
        
        fig, ax = plt.subplots(figsize=(10, max(6, len(samples) * 0.3)))
        sch.dendrogram(
            Z,
            labels=samples.tolist(),
            orientation='right',
            leaf_font_size=10
        )
        plt.title('Hierarchical Clustering Dendrogram', fontsize=16)
        plt.xlabel('Distance')
        plt.tight_layout()
        
        buf = BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        plt.close(fig)
        return buf
