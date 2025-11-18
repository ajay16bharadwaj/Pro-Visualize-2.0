import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import logging
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import zscore, linregress  

logger = logging.getLogger(__name__)

class DilutionSeriesVisualizer:
    """
    A class to generate standard plots for analyzing proteomics dilution series data.
    """

    def __init__(self, protein_df, metadata_df,
                 protein_id_col='Protein.Group', gene_col='Genes'):
        logger.info("Initializing DilutionSeriesVisualizer...")
        if not isinstance(protein_df, pd.DataFrame) or not isinstance(metadata_df, pd.DataFrame):
            raise ValueError("protein_df and metadata_df must be pandas DataFrames.")

        self.protein_df = protein_df.copy()
        self.metadata_df = metadata_df.copy()
        self.protein_id_col = protein_id_col
        self.gene_col = gene_col if gene_col in self.protein_df.columns else None
        
        # --- Validate Inputs ---
        required_meta_cols = ['Sample', 'Group', 'Concentration', 'Replicate']
        if not all(col in self.metadata_df.columns for col in required_meta_cols):
            raise ValueError(f"metadata_df must contain columns: {required_meta_cols}")
        if self.protein_id_col not in self.protein_df.columns:
            raise ValueError(f"Protein ID column '{self.protein_id_col}' not found.")

        # --- Identify Columns ---
        self.sample_cols = self.metadata_df['Sample'].unique().tolist()
        missing_samples = [s for s in self.sample_cols if s not in self.protein_df.columns]
        if missing_samples:
            raise ValueError(f"Sample columns missing in protein_df: {missing_samples}")

        self.annotation_cols = [self.protein_id_col]
        if self.gene_col:
            self.annotation_cols.append(self.gene_col)
        
        # --- Pre-calculate Group Order ---
        self.group_order = self.metadata_df.sort_values(by='Concentration')['Group'].unique().tolist()
        
        # --- Internal Cache ---
        self._log2_long_df = None
        self._mean_log2_stats_df = None
        self._cv_stats_df = None
        
        logger.info("Initialization complete.")

    # --- Private Helper Methods for Data Preparation ---
    def _get_long_data(self, log_transform=True):
        cache_attr = '_log2_long_df'
        value_col_name = 'Log2Intensity'

        if getattr(self, cache_attr) is None:
            logger.info(f"Preparing {value_col_name} long format data...")
            df_to_melt = self.protein_df[self.annotation_cols + self.sample_cols].copy()
            intensity_cols = df_to_melt.columns.difference(self.annotation_cols)
            df_numeric = df_to_melt[intensity_cols].replace(0, np.nan)
            df_log2 = np.log2(df_numeric)
            df_to_melt = pd.concat([df_to_melt[self.annotation_cols], df_log2], axis=1)

            df_long = pd.melt(df_to_melt, id_vars=self.annotation_cols, var_name='Sample', value_name=value_col_name)
            df_long.dropna(subset=[value_col_name], inplace=True)
            
            meta_to_merge = self.metadata_df[['Sample', 'Group', 'Concentration']].drop_duplicates()
            df_merged = pd.merge(df_long, meta_to_merge, on='Sample', how='left')
            setattr(self, cache_attr, df_merged)
        return getattr(self, cache_attr)

    def _get_cv_stats(self):
        if self._cv_stats_df is None:
            logger.info("Calculating CV% stats...")
            df_long_raw = pd.melt(self.protein_df, id_vars=self.annotation_cols, var_name='Sample', value_name='Intensity')
            df_merged = pd.merge(df_long_raw, self.metadata_df, on='Sample', how='left')
            
            cv_stats = df_merged.groupby([self.protein_id_col, 'Group'])['Intensity'].agg(Mean='mean', StdDev='std').reset_index()
            cv_stats['CV_Percent'] = (cv_stats['StdDev'] / cv_stats['Mean']) * 100
            
            group_to_conc = self.metadata_df[['Group', 'Concentration']].drop_duplicates().set_index('Group')
            self._cv_stats_df = cv_stats.join(group_to_conc, on='Group')
        return self._cv_stats_df
    
    def _get_mean_log2_stats(self):
        if self._mean_log2_stats_df is None:
            log2_long_data = self._get_long_data(log_transform=True)
            mean_stats = log2_long_data.groupby([self.protein_id_col, 'Concentration'])['Log2Intensity'].mean().reset_index()
            if self.gene_col:
                annotations = self.protein_df[[self.protein_id_col, self.gene_col]].drop_duplicates()
                mean_stats = pd.merge(mean_stats, annotations, on=self.protein_id_col, how='left')
            mean_stats = mean_stats[mean_stats['Concentration'] > 0]
            mean_stats['Log2Concentration'] = np.log2(mean_stats['Concentration'])
            self._mean_log2_stats_df = mean_stats
        return self._mean_log2_stats_df

    def plot_intensity_distribution(self, plot_type='box', **kwargs):
        """Generates Plot 1: Intensity Distribution."""
        plot_data = self._get_long_data(log_transform=True)
        
        plot_args = {
            "x": 'Group', "y": 'Log2Intensity', "color": 'Group',
            "title": "Log2 Intensity Distribution by Group",
            "category_orders": {'Group': self.group_order},
            "labels": {'Group': 'Concentration Group', 'Log2Intensity': 'Log2(Intensity)'}
        }
        # Update with user kwargs (template, color_discrete_map, etc.)
        plot_args.update(kwargs)

        if plot_type == 'violin':
            # Explicitly add box=True for violin plots
            fig = px.violin(plot_data, **plot_args, box=True)
        else:
            fig = px.box(plot_data, **plot_args)
            
        fig.update_layout(showlegend=False)
        return fig

    def plot_protein_trends(self, proteins_to_plot=None, n_top_proteins=5, **kwargs):
        """Generates Plot 3: Individual Protein Trends with R² values."""
        # Extract template at the start
        template = kwargs.get('template', 'plotly_white')
        
        mean_intensity_stats = self._get_mean_log2_stats()
        if proteins_to_plot:
            plot_data = mean_intensity_stats[mean_intensity_stats[self.protein_id_col].isin(proteins_to_plot)]
        else:
            protein_avg_intensity = mean_intensity_stats.groupby(self.protein_id_col)['Log2Intensity'].mean()
            top_protein_ids = protein_avg_intensity.nlargest(n_top_proteins).index.tolist()
            plot_data = mean_intensity_stats[mean_intensity_stats[self.protein_id_col].isin(top_protein_ids)]

        color_col = self.gene_col if self.gene_col and not plot_data[self.gene_col].isnull().all() else self.protein_id_col
        
        fig = px.line(plot_data, x='Log2Concentration', y='Log2Intensity', color=color_col, markers=True, **kwargs)
        
        # Calculate R2 and update legend names
        for trace in fig.data:
            trace_data = plot_data[plot_data[color_col] == trace.name]
            clean_data = trace_data.dropna(subset=['Log2Concentration', 'Log2Intensity'])
            
            if len(clean_data) >= 2:  # Consider changing to >= 3 for more reliable R²
                try:
                    slope, intercept, r_value, p_value, std_err = linregress(
                        clean_data['Log2Concentration'], 
                        clean_data['Log2Intensity']
                    )
                    r_squared = r_value ** 2
                    trace.name = f"{trace.name} (R²={r_squared:.3f})"
                except Exception:
                    pass

        conc_values = sorted(plot_data['Concentration'].unique())
        log2_conc_values = np.log2(conc_values)
        
        # ✅ FIX: Add template parameter
        fig.update_layout(
            title="Protein Intensity Trend across Dilution Series",
            xaxis=dict(tickmode='array', tickvals=log2_conc_values, ticktext=[str(c) for c in conc_values]),
            xaxis_title="Concentration (ng)",
            yaxis_title="Mean Log2(Intensity)",
            template=template  # ← ADD THIS
        )
        return fig

    def plot_heatmap_trends(self, min_concentrations_present=4, max_proteins_to_plot=500, apply_zscore=True, **kwargs):
        """Generates Plot 4: Heatmap of Trends."""
        mean_intensity_stats = self._get_mean_log2_stats()
        heatmap_matrix = mean_intensity_stats.pivot(index=self.protein_id_col, columns='Concentration', values='Log2Intensity').sort_index(axis=1)
        heatmap_filtered = heatmap_matrix.dropna(thresh=min_concentrations_present)

        if len(heatmap_filtered) > max_proteins_to_plot:
            row_variances = heatmap_filtered.var(axis=1, skipna=True)
            top_variance_proteins = row_variances.nlargest(max_proteins_to_plot).index
            heatmap_filtered = heatmap_filtered.loc[top_variance_proteins]

        if apply_zscore:
            matrix_for_plot = heatmap_filtered.apply(lambda x: zscore(x.dropna()), axis=1, result_type='expand').fillna(0)
            color_axis_label = "Z-Score (Log2 Intensity)"
            # Use a diverging blue-to-red scale for Z-scores
            color_scale = 'RdBu_r'
        else:
            matrix_for_plot = heatmap_filtered.fillna(heatmap_filtered.min().min())
            color_axis_label = "Mean Log2 Intensity"
            # Default color scale for non-normalized data
            color_scale = 'Viridis'

        fig = px.imshow(matrix_for_plot, aspect="auto",
                        labels=dict(x="Concentration", y="Protein", color=color_axis_label),
                        title="Heatmap of Protein Intensity Trends",
                        color_continuous_scale=color_scale,
                        **kwargs)

        fig.update_layout(yaxis={'visible': False, 'showticklabels': False})
        fig.update_xaxes(side="bottom", type='category')
        return fig

    def plot_cv_distribution(self, y_limit_percentile=98.0, **kwargs):
        cv_stats_plot = self._get_cv_stats()
        fig = px.box(cv_stats_plot, x='Group', y='CV_Percent', color='Group', title="Protein CV% Distribution",
                     category_orders={'Group': self.group_order}, 
                     labels={'Group': 'Group', 'CV_Percent': 'CV (%)'}, 
                     **kwargs)
        upper_limit = np.percentile(cv_stats_plot['CV_Percent'].dropna(), y_limit_percentile)
        fig.update_yaxes(range=[0, upper_limit * 1.1])
        fig.update_layout(showlegend=False)
        return fig
        
    def plot_protein_counts_per_sample(self, **kwargs):
        protein_counts = self.protein_df[self.sample_cols].notna().sum().reset_index()
        protein_counts.columns = ['Sample', 'ProteinCount']
        plot_data = pd.merge(protein_counts, self.metadata_df, on='Sample', how='left')
        plot_data = plot_data.sort_values(by=['Concentration', 'Replicate'])
        sample_order = plot_data['Sample'].tolist()
        
        fig = px.bar(plot_data, x='Sample', y='ProteinCount', color='Group', title="Quantified Proteins per Sample",
                     category_orders={'Sample': sample_order}, text='ProteinCount', **kwargs)
        fig.update_traces(textposition='outside')
        return fig

    def plot_pca(self, color_by='Group', symbol_by='Replicate', **kwargs):
        df_log2_wide = np.log2(self.protein_df[self.sample_cols].replace(0, np.nan))
        df_imputed = df_log2_wide.dropna(axis=0)
        df_pca_input = df_imputed.transpose()
        
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df_pca_input)

        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(scaled_data)
        explained_variance = pca.explained_variance_ratio_ * 100

        pca_result_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=df_pca_input.index).reset_index().rename(columns={'index': 'Sample'})
        pca_plot_df = pd.merge(pca_result_df, self.metadata_df, on='Sample', how='left')

        pc1_label = f"PC1 ({explained_variance[0]:.1f}%)"
        pc2_label = f"PC2 ({explained_variance[1]:.1f}%)"

        fig = px.scatter(pca_plot_df, x='PC1', y='PC2', color=color_by, symbol=symbol_by,
                         title="PCA of Samples", labels={'PC1': pc1_label, 'PC2': pc2_label}, **kwargs)
        return fig
    
    def plot_relative_abundance_ratios(self, add_expected_lines=True, **kwargs):
        """
        Generates a box plot of protein log2 intensity ratios relative to the
        mean of the lowest concentration group.
        """
        logger.info("Generating Log2 Relative Abundance Ratio plot...")
        try:
            mean_log2_stats = self._get_mean_log2_stats()
            if mean_log2_stats is None or mean_log2_stats.empty:
                raise ValueError("Mean log2 intensity data is not available.")
            
            min_concentration = self.metadata_df['Concentration'].min()
            if min_concentration <= 0:
                raise ValueError("Minimum concentration must be positive.")

            base_log2_map = mean_log2_stats[mean_log2_stats['Concentration'] == min_concentration].groupby(self.protein_id_col)['Log2Intensity'].first().dropna()
            
            ratio_data = mean_log2_stats.copy()
            ratio_data['BaseMeanLog2Intensity'] = ratio_data[self.protein_id_col].map(base_log2_map)
            ratio_data.dropna(subset=['BaseMeanLog2Intensity', 'Log2Intensity'], inplace=True)
            ratio_data = ratio_data[ratio_data['Concentration'] > min_concentration].copy()

            if ratio_data.empty:
                raise ValueError("No proteins found with valid base intensities to calculate ratios.")

            ratio_data['Log2Ratio'] = ratio_data['Log2Intensity'] - ratio_data['BaseMeanLog2Intensity']
            
            group_info = self.metadata_df[['Concentration', 'Group']].drop_duplicates()
            ratio_data = pd.merge(ratio_data, group_info, on='Concentration', how='left').dropna(subset=['Group'])

            base_group_name = self.metadata_df.loc[self.metadata_df['Concentration'] == min_concentration, 'Group'].unique()[0]
            groups_to_plot_order = [g for g in self.group_order if g != base_group_name and g in ratio_data['Group'].unique()]

            fig = px.box(ratio_data, x='Group', y='Log2Ratio', color='Group',
                         title="Protein Log2 Abundance Ratio vs. Lowest Concentration",
                         category_orders={'Group': groups_to_plot_order},
                         labels={'Group': f'Group (Ratio to {base_group_name})',
                                 'Log2Ratio': f'Log2(Intensity Ratio vs {min_concentration}ng)'},
                         **kwargs)

            if add_expected_lines:
                ratio_data['ExpectedLog2Ratio'] = np.log2(ratio_data['Concentration'] / min_concentration)
                expected_ratios_df = ratio_data[['Group', 'ExpectedLog2Ratio']].drop_duplicates()
                for _, row in expected_ratios_df.iterrows():
                    fig.add_hline(y=row['ExpectedLog2Ratio'], line_dash="dash", line_color="grey",
                                  annotation_text=f"Expected (Log2: {row['ExpectedLog2Ratio']:.2f})",
                                  annotation_position="bottom right")
            
            fig.add_hline(y=0, line_dash="solid", line_color="black", line_width=1)
            fig.update_layout(showlegend=False)
            return fig
        except Exception as e:
            logger.error(f"Error in plot_relative_abundance_ratios: {e}", exc_info=True)
            raise

    def plot_completeness_overview(self, identifier_col='Protein.Group', use_log_scale=False, cv_threshold=20.0, **kwargs):
        """
        Generates stacked bar charts showing detection completeness.
        """
        logger.info(f"Generating completeness overview for {identifier_col}...")
        
        # Validate the identifier column exists
        if identifier_col not in self.protein_df.columns:
            raise ValueError(f"Column '{identifier_col}' not found in protein data.")
        
        # Extract template specifically since go.Figure doesn't use kwargs automatically like px
        template = kwargs.get('template', 'plotly_white')

        # Get the list of groups in order
        groups = self.group_order
        
        # Prepare data structure for stacking
        total_list = []
        avg_list = []
        complete_list = []
        high_quality_list = []
        group_labels = []
        
        # --- 1. Data Calculation Loop ---
        for group in groups:
            group_samples = self.metadata_df[self.metadata_df['Group'] == group]['Sample'].tolist()
            group_samples = [s for s in group_samples if s in self.protein_df.columns]
            
            if not group_samples:
                continue
            
            # Subset data for this group
            group_data = self.protein_df[[identifier_col] + group_samples].copy()
            
            # Replace zeros with NaN for proper counting
            for col in group_samples:
                group_data[col] = group_data[col].replace(0, np.nan)
            
            # Calculate metrics
            total_detected = group_data[identifier_col][
                group_data[group_samples].notna().any(axis=1)
            ].nunique()
            
            avg_per_sample = group_data[group_samples].notna().sum().mean()
            
            complete_detection_mask = group_data[group_samples].notna().all(axis=1)
            complete_detection_count = complete_detection_mask.sum()
            
            complete_data = group_data[complete_detection_mask].copy()
            
            if len(complete_data) > 0:
                complete_data['mean_intensity'] = complete_data[group_samples].mean(axis=1)
                complete_data['std_intensity'] = complete_data[group_samples].std(axis=1)
                complete_data['cv_percent'] = (
                    complete_data['std_intensity'] / complete_data['mean_intensity']
                ) * 100
                
                high_quality_count = (complete_data['cv_percent'] < cv_threshold).sum()
            else:
                high_quality_count = 0
            
            group_labels.append(group)
            total_list.append(int(total_detected))
            avg_list.append(int(round(avg_per_sample)))
            complete_list.append(int(complete_detection_count))
            high_quality_list.append(int(high_quality_count))
        
        # --- 2. Calculate Differential Slices for Stacked Chart ---
        total_arr = np.array(total_list)
        avg_arr = np.array(avg_list)
        complete_arr = np.array(complete_list)
        high_quality_arr = np.array(high_quality_list)

        # Absolute Counts for Left Plot
        slice_green = high_quality_arr
        slice_blue = np.maximum(0, complete_arr - high_quality_arr)
        slice_pink = np.maximum(0, avg_arr - complete_arr)
        slice_teal = np.maximum(0, total_arr - avg_arr)

        # Percentage Values for Right Plot
        # Use numpy to safely divide, avoiding division by zero warnings
        with np.errstate(divide='ignore', invalid='ignore'):
            slice_green_pct = np.nan_to_num(slice_green / total_arr) * 100
            slice_blue_pct = np.nan_to_num(slice_blue / total_arr) * 100
            slice_pink_pct = np.nan_to_num(slice_pink / total_arr) * 100
            slice_teal_pct = np.nan_to_num(slice_teal / total_arr) * 100

        # Helper to clean text labels
        def get_clean_labels(arr):
            return [str(int(x)) if x > 0 else "" for x in arr]
        
        # --- 3. Create the Figure with Subplots ---
        plot_name = identifier_col.replace('.', ' ').replace('_', ' ').title()
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=(f"Total {plot_name}s (Absolute)", f"% {plot_name}s (to total)"),
            horizontal_spacing=0.1
        )
        
        # --- 4. Plot 1: Stacked Bar Chart - ABSOLUTE VALUES (Left) ---
        
        # High Quality (Green)
        fig.add_trace(go.Bar(
            name=f'All Reps + CV<{cv_threshold}%', x=group_labels, y=slice_green,
            marker_color='#2ca02c', hovertemplate='<b>%{x}</b><br>High Quality: %{y}<extra></extra>',
            text=get_clean_labels(slice_green), textposition='inside',
            legendgroup='High Quality', showlegend=True 
        ), row=1, col=1)
        
        # Complete (Blue)
        fig.add_trace(go.Bar(
            name='All Reps', x=group_labels, y=slice_blue,
            marker_color='#1f77b4', hovertemplate='<b>%{x}</b><br>Complete (High CV): %{y}<extra></extra>',
            text=get_clean_labels(slice_blue), textposition='inside',
            legendgroup='Complete', showlegend=True
        ), row=1, col=1)
        
        # Avg per Sample (Pink)
        fig.add_trace(go.Bar(
            name='Average/Sample', x=group_labels, y=slice_pink,
            marker_color='#ff7f0e', hovertemplate='<b>%{x}</b><br>Incomplete (Avg): %{y}<extra></extra>',
            text=get_clean_labels(slice_pink), textposition='inside',
            legendgroup='Average', showlegend=True
        ), row=1, col=1)
        
        # Total Detected (Teal)
        fig.add_trace(go.Bar(
            name='Total', x=group_labels, y=slice_teal,
            marker_color='#d62728', hovertemplate='<b>%{x}</b><br>Incomplete (Low Freq): %{y}<extra></extra>',
            text=get_clean_labels(slice_teal), textposition='inside',
            legendgroup='Total', showlegend=True
        ), row=1, col=1)


        # --- 5. Plot 2: Stacked Bar Chart - PERCENTAGE VALUES (Right) ---
        
        # High Quality (Green) - PCT
        fig.add_trace(go.Bar(
            name=f'All Reps + CV<{cv_threshold}%', x=group_labels, y=slice_green_pct,
            marker_color='#2ca02c', 
            hovertemplate='<b>%{x}</b><br>High Quality: %{y:.1f}%<extra></extra>',
            text=slice_green_pct,
            texttemplate='%{text:.0f}%', textposition='inside', textfont=dict(size=10, color='white'),
            legendgroup='High Quality', showlegend=False
        ), row=1, col=2)
        
        # Complete (Blue) - PCT
        fig.add_trace(go.Bar(
            name='All Reps', x=group_labels, y=slice_blue_pct,
            marker_color='#1f77b4',
            hovertemplate='<b>%{x}</b><br>Complete (High CV): %{y:.1f}%<extra></extra>',
            text=slice_blue_pct,
            texttemplate='%{text:.0f}%', textposition='inside', textfont=dict(size=10, color='white'),
            legendgroup='Complete', showlegend=False
        ), row=1, col=2)
        
        # Avg per Sample (Pink) - PCT
        fig.add_trace(go.Bar(
            name='Average/Sample', x=group_labels, y=slice_pink_pct,
            marker_color='#ff7f0e',
            hovertemplate='<b>%{x}</b><br>Incomplete (Avg): %{y:.1f}%<extra></extra>',
            text=slice_pink_pct,
            texttemplate='%{text:.0f}%', textposition='inside', textfont=dict(size=10, color='white'),
            legendgroup='Average', showlegend=False
        ), row=1, col=2)
        
        # Total Detected (Teal) - PCT
        fig.add_trace(go.Bar(
            name='Total', x=group_labels, y=slice_teal_pct,
            marker_color='#d62728',
            hovertemplate='<b>%{x}</b><br>Incomplete (Low Freq): %{y:.1f}%<extra></extra>',
            text=slice_teal_pct,
            texttemplate='%{text:.0f}%', textposition='inside', textfont=dict(size=10, color='white'),
            legendgroup='Total', showlegend=False
        ), row=1, col=2)

        # --- 6. Update Layouts ---
        fig.update_layout(
            height=500,
            template=template,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            barmode='stack', 
            yaxis=dict(title='Absolute Count', type='log' if use_log_scale else 'linear'),
            xaxis=dict(title='Concentration Group'),
            yaxis2=dict(title='% of Total', range=[0, 101]), 
            xaxis2=dict(title='Concentration Group')
        )

        # Optional: Set a uniform text size for the left plot as well
        fig.for_each_trace(lambda t: t.update(textfont=dict(color='white')), row=1, col=1)
        
        logger.info("Completeness overview plot generated successfully.")
        return fig