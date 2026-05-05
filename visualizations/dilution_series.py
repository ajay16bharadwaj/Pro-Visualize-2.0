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
    Generates standard plots for analyzing proteomics dilution series data.
    """

    def __init__(self, protein_df, metadata_df,
                 protein_id_col='Protein.Group', gene_col='Genes',
                 sample_col='Sample', concentration_col='Concentration',
                 replicate_col='Replicate', group_col='Group'):
        logger.info("Initializing DilutionSeriesVisualizer...")
        if not isinstance(protein_df, pd.DataFrame) or not isinstance(metadata_df, pd.DataFrame):
            raise ValueError("protein_df and metadata_df must be pandas DataFrames.")

        self.protein_df = protein_df.copy()
        self.metadata_df = metadata_df.copy()
        self.protein_id_col = protein_id_col
        self.gene_col = gene_col if gene_col in protein_df.columns else None

        # Validate user-specified metadata columns exist before renaming
        user_meta_cols = [sample_col, concentration_col, replicate_col, group_col]
        missing = [c for c in user_meta_cols if c not in self.metadata_df.columns]
        if missing:
            raise ValueError(f"metadata_df is missing columns: {missing}")
        if self.protein_id_col not in self.protein_df.columns:
            raise ValueError(f"Protein ID column '{self.protein_id_col}' not found in protein data.")

        # Normalize metadata to internal column names so all downstream code is consistent
        internal_names = {'Sample': sample_col, 'Concentration': concentration_col,
                          'Replicate': replicate_col, 'Group': group_col}
        col_renames = {v: k for k, v in internal_names.items() if v != k and v in self.metadata_df.columns}
        if col_renames:
            self.metadata_df = self.metadata_df.rename(columns=col_renames)

        # Validate sample alignment
        self.sample_cols = self.metadata_df['Sample'].unique().tolist()
        missing_samples = [s for s in self.sample_cols if s not in self.protein_df.columns]
        if missing_samples:
            raise ValueError(f"Sample columns missing in protein_df: {missing_samples}")

        self.annotation_cols = [self.protein_id_col]
        if self.gene_col:
            self.annotation_cols.append(self.gene_col)

        self.group_order = self.metadata_df.sort_values(by='Concentration')['Group'].unique().tolist()

        # Internal caches
        self._log2_long_df = None
        self._mean_log2_stats_df = None
        self._cv_stats_df = None
        self._r2_df = None

        logger.info("Initialization complete.")

    # ---------------------------------------------------------------------------
    # Private helpers
    # ---------------------------------------------------------------------------

    def _get_long_data(self, log_transform=True):
        if self._log2_long_df is None:
            df_to_melt = self.protein_df[self.annotation_cols + self.sample_cols].copy()
            intensity_cols = df_to_melt.columns.difference(self.annotation_cols)
            df_numeric = df_to_melt[intensity_cols].replace(0, np.nan)
            df_log2 = np.log2(df_numeric)
            df_to_melt = pd.concat([df_to_melt[self.annotation_cols], df_log2], axis=1)
            df_long = pd.melt(df_to_melt, id_vars=self.annotation_cols, var_name='Sample', value_name='Log2Intensity')
            df_long.dropna(subset=['Log2Intensity'], inplace=True)
            meta_to_merge = self.metadata_df[['Sample', 'Group', 'Concentration']].drop_duplicates()
            self._log2_long_df = pd.merge(df_long, meta_to_merge, on='Sample', how='left')
        return self._log2_long_df

    def _get_cv_stats(self):
        if self._cv_stats_df is None:
            df_long_raw = pd.melt(self.protein_df, id_vars=self.annotation_cols, var_name='Sample', value_name='Intensity')
            df_merged = pd.merge(df_long_raw, self.metadata_df, on='Sample', how='left')
            cv_stats = df_merged.groupby([self.protein_id_col, 'Group'])['Intensity'].agg(
                Mean='mean', StdDev='std').reset_index()
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

    def _classify_deviation_color(self, deviation: float,
                                   good_thresh: float = 0.2,
                                   warn_thresh: float = 0.5) -> str:
        abs_dev = abs(deviation)
        if abs_dev < good_thresh:
            return '#2ca02c'
        elif abs_dev < warn_thresh:
            return '#ff7f0e'
        else:
            return '#d62728'

    # ---------------------------------------------------------------------------
    # Sanity checks
    # ---------------------------------------------------------------------------

    def run_sanity_checks(self) -> list[dict]:
        """
        Returns a list of {level, message} dicts describing data quality issues.
        Levels: 'error', 'warning', 'info'.
        """
        issues = []
        meta = self.metadata_df

        # 1. Non-positive concentrations
        if (meta['Concentration'] <= 0).any():
            n = (meta['Concentration'] <= 0).sum()
            issues.append({"level": "error",
                           "message": f"{n} rows in metadata have non-positive Concentration values. "
                                      "Log-scale fitting will fail."})

        # 2. Duplicate concentration-replicate pairs (pipetting error signal)
        dup_check = meta.groupby(['Replicate', 'Concentration']).size()
        dups = dup_check[dup_check > 1]
        if not dups.empty:
            pairs = [f"Replicate '{r}' at {c}" for (r, c) in dups.index]
            issues.append({"level": "warning",
                           "message": f"Duplicate Concentration×Replicate pairs detected: {', '.join(pairs[:5])}. "
                                      "This may indicate metadata errors."})

        # 3. Concentration series looks non-geometric (warn if ratios vary > 2×)
        conc_vals = sorted(meta['Concentration'].unique())
        if len(conc_vals) >= 2:
            ratios = [conc_vals[i + 1] / conc_vals[i] for i in range(len(conc_vals) - 1)]
            if max(ratios) / min(ratios) > 2.0:
                issues.append({"level": "warning",
                               "message": "Concentration ratios between steps are uneven "
                                          f"(min ratio {min(ratios):.2f}×, max {max(ratios):.2f}×). "
                                          "Expected lines in the ratio plot assume a geometric series."})

        # 4. Negative slope proteins (requires R² computation)
        try:
            r2_df = self.get_r2_table()
            neg = (r2_df['slope'] < 0).sum()
            if neg > 0:
                pct = 100.0 * neg / len(r2_df)
                issues.append({"level": "warning",
                               "message": f"{neg} of {len(r2_df)} proteins ({pct:.1f}%) have a negative "
                                          "linear fit slope — intensity decreases with concentration. "
                                          "LOD/LOQ estimates will be unreliable for these proteins."})
            low_r2 = (r2_df['r_squared'] < 0.8).sum()
            if low_r2 > 0:
                pct = 100.0 * low_r2 / len(r2_df)
                issues.append({"level": "info",
                               "message": f"{low_r2} proteins ({pct:.1f}%) have R² < 0.8 — "
                                          "poor linearity across the dilution series."})
        except Exception as e:
            logger.warning(f"Sanity check skipped linearity step: {e}")

        return issues

    # ---------------------------------------------------------------------------
    # Data export helpers
    # ---------------------------------------------------------------------------

    def get_r2_table(self) -> pd.DataFrame:
        """Per-protein linear fit statistics (log2 intensity ~ log2 concentration)."""
        if self._r2_df is not None:
            return self._r2_df

        mean_stats = self._get_mean_log2_stats()
        rows = []
        for protein, grp in mean_stats.groupby(self.protein_id_col):
            clean = grp.dropna(subset=['Log2Concentration', 'Log2Intensity'])
            if len(clean) < 3:
                continue
            try:
                slope, intercept, r_val, p_val, std_err = linregress(
                    clean['Log2Concentration'], clean['Log2Intensity'])
                predicted = slope * clean['Log2Concentration'] + intercept
                residuals = clean['Log2Intensity'] - predicted
                residual_std = residuals.std()
                row = {
                    self.protein_id_col: protein,
                    'slope': round(slope, 4),
                    'intercept': round(intercept, 4),
                    'r_squared': round(r_val ** 2, 4),
                    'p_value': round(p_val, 6),
                    'residual_std': round(residual_std, 4),
                    'n_points': len(clean),
                }
                if self.gene_col and self.gene_col in clean.columns:
                    row[self.gene_col] = clean[self.gene_col].iloc[0]
                rows.append(row)
            except Exception:
                pass

        self._r2_df = pd.DataFrame(rows)
        return self._r2_df

    def get_cv_by_concentration_matrix(self) -> pd.DataFrame:
        """CV% per protein × concentration group, ready for CSV export."""
        cv_stats = self._get_cv_stats()
        matrix = cv_stats.pivot_table(
            index=self.protein_id_col, columns='Group',
            values='CV_Percent', aggfunc='mean')
        col_order = [g for g in self.group_order if g in matrix.columns]
        return matrix[col_order].round(2)

    def get_completeness_summary(self, cv_threshold: float = 20.0) -> pd.DataFrame:
        """Group-level detection completeness table, ready for CSV export."""
        rows = []
        for group in self.group_order:
            grp_samples = self.metadata_df[self.metadata_df['Group'] == group]['Sample'].tolist()
            grp_samples = [s for s in grp_samples if s in self.protein_df.columns]
            if not grp_samples:
                continue
            grp_data = self.protein_df[[self.protein_id_col] + grp_samples].copy()
            for col in grp_samples:
                grp_data[col] = grp_data[col].replace(0, np.nan)
            total = grp_data[self.protein_id_col][grp_data[grp_samples].notna().any(axis=1)].nunique()
            avg = grp_data[grp_samples].notna().sum().mean()
            complete_mask = grp_data[grp_samples].notna().all(axis=1)
            complete = int(complete_mask.sum())
            complete_data = grp_data[complete_mask].copy()
            if len(complete_data) > 0:
                complete_data['_mean'] = complete_data[grp_samples].mean(axis=1)
                complete_data['_std'] = complete_data[grp_samples].std(axis=1)
                complete_data['_cv'] = complete_data['_std'] / complete_data['_mean'] * 100
                high_quality = int((complete_data['_cv'] < cv_threshold).sum())
            else:
                high_quality = 0
            conc = self.metadata_df[self.metadata_df['Group'] == group]['Concentration'].iloc[0]
            rows.append({
                'Group': group, 'Concentration': conc,
                'Total Detected': int(total), 'Avg per Sample': round(avg, 1),
                'Complete (All Reps)': complete,
                f'High Quality (CV<{cv_threshold}%)': high_quality,
            })
        return pd.DataFrame(rows)

    def get_lod_loq_table(self) -> pd.DataFrame:
        """
        LOD/LOQ per protein using slope-based estimation from the log-log linear fit.

        LOD = 2^(x_min + 3.3 * σ / slope)
        LOQ = 2^(x_min + 10.0 * σ / slope)

        where x_min = log2(minimum tested concentration), σ = residual std of the fit.
        Results with slope ≤ 0 are set to NaN (physically undefined).
        """
        r2_df = self.get_r2_table()
        if r2_df.empty:
            return pd.DataFrame()

        x_min = np.log2(self.metadata_df['Concentration'].min())
        x_max = np.log2(self.metadata_df['Concentration'].max())

        rows = []
        for _, row in r2_df.iterrows():
            slope = row['slope']
            sigma = row['residual_std']
            if slope > 0:
                lod_log2_conc = x_min + (3.3 * sigma / slope)
                loq_log2_conc = x_min + (10.0 * sigma / slope)
                lod = round(2 ** lod_log2_conc, 4)
                loq = round(2 ** loq_log2_conc, 4)
                in_range = x_min <= lod_log2_conc <= x_max
            else:
                lod = loq = float('nan')
                in_range = False

            entry = {
                self.protein_id_col: row[self.protein_id_col],
                'R²': row['r_squared'],
                'slope': row['slope'],
                'residual_std': row['residual_std'],
                'LOD (ng)': lod,
                'LOQ (ng)': loq,
                'LOD in range': in_range,
            }
            if self.gene_col and self.gene_col in row.index:
                entry[self.gene_col] = row[self.gene_col]
            rows.append(entry)

        return pd.DataFrame(rows).sort_values('R²', ascending=False)

    # ---------------------------------------------------------------------------
    # Plot methods — existing
    # ---------------------------------------------------------------------------

    def plot_intensity_distribution(self, plot_type='box', **kwargs):
        plot_data = self._get_long_data(log_transform=True)
        plot_args = {
            "x": 'Group', "y": 'Log2Intensity', "color": 'Group',
            "title": "Log2 Intensity Distribution by Group",
            "category_orders": {'Group': self.group_order},
            "labels": {'Group': 'Concentration Group', 'Log2Intensity': 'Log2(Intensity)'}
        }
        plot_args.update(kwargs)
        if plot_type == 'violin':
            fig = px.violin(plot_data, **plot_args, box=True)
        else:
            fig = px.box(plot_data, **plot_args)
        fig.update_layout(showlegend=False)
        return fig

    def plot_protein_trends(self, proteins_to_plot=None, n_top_proteins=5, **kwargs):
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

        for trace in fig.data:
            trace_data = plot_data[plot_data[color_col] == trace.name]
            clean_data = trace_data.dropna(subset=['Log2Concentration', 'Log2Intensity'])
            if len(clean_data) >= 2:
                try:
                    slope, intercept, r_value, _, _ = linregress(
                        clean_data['Log2Concentration'], clean_data['Log2Intensity'])
                    trace.name = f"{trace.name} (R²={r_value ** 2:.3f})"
                except Exception:
                    pass

        conc_values = sorted(plot_data['Concentration'].unique())
        log2_conc_values = np.log2(conc_values)
        fig.update_layout(
            title="Protein Intensity Trend across Dilution Series",
            xaxis=dict(tickmode='array', tickvals=log2_conc_values, ticktext=[str(c) for c in conc_values]),
            xaxis_title="Concentration (ng)",
            yaxis_title="Mean Log2(Intensity)",
            template=template
        )
        return fig

    def plot_heatmap_trends(self, min_concentrations_present=4, max_proteins_to_plot=500, apply_zscore=True, **kwargs):
        mean_intensity_stats = self._get_mean_log2_stats()
        heatmap_matrix = mean_intensity_stats.pivot(
            index=self.protein_id_col, columns='Concentration', values='Log2Intensity').sort_index(axis=1)
        heatmap_filtered = heatmap_matrix.dropna(thresh=min_concentrations_present)
        if len(heatmap_filtered) > max_proteins_to_plot:
            row_variances = heatmap_filtered.var(axis=1, skipna=True)
            heatmap_filtered = heatmap_filtered.loc[row_variances.nlargest(max_proteins_to_plot).index]

        if apply_zscore:
            matrix_for_plot = heatmap_filtered.apply(lambda x: zscore(x.dropna()), axis=1, result_type='expand').fillna(0)
            color_axis_label, color_scale = "Z-Score (Log2 Intensity)", 'RdBu_r'
        else:
            matrix_for_plot = heatmap_filtered.fillna(heatmap_filtered.min().min())
            color_axis_label, color_scale = "Mean Log2 Intensity", 'Viridis'

        fig = px.imshow(matrix_for_plot, aspect="auto",
                        labels=dict(x="Concentration", y="Protein", color=color_axis_label),
                        title="Heatmap of Protein Intensity Trends",
                        color_continuous_scale=color_scale, **kwargs)
        fig.update_layout(yaxis={'visible': False, 'showticklabels': False})
        fig.update_xaxes(side="bottom", type='category')
        return fig

    def plot_cv_distribution(self, y_limit_percentile=98.0, **kwargs):
        cv_stats_plot = self._get_cv_stats()
        fig = px.box(cv_stats_plot, x='Group', y='CV_Percent', color='Group',
                     title="Protein CV% Distribution",
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
        fig = px.bar(plot_data, x='Sample', y='ProteinCount', color='Group',
                     title="Quantified Proteins per Sample",
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
        pca_result_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'],
                                     index=df_pca_input.index).reset_index().rename(columns={'index': 'Sample'})
        pca_plot_df = pd.merge(pca_result_df, self.metadata_df, on='Sample', how='left')
        fig = px.scatter(pca_plot_df, x='PC1', y='PC2', color=color_by, symbol=symbol_by,
                         title="PCA of Samples",
                         labels={'PC1': f"PC1 ({explained_variance[0]:.1f}%)",
                                 'PC2': f"PC2 ({explained_variance[1]:.1f}%)"}, **kwargs)
        return fig

    def plot_relative_abundance_ratios(self, add_expected_lines=True, show_deviations=True,
                                       good_thresh=0.2, warn_thresh=0.5, **kwargs):
        """
        Box plot of protein log2 intensity ratios vs the lowest concentration,
        with optional deviation overlay color-coded by configurable thresholds.
        """
        logger.info("Generating Log2 Relative Abundance Ratio plot...")
        mean_log2_stats = self._get_mean_log2_stats()
        if mean_log2_stats is None or mean_log2_stats.empty:
            raise ValueError("Mean log2 intensity data is not available.")

        min_concentration = self.metadata_df['Concentration'].min()
        if min_concentration <= 0:
            raise ValueError("Minimum concentration must be positive.")

        base_log2_map = (mean_log2_stats[mean_log2_stats['Concentration'] == min_concentration]
                         .groupby(self.protein_id_col)['Log2Intensity'].first().dropna())
        ratio_data = mean_log2_stats.copy()
        ratio_data['BaseMeanLog2Intensity'] = ratio_data[self.protein_id_col].map(base_log2_map)
        ratio_data.dropna(subset=['BaseMeanLog2Intensity', 'Log2Intensity'], inplace=True)
        ratio_data = ratio_data[ratio_data['Concentration'] > min_concentration].copy()
        if ratio_data.empty:
            raise ValueError("No proteins found with valid base intensities to calculate ratios.")

        ratio_data['Log2Ratio'] = ratio_data['Log2Intensity'] - ratio_data['BaseMeanLog2Intensity']
        group_info = self.metadata_df[['Concentration', 'Group']].drop_duplicates()
        ratio_data = pd.merge(ratio_data, group_info, on='Concentration', how='left').dropna(subset=['Group'])
        ratio_data['ExpectedLog2Ratio'] = np.log2(ratio_data['Concentration'] / min_concentration)

        base_group_name = self.metadata_df.loc[
            self.metadata_df['Concentration'] == min_concentration, 'Group'].unique()[0]
        groups_to_plot_order = [g for g in self.group_order
                                if g != base_group_name and g in ratio_data['Group'].unique()]

        fig = px.box(ratio_data, x='Group', y='Log2Ratio', color='Group',
                     title="Protein Log2 Abundance Ratio vs. Lowest Concentration",
                     category_orders={'Group': groups_to_plot_order},
                     labels={'Group': f'Group (Ratio to {base_group_name})',
                             'Log2Ratio': f'Log2(Intensity Ratio vs {min_concentration}ng)'},
                     **kwargs)

        if add_expected_lines:
            expected_ratios_df = ratio_data[['Group', 'ExpectedLog2Ratio']].drop_duplicates()
            for _, row in expected_ratios_df.iterrows():
                fig.add_hline(y=row['ExpectedLog2Ratio'], line_dash="dash", line_color="grey",
                              annotation_text=f"Exp: {row['ExpectedLog2Ratio']:.2f}",
                              annotation_position="bottom right")

        group_stats = None
        if show_deviations:
            ratio_data['Deviation'] = ratio_data['Log2Ratio'] - ratio_data['ExpectedLog2Ratio']
            group_stats = ratio_data.groupby('Group').agg(
                MeanDeviation=('Deviation', 'mean'),
                StdDeviation=('Log2Ratio', 'std'),
                Count=('Log2Ratio', 'count'),
                MeanObserved=('Log2Ratio', 'mean'),
                Expected=('ExpectedLog2Ratio', 'first')
            ).reindex(groups_to_plot_order).reset_index()
            group_stats['Color'] = group_stats['MeanDeviation'].apply(
                lambda d: self._classify_deviation_color(d, good_thresh, warn_thresh))

        fig.add_hline(y=0, line_dash="solid", line_color="black", line_width=1)

        if show_deviations and group_stats is not None:
            fig.add_trace(go.Scatter(
                x=group_stats['Group'], y=group_stats['MeanObserved'], mode='markers',
                marker=dict(size=14, symbol='diamond', color=group_stats['Color'],
                            line=dict(width=2, color='black')),
                error_y=dict(type='data', array=group_stats['StdDeviation'], visible=True,
                             color='rgba(0,0,0,0.5)', thickness=2, width=6),
                name='Mean Observed ± SD', showlegend=True,
                hovertemplate=(
                    '<b>%{x}</b><br>Mean Observed: %{y:.3f}<br>'
                    'Bias (Mean Dev): %{customdata[0]:+.3f}<br>SD: %{customdata[1]:.3f}<br><extra></extra>'
                ),
                customdata=group_stats[['MeanDeviation', 'StdDeviation']].values
            ))
            legend_group = 'Deviation Classification'
            for label, color in [
                (f'Low Bias (|dev| < {good_thresh})', '#2ca02c'),
                (f'Moderate Bias (|dev| < {warn_thresh})', '#ff7f0e'),
                (f'High Bias (|dev| ≥ {warn_thresh})', '#d62728'),
            ]:
                fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', name=label,
                                         marker=dict(size=10, symbol='diamond', color=color,
                                                     line=dict(width=2, color='black')),
                                         legendgroup=legend_group, showlegend=True))

        fig.update_layout(showlegend=True)
        return fig

    def plot_completeness_overview(self, identifier_col=None, use_log_scale=False,
                                   cv_threshold=20.0, **kwargs):
        """Stacked bar charts showing detection completeness."""
        if identifier_col is None:
            identifier_col = self.protein_id_col
        if identifier_col not in self.protein_df.columns:
            raise ValueError(f"Column '{identifier_col}' not found in protein data.")

        template = kwargs.get('template', 'plotly_white')
        groups = self.group_order
        total_list, avg_list, complete_list, high_quality_list, group_labels = [], [], [], [], []

        for group in groups:
            group_samples = self.metadata_df[self.metadata_df['Group'] == group]['Sample'].tolist()
            group_samples = [s for s in group_samples if s in self.protein_df.columns]
            if not group_samples:
                continue
            group_data = self.protein_df[[identifier_col] + group_samples].copy()
            for col in group_samples:
                group_data[col] = group_data[col].replace(0, np.nan)
            total_detected = group_data[identifier_col][group_data[group_samples].notna().any(axis=1)].nunique()
            avg_per_sample = group_data[group_samples].notna().sum().mean()
            complete_detection_mask = group_data[group_samples].notna().all(axis=1)
            complete_detection_count = complete_detection_mask.sum()
            complete_data = group_data[complete_detection_mask].copy()
            if len(complete_data) > 0:
                complete_data['mean_intensity'] = complete_data[group_samples].mean(axis=1)
                complete_data['std_intensity'] = complete_data[group_samples].std(axis=1)
                complete_data['cv_percent'] = (complete_data['std_intensity'] / complete_data['mean_intensity']) * 100
                high_quality_count = (complete_data['cv_percent'] < cv_threshold).sum()
            else:
                high_quality_count = 0
            group_labels.append(group)
            total_list.append(int(total_detected))
            avg_list.append(int(round(avg_per_sample)))
            complete_list.append(int(complete_detection_count))
            high_quality_list.append(int(high_quality_count))

        total_arr = np.array(total_list)
        avg_arr = np.array(avg_list)
        complete_arr = np.array(complete_list)
        high_quality_arr = np.array(high_quality_list)

        slice_green = high_quality_arr
        slice_blue = np.maximum(0, complete_arr - high_quality_arr)
        slice_pink = np.maximum(0, avg_arr - complete_arr)
        slice_teal = np.maximum(0, total_arr - avg_arr)

        with np.errstate(divide='ignore', invalid='ignore'):
            slice_green_pct = np.nan_to_num(slice_green / total_arr) * 100
            slice_blue_pct = np.nan_to_num(slice_blue / total_arr) * 100
            slice_pink_pct = np.nan_to_num(slice_pink / total_arr) * 100
            slice_teal_pct = np.nan_to_num(slice_teal / total_arr) * 100

        def get_clean_labels(arr):
            return [str(int(x)) if x > 0 else "" for x in arr]

        plot_name = identifier_col.replace('.', ' ').replace('_', ' ').title()
        fig = make_subplots(rows=1, cols=2,
                            subplot_titles=(f"Total {plot_name}s (Absolute)", f"% {plot_name}s (to total)"),
                            horizontal_spacing=0.1)

        bar_specs = [
            (f'All Reps + CV<{cv_threshold}%', slice_green, '#2ca02c', 'High Quality', True),
            ('All Reps', slice_blue, '#1f77b4', 'Complete', True),
            ('Average/Sample', slice_pink, '#ff7f0e', 'Average', True),
            ('Total', slice_teal, '#d62728', 'Total', True),
        ]
        for name, y, color, leg_grp, show_leg in bar_specs:
            fig.add_trace(go.Bar(name=name, x=group_labels, y=y, marker_color=color,
                                 text=get_clean_labels(y), textposition='inside',
                                 legendgroup=leg_grp, showlegend=show_leg), row=1, col=1)

        pct_specs = [
            (f'All Reps + CV<{cv_threshold}%', slice_green_pct, '#2ca02c', 'High Quality'),
            ('All Reps', slice_blue_pct, '#1f77b4', 'Complete'),
            ('Average/Sample', slice_pink_pct, '#ff7f0e', 'Average'),
            ('Total', slice_teal_pct, '#d62728', 'Total'),
        ]
        for name, y, color, leg_grp in pct_specs:
            fig.add_trace(go.Bar(name=name, x=group_labels, y=y, marker_color=color,
                                 text=y, texttemplate='%{text:.0f}%', textposition='inside',
                                 textfont=dict(size=10, color='white'),
                                 legendgroup=leg_grp, showlegend=False), row=1, col=2)

        fig.update_layout(
            height=500, template=template, barmode='stack',
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            yaxis=dict(title='Absolute Count', type='log' if use_log_scale else 'linear'),
            xaxis=dict(title='Concentration Group'),
            yaxis2=dict(title='% of Total', range=[0, 101]),
            xaxis2=dict(title='Concentration Group')
        )
        return fig

    # ---------------------------------------------------------------------------
    # New plot methods — P4
    # ---------------------------------------------------------------------------

    def plot_r2_histogram(self, n_bins=20, **kwargs):
        """Histogram of per-protein R² values from the log-log linear fit."""
        r2_df = self.get_r2_table()
        if r2_df.empty:
            raise ValueError("No proteins with sufficient data points for R² calculation.")

        fig = px.histogram(r2_df, x='r_squared', nbins=n_bins,
                           title=f"Distribution of Linearity (R²) — {len(r2_df)} proteins",
                           labels={'r_squared': 'R² (Log2 Intensity ~ Log2 Concentration)',
                                   'count': 'Number of Proteins'},
                           **kwargs)
        fig.add_vline(x=0.8, line_dash="dash", line_color="orange",
                      annotation_text="R²=0.8", annotation_position="top right")
        fig.add_vline(x=0.95, line_dash="dash", line_color="green",
                      annotation_text="R²=0.95", annotation_position="top left")
        median_r2 = r2_df['r_squared'].median()
        fig.update_layout(
            annotations=[
                dict(x=0.02, y=0.97, xref="paper", yref="paper", showarrow=False,
                     text=f"Median R² = {median_r2:.3f}", bgcolor="white",
                     bordercolor="grey", borderwidth=1)
            ]
        )
        return fig

    def plot_lod_loq(self, top_n=50, **kwargs):
        """
        Scatter plot of LOD vs R² per protein, with LOQ shown as hover text.

        Only proteins with valid (positive slope, sufficient data) fits are shown.
        Vertical lines mark the minimum and maximum tested concentrations.
        """
        lod_loq_df = self.get_lod_loq_table()
        if lod_loq_df.empty:
            raise ValueError("Could not compute LOD/LOQ — no proteins with valid fits.")

        valid = lod_loq_df.dropna(subset=['LOD (ng)', 'LOQ (ng)'])
        if valid.empty:
            raise ValueError("All proteins have non-positive slope — LOD/LOQ undefined.")

        hover_col = self.gene_col if self.gene_col and self.gene_col in valid.columns else self.protein_id_col
        min_conc = self.metadata_df['Concentration'].min()
        max_conc = self.metadata_df['Concentration'].max()

        # Show top_n by R² for readability
        plot_data = valid.nlargest(min(top_n, len(valid)), 'R²').copy()
        plot_data['LOD (ng)'] = plot_data['LOD (ng)'].clip(lower=min_conc * 0.01)

        fig = px.scatter(
            plot_data, x='LOD (ng)', y='R²', color='slope',
            hover_name=hover_col,
            hover_data={self.protein_id_col: True, 'LOQ (ng)': ':.4f',
                        'slope': ':.3f', 'residual_std': ':.3f', 'LOD in range': True},
            color_continuous_scale='RdYlGn',
            title=f"LOD vs Linearity (R²) — top {len(plot_data)} proteins by R²",
            labels={'LOD (ng)': 'LOD (ng, log scale)', 'R²': 'R² of Log-Log Fit',
                    'slope': 'Slope'},
            **kwargs
        )
        fig.update_xaxes(type='log')
        fig.add_vline(x=min_conc, line_dash="dash", line_color="blue",
                      annotation_text=f"Min tested ({min_conc} ng)", annotation_position="top right")
        fig.add_vline(x=max_conc, line_dash="dash", line_color="red",
                      annotation_text=f"Max tested ({max_conc} ng)", annotation_position="top left")
        fig.update_layout(coloraxis_colorbar=dict(title="Slope"))
        return fig
