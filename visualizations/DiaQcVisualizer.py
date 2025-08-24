import pandas as pd
import logging
from pathlib import Path
from typing import Dict, List
import plotly.graph_objects as go
import plotly.express as px

# Set up a logger for this module
logger = logging.getLogger(__name__)

class DiaQcVisualizer:
    """
    Handles DIA-NN QC data loading, processing, and visualization with an
    intermediate step for metadata editing.
    """

    def __init__(self, filepath: str):
        """
        Initializes the visualizer by loading data and extracting metadata.
        """
        self.filepath = Path(filepath)
        logger.info(f"Initializing DiaQcVisualizer for file: {self.filepath.name}")
        
        # 1. Load the raw data
        self.raw_data = self._load_data()
        self._validate_columns(self.raw_data)
        
        # 2. Extract metadata into a separate, editable dataframe
        self.metadata = self._extract_metadata(self.raw_data)
        
        # 3. Prepare a cache for the final merged data
        self._processed_data = None
        
        logger.info("Initialization complete. Raw data loaded and metadata extracted.")

    def get_metadata(self) -> pd.DataFrame:
        """
        Returns the extracted metadata DataFrame for viewing or editing.
        
        Returns:
            A pandas DataFrame with one row per unique run and extracted metadata columns.
        """
        return self.metadata

    def set_metadata(self, edited_metadata_df: pd.DataFrame):
        """
        Updates the internal metadata with a user-edited version.
        
        Args:
            edited_metadata_df: The DataFrame containing the user's modifications.
        """
        logger.info("Updating internal metadata with user-provided edits.")
        self.metadata = edited_metadata_df
        # Clear the cache for processed data, so it gets re-merged on the next call
        self._processed_data = None

    def get_processed_data(self) -> pd.DataFrame:
        """
        Returns the full dataset by merging raw data with the (potentially edited) metadata.
        
        This method caches the result to avoid re-merging on subsequent calls unless
        the metadata is updated.

        Returns:
            A pandas DataFrame containing the merged data, ready for plotting.
        """
        if self._processed_data is not None:
            logger.info("Returning cached processed data.")
            return self._processed_data

        logger.info("Merging raw data with metadata to create processed dataset.")
        # Merge the full raw data with the metadata table on the 'Run' column
        processed_df = pd.merge(self.raw_data, self.metadata, on='Run', how='left')
        
        # Cache the result
        self._processed_data = processed_df
        
        return self._processed_data

    def _load_data(self) -> pd.DataFrame:
        # (This method remains the same as before)
        if not self.filepath.is_file():
            raise FileNotFoundError(f"Error: The file '{self.filepath}' was not found.")
        if self.filepath.suffix.lower() != '.parquet':
            raise ValueError("Please use a .parquet file.")
        try:
            return pd.read_parquet(self.filepath)
        except Exception as e:
            raise IOError(f"Failed to read the file '{self.filepath}'. Reason: {e}")

    def _validate_columns(self, df: pd.DataFrame):
        # (This method remains the same as before)
        if 'Run' not in df.columns:
            raise KeyError("The data is missing the required 'Run' column.")

    def _extract_metadata(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Creates a summary DataFrame with one row per unique run and extracted metadata.
        """
        logger.info("Extracting unique run metadata for editing...")
        if 'Run' not in df.columns:
            return pd.DataFrame()
            
        # Create a new DataFrame with just the unique run names
        unique_runs_df = df[['Run']].drop_duplicates().reset_index(drop=True)
        runs = unique_runs_df['Run'].astype(str)

        # --- Generic patterns for extraction ---
        date_pattern = r'^(\d{8}|\d{6})'
        amount_pattern = r'(\d+ng)'
        well_pattern = r'_([A-Z]{2,3}\d{1,2})_'
        injection_id_pattern = r'_(\d+)$'
        lr_label_pattern = r'_(L|R)_'
        project_pattern = r'_([A-Za-z]+)_(?:[LR]_)?\d+ng'
        plate_info_pattern = r'([Pp]late[\s\w_]+?)_([A-Z]{2,3}\d{1,2})'

        # --- Perform extractions on the unique runs ---
        unique_runs_df['Project'] = runs.str.extract(project_pattern, expand=False).fillna('Unknown')
        unique_runs_df['Amount'] = runs.str.extract(amount_pattern, expand=False).fillna('Unknown')
        unique_runs_df['LR_Label'] = runs.str.extract(lr_label_pattern, expand=False).fillna('N/A')
        unique_runs_df['Plate_Info'] = runs.str.extract(plate_info_pattern, expand=False)[0].str.strip().fillna('Unknown')
        unique_runs_df['Well_Position'] = runs.str.extract(well_pattern, expand=False).fillna('Unknown')
        unique_runs_df['Injection_ID'] = pd.to_numeric(
            runs.str.extract(injection_id_pattern, expand=False), errors='coerce'
        )
        
        date_strs = runs.str.extract(date_pattern, expand=False)
        is_short_date = (date_strs.str.len() == 6) & (date_strs.notna())
        if is_short_date.any():
            date_strs.loc[is_short_date] = '20' + date_strs.loc[is_short_date]
        unique_runs_df['Acquisition_Date'] = pd.to_datetime(date_strs, format='%Y%m%d', errors='coerce')
        
        # --- Calculate Injection Order ---
        if 'Injection_ID' in unique_runs_df.columns and unique_runs_df['Injection_ID'].notna().any():
            unique_runs_df = unique_runs_df.sort_values('Injection_ID').reset_index(drop=True)
            unique_runs_df['Injection_Order'] = unique_runs_df.index + 1
            logger.info("Successfully calculated 'Injection_Order'.")
        else:
            logger.warning("Could not determine injection order. Assigning default order.")
            unique_runs_df['Injection_Order'] = unique_runs_df.index + 1

        return unique_runs_df
    
    def summarize_metadata(self, metadata_df: pd.DataFrame) -> Dict:
        """
        Calculates and returns a dictionary of summary statistics from the metadata.

        Args:
            metadata_df: The metadata DataFrame to be summarized.

        Returns:
            A dictionary containing key summary metrics and dataframes.
        """
        summary = {}
        if metadata_df.empty:
            return {"error": "Metadata is empty."}

        try:
            # --- General Stats ---
            summary['total_samples'] = metadata_df['Run'].nunique()
            summary['unique_dates'] = metadata_df['Acquisition_Date'].nunique()
            summary['start_date'] = metadata_df['Acquisition_Date'].min().date()
            summary['end_date'] = metadata_df['Acquisition_Date'].max().date()

            # --- Grouped Stats (as DataFrames) ---
            summary['samples_per_date'] = metadata_df.groupby('Acquisition_Date')['Run'].nunique().reset_index()
            summary['project_split'] = metadata_df.groupby('Project')['Run'].nunique().reset_index()
            summary['amount_split'] = metadata_df.groupby('Amount')['Run'].nunique().reset_index()
            
            # --- Special Run Designations ---
            # These columns are generated during the initial extraction, but might be edited.
            # Check if the boolean flag columns exist before summarizing.
            summary['check_runs'] = metadata_df['Is_Check_Sample'].sum() if 'Is_Check_Sample' in metadata_df.columns else 0
            summary['test_runs'] = metadata_df['Is_Test_Sample'].sum() if 'Is_Test_Sample' in metadata_df.columns else 0
            summary['bad_inj_runs'] = metadata_df['Is_Bad_Injection'].sum() if 'Is_Bad_Injection' in metadata_df.columns else 0
            summary['new_col_runs'] = metadata_df['New_Column_Run'].sum() if 'New_Column_Run' in metadata_df.columns else 0

        except Exception as e:
            logger.error(f"Error during metadata summarization: {e}", exc_info=True)
            return {"error": f"An error occurred during summarization: {e}"}

        return summary
    
    def find_stable_peptides(self, top_n: int = 5, min_detection_rate: float = 0.90) -> List[str]:
        """
        Finds high-quality, consistently detected peptides ideal for QC monitoring.

        This method identifies precursors that are detected in a high percentage of runs
        and then selects the most abundant ones from that stable set.

        Args:
            top_n (int): The number of sentinel peptides to return.
            min_detection_rate (float): The minimum fraction of runs a peptide must be
                                        present in to be considered stable (e.g., 0.9 means 90%).

        Returns:
            A list of the top N stable precursor IDs.

        Raises:
            KeyError: If required columns are not found in the data.
            ValueError: If no peptides meet the detection rate criteria.
        """
        logger.info(f"Finding top {top_n} sentinel peptides with >= {min_detection_rate:.0%} detection rate.")
        
        df = self.get_processed_data()

        # 1. Validate required columns
        required_cols = ['Run', 'Precursor.Id', 'Precursor.Quantity']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise KeyError(f"Missing required columns for finding stable peptides: {', '.join(missing_cols)}")

        # 2. Calculate total number of unique runs from the metadata
        num_runs = self.metadata['Run'].nunique()
        if num_runs == 0:
            return []

        # 3. Calculate detection rate and mean quantity for each precursor
        precursor_stats = df.groupby('Precursor.Id').agg(
            Detection_Count=('Run', 'nunique'),
            Mean_Quantity=('Precursor.Quantity', 'mean')
        )
        precursor_stats['Detection_Rate'] = precursor_stats['Detection_Count'] / num_runs

        # 4. Filter for consistently detected peptides
        consistent_peptides = precursor_stats[precursor_stats['Detection_Rate'] >= min_detection_rate]

        if consistent_peptides.empty:
            raise ValueError(f"No peptides found detected in at least {min_detection_rate:.0%} of runs. Please lower the detection rate threshold.")
            
        # 5. Select the top N most abundant peptides from the stable set
        sentinels = consistent_peptides.nlargest(top_n, 'Mean_Quantity').index.tolist()
        
        logger.info(f"Identified {len(sentinels)} sentinel peptides.")
        return sentinels
    
    def _filter_by_q_value(self, q_value_cutoff: float = 0.01) -> pd.DataFrame:
        """Filters the processed data for high-confidence precursors."""
        df = self.get_processed_data()
        
        required_cols = ['Decoy', 'Q.Value']
        if not all(col in df.columns for col in required_cols):
            raise KeyError("Data is missing 'Decoy' or 'Q.Value' columns for filtering.")
            
        filtered_df = df[
            (df['Decoy'] == 0) & 
            (df['Q.Value'] < q_value_cutoff)
        ].copy()
        
        logger.info(f"Filtered data to {len(filtered_df)} precursors with Q.Value < {q_value_cutoff}.")
        return filtered_df

    def plot_control_chart(self, peptide_id: str, metric_col: str, y_axis_title: str, q_value_cutoff: float):
        """Generates a Levey-Jennings style control chart for a specific peptide and metric."""
        df = self._filter_by_q_value(q_value_cutoff)
        peptide_df = df[df['Precursor.Id'] == peptide_id].sort_values('Injection_Order').copy()

        if peptide_df.empty:
            raise ValueError(f"No data found for peptide {peptide_id} at Q-value < {q_value_cutoff}.")

        mean = peptide_df[metric_col].mean()
        std = peptide_df[metric_col].std()
        
        upper_3_std, lower_3_std = mean + 3 * std, mean - 3 * std
        upper_2_std, lower_2_std = mean + 2 * std, mean - 2 * std
        upper_1_std, lower_1_std = mean + 1 * std, mean - 1 * std

        fig = go.Figure()
        fig.add_hrect(y0=lower_3_std, y1=upper_3_std, fillcolor='rgba(211,211,211,0.4)', line_width=0)
        fig.add_hrect(y0=lower_2_std, y1=upper_2_std, fillcolor='rgba(173,216,230,0.4)', line_width=0)
        fig.add_hrect(y0=lower_1_std, y1=upper_1_std, fillcolor='rgba(144,238,144,0.4)', line_width=0)
        fig.add_hline(y=mean, line_width=2, line_dash="dash", line_color="red")
        
        fig.add_trace(go.Scatter(
            x=peptide_df['Acquisition_Date'],
            y=peptide_df[metric_col],
            mode='lines+markers', name='Value',
            marker=dict(color='black', size=6),
            hoverinfo='text',
            text=[f"Run: {run}<br>Date: {date.date()}<br>Value: {val:.2f}" 
                  for run, date, val in zip(peptide_df['Run'], peptide_df['Acquisition_Date'], peptide_df[metric_col])]
        ))

        fig.update_layout(
            title=f'Control Chart for: {peptide_id}<br>Metric: {y_axis_title}',
            xaxis_title='Acquisition Date', yaxis_title=y_axis_title,
            title_x=0.5, showlegend=False
        )
        return fig

    def plot_rt_drift(self, sentinel_peptides: List[str], q_value_cutoff: float):
        """Plots the normalized RT drift for a list of sentinel peptides."""
        df = self._filter_by_q_value(q_value_cutoff)
        validation_df = df[df['Precursor.Id'].isin(sentinel_peptides)].copy()

        if validation_df.empty:
            raise ValueError("No data found for the selected sentinel peptides.")

        # --- THE FIX ---
        # Explicitly sort the DataFrame by injection order before plotting.
        # This ensures the lines are drawn chronologically.
        validation_df = validation_df.sort_values('Injection_Order')
            
        validation_df['RT_normalized'] = validation_df['RT'] - validation_df.groupby('Precursor.Id')['RT'].transform('mean')

        fig = px.line(
            validation_df, x='Injection_Order', y='RT_normalized', color='Precursor.Id',
            title='Systematic Retention Time Drift of Sentinel Peptides',
            labels={'Injection_Order': 'Injection Order', 'RT_normalized': 'RT Deviation from Mean (min)'},
            markers=True, hover_data={'RT': ':.2f'},
            color_discrete_sequence=px.colors.qualitative.Alphabet
        )
        fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="black")
        fig.update_layout(title_x=0.5)
        return fig

    def plot_peak_width_distribution(self, q_value_cutoff: float):
        """Creates a boxplot of chromatographic peak width (FWHM) by date."""
        df = self._filter_by_q_value(q_value_cutoff)
        fig = px.box(
            df.sort_values('Acquisition_Date'), x='Acquisition_Date', y='FWHM',
            title='Chromatographic Peak Width (FWHM) by Acquisition Date',
            labels={'Acquisition_Date': 'Acquisition Date', 'FWHM': 'Peak Width (min)'}
        )
        fig.update_layout(title_x=0.5)
        return fig
    
    def plot_peptide_elution_distribution(self, q_value_cutoff: float, num_bins: int):
        """
        Generates a histogram showing the distribution of unique peptide elutions over RT.

        Args:
            q_value_cutoff (float): The Q-value threshold for including precursors.
            num_bins (int): The number of bins to use for the histogram.

        Returns:
            A Plotly histogram figure object.
        """
        df = self._filter_by_q_value(q_value_cutoff)
        
        # Get a representative distribution using only unique precursors
        unique_precursors_df = df.drop_duplicates(subset=['Precursor.Id'])

        if unique_precursors_df.empty:
            raise ValueError("No unique precursors found after filtering.")

        fig = px.histogram(
            unique_precursors_df,
            x='RT',
            nbins=num_bins,
            title='Distribution of Peptide Elutions over Retention Time',
            labels={'RT': 'Retention Time (min)'}
        )

        fig.update_layout(
            title_x=0.5,
            bargap=0.1,
            yaxis_title="Number of Unique Peptides"
        )
        return fig
    

    def plot_rt_prediction_error(self, q_value_cutoff: float, moving_avg_window: int):
        """Plots the moving average of the RT prediction error with a symmetrical y-axis."""
        df = self._filter_by_q_value(q_value_cutoff)
        
        if 'Predicted.RT' not in df.columns:
            raise KeyError("Data is missing 'Predicted.RT' column required for this plot.")
            
        df['RT_Prediction_Error'] = df['RT'] - df['Predicted.RT']
        
        run_error_df = df.groupby('Injection_Order')['RT_Prediction_Error'].median().reset_index()
        run_error_df['Error_Rolling_Avg'] = run_error_df['RT_Prediction_Error'].rolling(
            window=moving_avg_window, center=True, min_periods=1
        ).mean()

        # --- THE FIX ---
        # 1. Find the maximum absolute error to set a symmetrical scale.
        max_abs_error = run_error_df['Error_Rolling_Avg'].abs().max()
        
        # 2. Add a 10% buffer so points aren't on the edge of the plot.
        y_range_limit = max_abs_error * 1.1 if max_abs_error > 0 else 1.0

        fig = px.line(
            run_error_df.sort_values('Injection_Order'), x='Injection_Order', y='Error_Rolling_Avg',
            title=f'RT Prediction Error ({moving_avg_window}-Run Moving Average)',
            labels={'Injection_Order': 'Injection Order', 'Error_Rolling_Avg': 'Prediction Error (min)'},
            markers=True
        )
        fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="black")
        
        # 3. Apply the new symmetrical y-axis range to the layout.
        fig.update_layout(
            title_x=0.5,
            yaxis_range=[-y_range_limit, y_range_limit]
        )
        return fig
    
    def plot_im_drift(self, sentinel_peptides: List[str], q_value_cutoff: float):
        """Plots the normalized Ion Mobility (IM) drift for a list of sentinel peptides."""
        df = self._filter_by_q_value(q_value_cutoff)
        
        if 'IM' not in df.columns:
            raise KeyError("Data is missing 'IM' column required for this plot.")
            
        validation_df = df[df['Precursor.Id'].isin(sentinel_peptides)].copy()

        if validation_df.empty:
            raise ValueError("No data found for the selected sentinel peptides.")
        
        # Sort chronologically to ensure lines are drawn correctly
        validation_df = validation_df.sort_values('Injection_Order')
            
        # Normalize the IM for each peptide to compare trends
        validation_df['IM_normalized'] = validation_df['IM'] - validation_df.groupby('Precursor.Id')['IM'].transform('mean')

        fig = px.line(
            validation_df, x='Injection_Order', y='IM_normalized', color='Precursor.Id',
            title='Systematic Ion Mobility Shift of Sentinel Peptides',
            labels={'Injection_Order': 'Injection Order', 'IM_normalized': 'IM Deviation from Mean (1/K0)'},
            markers=True,
            color_discrete_sequence=px.colors.qualitative.Alphabet
        )
        fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="black")
        fig.update_layout(title_x=0.5)
        return fig
    
    def plot_mass_accuracy_distribution(self, q_value_cutoff: float, y_range: list):
        """Creates a boxplot of mass accuracy (Mass.Evidence) by acquisition date."""
        df = self._filter_by_q_value(q_value_cutoff)
        
        if 'Mass.Evidence' not in df.columns:
            raise KeyError("Data is missing 'Mass.Evidence' column required for this plot.")

        fig = px.box(
            df.sort_values('Acquisition_Date'),
            x='Acquisition_Date',
            y='Mass.Evidence',
            title='Mass Accuracy Distribution by Acquisition Date',
            labels={'Acquisition_Date': 'Acquisition Date', 'Mass.Evidence': 'Mass Error (ppm)'}
        )
        fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="black")
        # Set y-axis range dynamically based on user input
        fig.update_yaxes(range=y_range)
        fig.update_layout(title_x=0.5)
        return fig

    def plot_sentinel_mass_accuracy(self, sentinel_peptides: List[str], q_value_cutoff: float):
        """Plots the mass accuracy for a list of sentinel peptides over time."""
        df = self._filter_by_q_value(q_value_cutoff)
        
        if 'Mass.Evidence' not in df.columns:
            raise KeyError("Data is missing 'Mass.Evidence' column required for this plot.")
            
        validation_df = df[df['Precursor.Id'].isin(sentinel_peptides)].copy()

        if validation_df.empty:
            raise ValueError("No data found for the selected sentinel peptides.")
            
        fig = px.line(
            validation_df.sort_values('Injection_Order'), 
            x='Injection_Order', y='Mass.Evidence', color='Precursor.Id',
            title='Mass Accuracy for Sentinel Peptides over Time',
            labels={'Injection_Order': 'Injection Order', 'Mass.Evidence': 'Mass Error (ppm)'},
            markers=True,
            color_discrete_sequence=px.colors.qualitative.Alphabet
        )
        fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="black")
        fig.update_layout(title_x=0.5)
        return fig

    def plot_mass_error_trend(self, q_value_cutoff: float, moving_avg_window: int):
        """Plots the smoothed, rolling average of the median mass error per run."""
        df = self._filter_by_q_value(q_value_cutoff)
        
        if 'Mass.Evidence' not in df.columns:
            raise KeyError("Data is missing 'Mass.Evidence' column required for this plot.")
            
        run_mass_error = df.groupby('Injection_Order')['Mass.Evidence'].median().reset_index()
        run_mass_error['Error_Rolling_Avg'] = run_mass_error['Mass.Evidence'].rolling(
            window=moving_avg_window, center=True, min_periods=1
        ).mean()

        fig = px.line(
            run_mass_error.sort_values('Injection_Order'),
            x='Injection_Order', y='Error_Rolling_Avg',
            title=f'Mass Error Trend ({moving_avg_window}-Run Rolling Average)',
            labels={'Injection_Order': 'Injection Order', 'Error_Rolling_Avg': 'Median Mass Error (ppm)'},
            markers=True
        )
        fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="black")
        fig.update_layout(title_x=0.5)
        return fig

    
    