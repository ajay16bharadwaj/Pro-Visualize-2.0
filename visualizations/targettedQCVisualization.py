# pro_visualize/visualizations/targettedQCVisualization.py

import pandas as pd
import re
import logging
from pathlib import Path
from typing import Dict, Optional, Any
import plotly.graph_objects as go
import plotly.express as px

# Set up a logger for this module
# The configuration (level, handler) should be set in the main app
logger = logging.getLogger(__name__)

# pro_visualize/qc/targeted_qc.py

import pandas as pd
import re
import logging
from pathlib import Path
from typing import Dict, Optional, Any
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# Set up a logger for this module
# The configuration (level, handler) should be set in the main app
logger = logging.getLogger(__name__)

class TargetedQcVisualizer:
    """
    Handles QC data loading, processing, and visualization for targeted proteomics.

    This class is initialized with a filepath to a report file (CSV or TSV).
    It provides methods to extract metadata and generate various QC plots
    using Plotly, with customizable configurations.

    Attributes:
        filepath (Path): The path to the input data file.
        data (pd.DataFrame): The loaded and validated proteomics data.
    """
    
    def __init__(self, filepath: str):
        """
        Initializes the visualizer by loading and validating the data file.
        """
        self.filepath = Path(filepath)
        logger.info(f"Initializing TargetedQcVisualizer for file: {self.filepath.name}")
        self.data = self._load_data()
        self._validate_columns()
        
        # Add this line to initialize the metadata cache
        self.metadata = None
        
        logger.info("Initialization complete. Data loaded and validated successfully.")

    def _load_data(self) -> pd.DataFrame:
        """
        Loads data from the specified file, auto-detecting the separator.
        """
        if not self.filepath.is_file():
            logger.error(f"File not found at path: {self.filepath}")
            raise FileNotFoundError(f"Error: The file '{self.filepath}' was not found.")

        file_suffix = self.filepath.suffix.lower()
        
        if file_suffix == '.csv':
            separator = ','
        elif file_suffix == '.tsv':
            separator = '\t'
        else:
            logger.error(f"Unsupported file type: {file_suffix}")
            raise ValueError(f"Unsupported file type '{file_suffix}'. Please use .csv or .tsv.")
        
        logger.info(f"Loading '{self.filepath.name}' with separator: '{separator}'")
        try:
            df = pd.read_csv(self.filepath, sep=separator)
            logger.info(f"Successfully loaded {len(df)} rows from the file.")
            return df
        except Exception as e:
            logger.error(f"Failed to read file '{self.filepath.name}'. Reason: {e}", exc_info=True)
            raise IOError(f"Failed to read the file '{self.filepath}'. Reason: {e}")

    def _validate_columns(self):
        """
        Checks if required columns for metadata extraction are present.
        """
        required = ['Replicate Name', 'Acquired Time']
        missing = [col for col in required if col not in self.data.columns]
        if missing:
            logger.error(f"Missing required columns in the input file: {missing}")
            raise KeyError(f"The data is missing required columns: {', '.join(missing)}")
        logger.debug("Required columns ['Replicate Name', 'Acquired Time'] found.")

    def extract_metadata(self) -> pd.DataFrame:
        """
        Parses 'Replicate Name' and 'Acquired Time' to create a metadata table.

        This method checks for a cached version of the metadata first. If not
        found, it generates, caches, and then returns the metadata.
        """
        # Check if metadata has already been generated and cached
        if self.metadata is not None:
            logger.info("Returning cached metadata.")
            return self.metadata

        logger.info("No cached metadata found. Starting extraction process.")
        try:
            # 1. Get unique replicates and their acquisition time
            meta_df = self.data[['Replicate Name', 'Acquired Time']].drop_duplicates().reset_index(drop=True)
            logger.debug(f"Found {len(meta_df)} unique replicates.")
            
            # (The rest of the calculation logic remains exactly the same...)
            # 2. Convert 'Acquired Time'...
            logger.info("Converting 'Acquired Time' to datetime objects.")
            try:
                meta_df['Acquired Time'] = pd.to_datetime(
                    meta_df['Acquired Time'], format='%m/%d/%Y %I:%M:%S %p', errors='raise'
                )
                logger.info("Parsed 'Acquired Time' using the format '%%m/%%d/%%Y %%I:%%M:%%S %%p'.")
            except ValueError:
                logger.warning(
                    "Could not parse 'Acquired Time' with the standard format. "
                    "Falling back to automatic inference."
                )
                meta_df['Acquired Time'] = pd.to_datetime(meta_df['Acquired Time'], errors='coerce')

            # 3. Define regex pattern...
            pattern = re.compile(
                r'^(?P<Date>\d{8})_'
                r'(?P<Operator>[A-Z]+)_'
                r'(?P<SampleType>.+?)_'
                r'(?P<Amount>\d+ng)_'
                r'(?P<Details>.+?)_'
                r'Batch(?P<Batch>\d+)_'
                r'(?P<Injection>\d+)$'
            )
            # 4. Apply regex...
            logger.info("Parsing 'Replicate Name' to extract metadata attributes.")
            extracted_data = meta_df['Replicate Name'].str.extract(pattern)
            # 5. Combine...
            meta_df = pd.concat([meta_df, extracted_data], axis=1)
            # 6. Sort...
            meta_df = meta_df.sort_values(by='Acquired Time').reset_index(drop=True)
            # 7. Add new columns...
            logger.info("Adding 'Acquisition Date' and 'Run Order' columns.")
            meta_df['Acquisition Date'] = meta_df['Acquired Time'].dt.date
            meta_df['Run Order'] = meta_df.index + 1
            meta_df.loc[meta_df['Acquired Time'].isna(), 'Run Order'] = pd.NA
            logger.info("Creating short replicate names for plot axes.")
            # Default to the full name
            meta_df['Short Replicate Name'] = meta_df['Replicate Name'] 
            # Find rows where Batch and Injection were successfully parsed
            parsed_rows = meta_df['Batch'].notna() & meta_df['Injection'].notna()
            # Create the short name for just those rows
            meta_df.loc[parsed_rows, 'Short Replicate Name'] = \
                'Batch' + meta_df.loc[parsed_rows, 'Batch'].astype(str) + \
                '_' + meta_df.loc[parsed_rows, 'Injection'].astype(str)

            # Cache the result before returning
            logger.info("Metadata extraction complete. Caching result for future use.")
            self.metadata = meta_df
            
            return self.metadata

        except Exception as e:
            logger.error(f"An error occurred during metadata extraction: {e}", exc_info=True)
            raise RuntimeError(f"An error occurred during metadata extraction: {e}")

    def _create_default_plot_config(self) -> Dict[str, Any]:
        """Returns a default dictionary for Plotly figure styling."""
        return {
            "template": "plotly_white",
            "font_family": "Arial, sans-serif",
            "font_size": 12,
            "title_font_size": 18,
            "color_discrete_sequence": px.colors.qualitative.Plotly,
            "margin": {"l": 60, "r": 20, "t": 50, "b": 60}
        }

    def _get_top_n_peptides_data(self, top_n: int) -> pd.DataFrame:
        """
        Selects data for the most abundant peptides for RT visualization.
        
        Args:
            top_n (int): The number of top abundant peptides to select.

        Returns:
            pd.DataFrame: A filtered DataFrame containing only the data for the
                          top N peptides.
        
        Raises:
            KeyError: If required columns for ranking are missing.
        """
        rank_cols = ['Peptide', 'Total Area']
        missing = [col for col in rank_cols if col not in self.data.columns]
        if missing:
            msg = f"Cannot rank peptides. Missing required columns: {missing}"
            logger.error(msg)
            raise KeyError(msg)
        
        logger.info(f"Selecting top {top_n} most abundant peptides based on 'Total Area'.")
        # Find the most abundant peptides across all runs
        top_peptides = (
            self.data.groupby('Peptide')['Total Area']
            .sum()
            .nlargest(top_n)
            .index
        )
        
        # Filter the main dataframe for these peptides
        return self.data[self.data['Peptide'].isin(top_peptides)]


    def plot_retention_time_stability(
        self, 
        top_n_peptides: int = 30, 
        plot_config: Optional[Dict[str, Any]] = None
    ) -> go.Figure:
        """
        Generates a scatter plot to visualize retention time stability.

        This plot shows the retention time for the N most abundant peptides
        across all replicates, ordered chronologically. It helps assess the
        stability of the chromatography.

        Args:
            top_n_peptides (int): The number of most abundant peptides to display.
                                  Defaults to 30.
            plot_config (Optional[Dict[str, Any]]): A dictionary to override
                default plot styling (e.g., colors, fonts).

        Returns:
            go.Figure: A Plotly scatter plot figure object.
            
        Raises:
            KeyError: If required columns for plotting are missing.
        """
        # 1. Validate required columns for this specific plot
        plot_cols = ['Best Retention Time', 'Replicate Name', 'Peptide']
        missing = [col for col in plot_cols if col not in self.data.columns]
        if missing:
            msg = f"Cannot create plot. Missing required columns: {missing}"
            logger.error(msg)
            raise KeyError(msg)
            
        # 2. Get default config and update with user's config
        config = self._create_default_plot_config()
        if plot_config:
            config.update(plot_config)
            
        # 3. Prepare the data for plotting
        try:
            plot_data = self._get_top_n_peptides_data(top_n=top_n_peptides)
            metadata = self.extract_metadata()
        except (KeyError, RuntimeError) as e:
            logger.error(f"Failed to prepare data for RT stability plot: {e}")
            # Return an empty figure with an error message
            fig = go.Figure()
            fig.update_layout(title="Error Preparing Data", annotations=[{
                "text": "Could not prepare data for plot.<br>Check logs for details.",
                "showarrow": False
            }])
            return fig
            
        # 4. Merge with metadata to get Run Order and the new Short Replicate Name
        plot_data = pd.merge(
            plot_data,
            metadata[['Replicate Name', 'Run Order', 'Short Replicate Name']],
            on='Replicate Name',
            how='left'
        )
        plot_data = plot_data.sort_values(by=['Run Order', 'Peptide'])
        
        # Define the chronological order for the x-axis using the short names
        replicate_order = plot_data.drop_duplicates(
            subset=['Replicate Name']
        ).sort_values('Run Order')['Short Replicate Name']

        logger.info("Generating retention time stability plot.")
        
        # 5. Create the figure, using the short name for the x-axis
        fig = px.scatter(
            plot_data,
            x='Short Replicate Name',           # Use the new short name here
            y='Best Retention Time',
            color='Peptide',
            hover_name='Replicate Name',        # Show the full name on hover
            title='Retention Time Stability Across Runs',
            labels={
                'Short Replicate Name': 'Replicate', # Update the axis label
                'Best Retention Time': 'Retention Time (minutes)',
                'Peptide': 'Peptide'
            },
            template=config['template'],
            category_orders={'Short Replicate Name': replicate_order} # Use short name here too
        )

        # 6. Apply final layout customizations
        fig.update_layout(
            font=dict(family=config['font_family'], size=config['font_size']),
            title_font_size=config['title_font_size'],
            margin=config['margin'],
            xaxis_tickangle=-45
        )
        fig.update_traces(
            marker=dict(size=8, symbol='line-ns-open', line=dict(width=2)),
            selector=dict(mode='markers')
        )
        
        logger.info("Successfully generated retention time stability plot.")
        return fig

    def plot_peak_area(
        self,
        top_n_peptides: int = 15,
        plot_config: Optional[Dict[str, Any]] = None
    ) -> go.Figure:
        """
        Generates a line plot to visualize peak area stability.

        This plot shows the peak area for the N most abundant peptides
        across all replicates, ordered chronologically. It helps assess the
        stability of the MS signal.

        Args:
            top_n_peptides (int): The number of most abundant peptides to display.
                                  Defaults to 15.
            plot_config (Optional[Dict[str, Any]]): A dictionary to override
                default plot styling.

        Returns:
            go.Figure: A Plotly line plot figure object.
        """
        # 1. Validate required columns for this specific plot
        plot_cols = ['Total Area', 'Replicate Name', 'Peptide']
        missing = [col for col in plot_cols if col not in self.data.columns]
        if missing:
            msg = f"Cannot create plot. Missing required columns: {missing}"
            logger.error(msg)
            raise KeyError(msg)
        
        # 2. Get default config and update with user's config
        config = self._create_default_plot_config()
        if plot_config:
            config.update(plot_config)
            
        # 3. Prepare the data for plotting
        try:
            plot_data = self._get_top_n_peptides_data(top_n=top_n_peptides).copy()
            metadata = self.extract_metadata()
        except (KeyError, RuntimeError) as e:
            logger.error(f"Failed to prepare data for peak area plot: {e}")
            fig = go.Figure()
            fig.update_layout(title="Error Preparing Data", annotations=[{
                "text": "Could not prepare data for plot.<br>Check logs for details.",
                "showarrow": False
            }])
            return fig
            
        # 4. Scale the peak area and merge with metadata for sorting
        plot_data['Peak Area (10^6)'] = plot_data['Total Area'] / 1_000_000
        plot_data = pd.merge(
            plot_data,
            metadata[['Replicate Name', 'Run Order', 'Short Replicate Name']],
            on='Replicate Name',
            how='left'
        ).sort_values(by='Run Order')
        
        replicate_order = plot_data.drop_duplicates(
            'Replicate Name'
        ).sort_values('Run Order')['Short Replicate Name']

        logger.info("Generating peak area stability plot.")
        
        # 5. Create the figure using px.line
        fig = px.line(
            plot_data,
            x='Short Replicate Name',
            y='Peak Area (10^6)',
            color='Peptide',
            markers=True,  # Add markers to the lines
            hover_name='Replicate Name',
            title='Peak Area Across Runs',
            labels={
                'Short Replicate Name': 'Replicate',
                'Peak Area (10^6)': 'Peak Area (10^6)',
                'Peptide': 'Peptide'
            },
            template=config['template'],
            category_orders={'Short Replicate Name': replicate_order}
        )

        # 6. Apply final layout customizations
        fig.update_layout(
            font=dict(family=config['font_family'], size=config['font_size']),
            title_font_size=config['title_font_size'],
            margin=config['margin'],
            xaxis_tickangle=-45
        )
        
        logger.info("Successfully generated peak area stability plot.")
        return fig

    def calculate_peptide_stats(self) -> pd.DataFrame:
        """
        Calculates summary statistics for each peptide across all runs.

        This method computes the mean retention time (RT), the retention time
        Coefficient of Variation (RT_CV), and the peak area Coefficient of
        Variation (Area_CV) for each unique peptide.

        Returns:
            pd.DataFrame: A table of peptide statistics, sorted by RT.
            
        Raises:
            KeyError: If required data columns are missing.
        """
        # 1. Validate required columns for this calculation
        stat_cols = ['Peptide', 'Best Retention Time', 'Total Area']
        missing = [col for col in stat_cols if col not in self.data.columns]
        if missing:
            msg = f"Cannot calculate stats. Missing required columns: {missing}"
            logger.error(msg)
            raise KeyError(msg)
            
        logger.info("Calculating peptide summary statistics (RT, RT_CV, Area_CV).")
        
        # 2. Group by peptide and calculate mean and standard deviation
        stats_df = self.data.groupby('Peptide').agg(
            RT_mean=('Best Retention Time', 'mean'),
            RT_std=('Best Retention Time', 'std'),
            Area_mean=('Total Area', 'mean'),
            Area_std=('Total Area', 'std')
        )
        
        # 3. Calculate Coefficient of Variation (CV) as (std/mean) * 100
        # Handle division by zero by replacing resulting NaN/inf with 0
        stats_df['RT_CV'] = (stats_df['RT_std'] / stats_df['RT_mean']) * 100
        stats_df['Area_CV'] = (stats_df['Area_std'] / stats_df['Area_mean']) * 100
        stats_df.fillna(0, inplace=True) # Replace NaN values with 0

        # 4. Format the final table to match the desired output
        final_table = stats_df[['RT_mean', 'RT_CV', 'Area_CV']].reset_index()
        final_table = final_table.rename(columns={
            'Peptide': 'Peptide Sequence',
            'RT_mean': 'RT'
        })
        
        # 5. Sort by retention time and round for a clean display
        final_table = final_table.sort_values(by='RT').reset_index(drop=True)
        final_table = final_table.round({'RT': 1, 'RT_CV': 1, 'Area_CV': 1})
        
        logger.info("Successfully calculated peptide summary statistics.")
        return final_table

    def plot_rt_distribution(
        self,
        plot_config: Optional[Dict[str, Any]] = None
    ) -> go.Figure:
        """
        Generates an RT distribution "lollipop" plot of the most intense
        instance of each peptide.

        This plot visualizes the overall chromatographic separation and peptide
        intensity distribution.

        Args:
            plot_config (Optional[Dict, Any]]): A dictionary to override
                default plot styling.

        Returns:
            go.Figure: A Plotly figure object.
        """
        # 1. Validate required columns
        plot_cols = ['Peptide', 'Best Retention Time', 'Total Area', 'Protein']
        missing = [col for col in plot_cols if col not in self.data.columns]
        if missing:
            msg = f"Cannot create plot. Missing required columns: {missing}"
            logger.error(msg)
            raise KeyError(msg)

        # 2. Get default config and update with user's config
        config = self._create_default_plot_config()
        if plot_config:
            config.update(plot_config)

        # 3. Find the most intense measurement for each peptide
        logger.info("Finding the most intense measurement for each peptide.")
        best_peptides = self.data.loc[self.data.groupby('Peptide')['Total Area'].idxmax()].copy()
        
        # 4. Prepare data for plotting
        best_peptides['Intensity (10^6)'] = best_peptides['Total Area'] / 1_000_000
        best_peptides['Short Label'] = best_peptides['Peptide'].str[:3] # 3-letter label
        
        logger.info("Generating RT distribution plot.")
        
        # 5. Create the figure using graph_objects for a lollipop chart
        fig = go.Figure()
        # Create the "stems" using thin bars
        fig.add_trace(go.Bar(
            x=best_peptides['Best Retention Time'],
            y=best_peptides['Intensity (10^6)'],
            width=0.01,  # Make bars very thin to act as stems
            marker_color='#d3d3d3', # Light grey color for stems
            hoverinfo='none',      # Stems don't need hover text
            showlegend=False
        ))

        # Create the "lollipops" (markers) and text labels
        fig.add_trace(go.Scatter(
            x=best_peptides['Best Retention Time'],
            y=best_peptides['Intensity (10^6)'],
            mode='markers+text',
            text=best_peptides['Short Label'],
            textposition='top center',
            marker=dict(
                size=8,
                color=best_peptides['Intensity (10^6)'], # Color markers by intensity
                colorscale='viridis',
                showscale=False
            ),
            # Define rich hover text
            customdata=best_peptides[['Peptide', 'Protein', 'Total Area']],
            hovertemplate=(
                '<b>%{customdata[0]}</b><br>'
                'Protein: %{customdata[1]}<br>'
                'Retention Time: %{x:.2f} min<br>'
                'Peak Area: %{customdata[2]:.2e}'
                '<extra></extra>' # Hides the trace name
            ),
            showlegend=False
        ))

        # 6. Apply final layout customizations
        fig.update_layout(
            title_text='RT Distribution',
            xaxis_title='Retention Time (minutes)',
            yaxis_title='Intensity (10^6)',
            font=dict(family=config['font_family'], size=config['font_size']),
            title_font_size=config['title_font_size'],
            margin=config['margin'],
            template=config['template'],
            bargap=0 # Remove gap between bars
        )
        
        logger.info("Successfully generated RT distribution plot.")
        return fig

    def plot_cv_distributions(
        self,
        area_cv_threshold: float = 20.0,
        rt_cv_threshold: float = 2.0,
        plot_config: Optional[Dict[str, Any]] = None
    ) -> go.Figure:
        """
        Generates histograms of the Area and RT Coefficient of Variation (CV).

        This plot provides an experiment-wide overview of reproducibility for
        both peak area and retention time.

        Args:
            area_cv_threshold (float): A QC threshold line for the Area CV plot.
                                       Defaults to 20.0 (i.e., 20%).
            rt_cv_threshold (float): A QC threshold line for the RT CV plot.
                                     Defaults to 2.0 (i.e., 2%).
            plot_config (Optional[Dict, Any]]): A dictionary to override
                default plot styling.

        Returns:
            go.Figure: A Plotly figure object with two subplots.
        """
        # 1. Get default config and update with user's config
        config = self._create_default_plot_config()
        if plot_config:
            config.update(plot_config)
            
        # 2. Get the peptide statistics
        try:
            stats_df = self.calculate_peptide_stats()
        except (KeyError, RuntimeError) as e:
            logger.error(f"Failed to get stats for CV plot: {e}")
            fig = go.Figure()
            fig.update_layout(title="Error Preparing Data")
            return fig
            
        logger.info("Generating CV distribution histograms.")

        # 3. Create a figure with two subplots
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=("Area CV Distribution", "RT CV Distribution")
        )

        # 4. Add the Area CV histogram to the first subplot
        fig.add_trace(go.Histogram(
            x=stats_df['Area_CV'],
            name='Area CV',
            marker_color=config['color_discrete_sequence'][0],
            nbinsx=30
        ), row=1, col=1)

        # 5. Add the RT CV histogram to the second subplot
        fig.add_trace(go.Histogram(
            x=stats_df['RT_CV'],
            name='RT CV',
            marker_color=config['color_discrete_sequence'][1],
            nbinsx=30
        ), row=1, col=2)

        # 6. Add vertical lines for QC thresholds
        fig.add_vline(
            x=area_cv_threshold, line_width=2, line_dash="dash", line_color="red",
            annotation_text=f"Threshold: {area_cv_threshold}%",
            annotation_position="top right",
            row=1, col=1
        )
        fig.add_vline(
            x=rt_cv_threshold, line_width=2, line_dash="dash", line_color="red",
            annotation_text=f"Threshold: {rt_cv_threshold}%",
            annotation_position="top right",
            row=1, col=2
        )

        # 7. Update layout
        fig.update_layout(
            title_text='Peptide Measurement Reproducibility',
            template=config['template'],
            font=dict(family=config['font_family'], size=config['font_size']),
            showlegend=False
        )
        fig.update_xaxes(title_text="CV (%)", row=1, col=1)
        fig.update_xaxes(title_text="CV (%)", row=1, col=2)
        fig.update_yaxes(title_text="Number of Peptides", row=1, col=1)
        
        logger.info("Successfully generated CV distribution plot.")
        return fig

    def get_failing_peptides(
        self,
        area_cv_threshold: float = 20.0,
        rt_cv_threshold: float = 2.0
    ) -> pd.DataFrame:
        """
        Returns a DataFrame of peptides that fail the specified CV thresholds.

        This method provides a detailed list of peptides that do not meet the
        reproducibility criteria, complementing the CV distribution plots.

        Args:
            area_cv_threshold (float): The QC threshold for Area CV (%).
            rt_cv_threshold (float): The QC threshold for RT CV (%).

        Returns:
            pd.DataFrame: A table of peptides exceeding one or both thresholds.
        """
        logger.info(f"Filtering for peptides with Area CV > {area_cv_threshold}% or RT CV > {rt_cv_threshold}%.")
        
        # 1. Get the full statistics table
        stats_df = self.calculate_peptide_stats()

        # 2. Identify peptides failing each threshold
        fails_area = stats_df['Area_CV'] > area_cv_threshold
        fails_rt = stats_df['RT_CV'] > rt_cv_threshold

        # 3. Filter for peptides that fail at least one of the conditions
        failing_df = stats_df[fails_area | fails_rt].copy()

        # 4. Add a 'Reason' column for clarity
        reasons = []
        for index, row in failing_df.iterrows():
            reason_parts = []
            if row['Area_CV'] > area_cv_threshold:
                reason_parts.append("High Area CV")
            if row['RT_CV'] > rt_cv_threshold:
                reason_parts.append("High RT CV")
            reasons.append(" & ".join(reason_parts))
        
        failing_df['Reason'] = reasons
        
        if failing_df.empty:
            logger.info("No peptides failed the specified QC thresholds. Great!")
        else:
            logger.warning(f"Found {len(failing_df)} peptides failing QC thresholds.")
            
        return failing_df
    