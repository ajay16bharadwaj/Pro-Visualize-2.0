# Pro-Visualize

A comprehensive, code-free proteomics data visualization and analysis platform for mass spectrometry workflows.

## Overview

Pro-Visualize addresses the fragmented nature of mass spectrometry data analysis by providing a unified, interactive environment that bridges instrument output to biological insights. Built on Streamlit, this application eliminates the need for researchers to juggle multiple disconnected tools by integrating quality control, quantitative analysis, comparative statistics, and pathway enrichment within a single platform.

### Key Benefits

- **Unified Workflow**: From raw instrument data to biological interpretation in one application
- **Code-Free Interface**: Interactive visualizations and analysis without programming knowledge
- **Multi-Platform Support**: Compatible with data from multiple mass spectrometry platforms
- **Reproducible Analysis**: Standardized workflows ensure consistency across experiments
- **Publication-Ready Outputs**: High-quality, presentation-grade visualizations with customizable aesthetics

## Features

Pro-Visualize consists of four core analytical modules:

### 1. Quality Control

Comprehensive quality assessment for both DIA (Data-Independent Acquisition) and targeted proteomics experiments.

**Key Capabilities:**
- Sample-level quality metrics and distributions
- Technical replicate reproducibility analysis
- Control chart monitoring for longitudinal QC
- Missing value pattern analysis
- Coefficient of variation (CV) assessment
- Detection completeness across samples

**Visualizations:**
- Sample distribution plots (box plots, violin plots)
- Correlation heatmaps and matrices
- PCA plots for sample clustering
- Control charts with statistical limits
- Missing data heatmaps
- CV distribution plots

### 2. Dilution Series Analysis

Specialized tools for analyzing loading curve experiments to assess linearity and quantification range.

**Key Capabilities:**
- Linearity assessment across concentration ranges
- R² calculation and visualization
- Detection limit determination
- Dynamic range evaluation
- Loading curve optimization

**Visualizations:**
- Loading curves with regression statistics
- R² value distributions
- Dilution series heatmaps
- Linearity deviation plots

### 3. Quantification Analysis

In-depth analysis of data completeness, distributions, and quantitative patterns.

**Key Capabilities:**
- Data completeness/wholeness assessment
- Protein abundance distribution analysis
- Missing value pattern characterization
- Zero value handling and replacement
- FDR (False Discovery Rate) validation
- Normalization strategy comparison

**Visualizations:**
- Stacked bar charts for data completeness
- Protein rank order plots
- Abundance distribution histograms
- Missing data patterns
- Pre/post normalization comparisons

### 4. Comparative Analysis

Statistical comparison between experimental groups with pathway enrichment integration.

**Key Capabilities:**
- Differential expression analysis
- Multiple statistical test options (t-test, ANOVA, Mann-Whitney)
- Multiple testing correction (Bonferroni, Benjamini-Hochberg)
- Fold change calculation and filtering
- Pathway enrichment analysis
- Gene ontology term enrichment
- Pathway database integration

**Visualizations:**
- Volcano plots with customizable thresholds
- MA plots
- Heatmaps of differentially expressed proteins
- Pathway enrichment bar charts
- GO term bubble plots
- Protein-pathway network diagrams

## Input Data Formats

Pro-Visualize accepts data from various proteomics analysis pipelines. Below are the supported formats and requirements for each module.

### Quality Control Module

**File Format:** CSV or TSV (tab-separated values)

**Required Columns:**
- `Protein.Ids` or `ProteinId`: Unique protein identifiers
- Sample columns: Numerical abundance/intensity values (one column per sample)

**Optional Columns:**
- `Protein.Names`: Protein names for annotation
- `Gene.Names`: Gene symbols
- `Fasta.Headers`: Full FASTA header information

**Example Structure:**
```
Protein.Ids,Sample_1,Sample_2,Sample_3,Sample_4
P12345,1500.5,1620.3,1480.2,1550.8
P23456,2300.1,2280.5,2310.9,2290.3
P34567,890.2,920.5,0,910.1
```

**Notes:**
- Missing values should be represented as 0, NaN, or empty cells
- Technical replicates should be clearly labeled with consistent naming
- Minimum 2 samples required; 3+ recommended for statistical analysis

### Dilution Series Module

**File Format:** CSV or TSV

**Required Columns:**
- `Protein.Ids` or `ProteinId`: Unique protein identifiers
- Dilution point columns: Named to indicate concentration/dilution level (e.g., `1X`, `2X`, `5X`, `10X`)

**Naming Convention:**
- Dilution columns should include numerical values indicating relative concentration
- Multiple replicates per dilution point supported

**Example Structure:**
```
Protein.Ids,1X_Rep1,1X_Rep2,2X_Rep1,2X_Rep2,5X_Rep1,5X_Rep2
P12345,1000.0,1020.5,2050.3,2100.1,5100.0,5200.5
P23456,1500.2,1480.9,3000.5,2980.3,7500.1,7600.8
```

**Notes:**
- At least 3 dilution points required for linearity assessment
- 5+ dilution points recommended for comprehensive analysis
- Technical replicates improve R² calculation reliability

### Quantification Module

**File Format:** CSV or TSV

**Required Columns:**
- `Protein.Ids` or `ProteinId`: Unique protein identifiers
- Sample columns: Numerical abundance values
- `Q.Value` or `FDR` (optional): False discovery rate for filtering

**Example Structure:**
```
Protein.Ids,Gene.Names,Sample_A1,Sample_A2,Sample_B1,Sample_B2,Q.Value
P12345,GENE1,1500.5,1620.3,1480.2,1550.8,0.001
P23456,GENE2,2300.1,2280.5,2310.9,2290.3,0.005
P34567,GENE3,890.2,NaN,0,910.1,0.015
```

**Notes:**
- Supports mixed missing value representations (0, NaN, empty)
- FDR/Q-value column used for statistical filtering if present
- Gene names enhance downstream interpretation

### Comparative Analysis Module

**File Format:** CSV or TSV

**Required Columns:**
- `Protein.Ids` or `ProteinId`: Unique protein identifiers
- `Gene.Names` or `Gene`: Gene symbols (required for pathway enrichment)
- Group columns: Numerical abundance values labeled by experimental group

**Optional Columns:**
- `Protein.Names`: Protein descriptions
- `Q.Value` or `FDR`: Existing FDR values
- `log2FC` or `FoldChange`: Pre-calculated fold changes

**Group Definition File (Optional):**
A separate CSV file mapping samples to experimental groups:
```
Sample,Group
Sample_A1,Control
Sample_A2,Control
Sample_A3,Control
Sample_B1,Treatment
Sample_B2,Treatment
Sample_B3,Treatment
```

**Example Structure:**
```
Protein.Ids,Gene.Names,Control_1,Control_2,Control_3,Treatment_1,Treatment_2,Treatment_3
P12345,GENE1,1500.5,1620.3,1550.8,2100.3,2200.5,2150.8
P23456,GENE2,2300.1,2280.5,2290.3,1800.2,1750.5,1820.1
```

**Notes:**
- Minimum 2 samples per group required; 3+ recommended
- Gene names are essential for pathway enrichment functionality
- Supports both raw and log-transformed data

## Data Preprocessing

Pro-Visualize includes comprehensive preprocessing pipelines:

### Missing Value Handling
- Automatic detection of zero vs. true missing values
- Multiple imputation strategies available
- Below detection limit (BDL) identification

### Normalization Options
- Median normalization
- Total intensity normalization
- Quantile normalization
- Z-score transformation

### Quality Filtering
- FDR/Q-value thresholding
- Minimum detection frequency filtering
- CV-based filtering
- Intensity threshold filtering

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager
- Git (for cloning repository)

### Option 1: Clone from GitHub

```bash
# Clone the repository
git clone https://github.com/yourusername/pro-visualize.git
cd pro-visualize

# Create a virtual environment (recommended)
python -m venv venv

# Activate virtual environment
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### Option 2: Install Dependencies Manually

```bash
# Core dependencies
pip install streamlit pandas numpy scipy
pip install plotly matplotlib seaborn
pip install scikit-learn
pip install openpyxl xlrd  # For Excel file support
```

### Requirements File

Create a `requirements.txt` file with the following dependencies:

```
streamlit>=1.28.0
pandas>=2.0.0
numpy>=1.24.0
scipy>=1.10.0
plotly>=5.14.0
matplotlib>=3.7.0
seaborn>=0.12.0
scikit-learn>=1.3.0
openpyxl>=3.1.0
xlrd>=2.0.1
statsmodels>=0.14.0
```

## Running the Application

### Local Development

```bash
# Ensure you're in the project directory and virtual environment is activated
streamlit run app.py
```

The application will open automatically in your default web browser at `http://localhost:8501`

### Configuration

Create a `.streamlit/config.toml` file for customization:

```toml
[theme]
primaryColor = "#1f77b4"
backgroundColor = "#ffffff"
secondaryBackgroundColor = "#f0f2f6"
textColor = "#262730"

[server]
port = 8501
enableCORS = false
enableXsrfProtection = true

[browser]
gatherUsageStats = false
```

## Usage Workflow

### Basic Analysis Pipeline

1. **Launch Application**: Run the Streamlit app
2. **Select Module**: Choose from Quality Control, Dilution Series, Quantification, or Comparative Analysis
3. **Upload Data**: Use the file uploader to load your proteomics data
4. **Configure Parameters**: Adjust analysis parameters and thresholds
5. **Generate Visualizations**: View interactive plots and analyses
6. **Export Results**: Download plots and statistical results

### Quality Control Workflow

```
Upload Data → Set FDR Threshold → Review Distributions → 
Check Correlations → Examine Missing Patterns → 
Assess CV Values → Export QC Report
```

### Comparative Analysis Workflow

```
Upload Data → Define Groups → Select Statistical Test → 
Set Thresholds (FC, p-value) → Generate Volcano Plot → 
Identify Significant Proteins → Run Pathway Enrichment → 
Export Results
```

## Project Structure

```
pro-visualize/
├── app.py                          # Main Streamlit application
├── requirements.txt                # Python dependencies
├── README.md                       # This file
├── .streamlit/
│   └── config.toml                # Streamlit configuration
├── visualizations/                # Backend visualizer classes
│   ├── qc_visualizer.py          # Quality control plots
│   ├── dilution_visualizer.py    # Dilution series analysis
│   ├── quant_visualizer.py       # Quantification plots
│   └── comparative_visualizer.py # Statistical comparisons
├── ui/                            # Frontend UI modules
│   ├── qc_ui.py
│   ├── dilution_ui.py
│   ├── quant_ui.py
│   └── comparative_ui.py
├── preprocessing/                 # Data preprocessing utilities
│   ├── data_loader.py
│   ├── normalization.py
│   └── missing_value_handler.py
└── utils/                         # Shared utilities
    ├── statistics.py
    ├── plotting_utils.py
    └── color_schemes.py
```

## Key Features in Detail

### Visualization Quality

All plots are designed for publication and presentation:
- **Colorblind-friendly palettes**: Default color schemes accessible to all users
- **High-resolution outputs**: Vector and raster formats for publication
- **Customizable aesthetics**: Adjust fonts, sizes, and styles
- **Interactive elements**: Hover information, zoom, pan capabilities
- **Statistical annotations**: R², p-values, confidence intervals displayed

### Statistical Rigor

- Multiple testing correction methods
- Assumption checking (normality, homogeneity of variance)
- Effect size calculations
- Confidence interval estimation
- Robust to missing data patterns

### Reproducibility Features

- Parameter logging for all analyses
- Session state preservation
- Export of analysis configurations
- Standardized data preprocessing pipelines

## Technical Notes

### Performance Considerations

- Optimized for datasets up to 10,000 proteins
- Recommended maximum file size: 100 MB
- Large dataset handling via sampling and chunking
- Caching enabled for improved responsiveness

### Browser Compatibility

- Chrome/Edge (recommended)
- Firefox
- Safari
- Minimum screen resolution: 1280x720

## Troubleshooting

### Common Issues

**Issue**: Upload fails with large files  
**Solution**: Increase Streamlit's max upload size in `.streamlit/config.toml`:
```toml
[server]
maxUploadSize = 200
```

**Issue**: Memory errors with large datasets  
**Solution**: Filter data before upload or use sampling options in the application

**Issue**: Plots not displaying  
**Solution**: Ensure JavaScript is enabled in browser; try refreshing the page

**Issue**: Missing dependencies  
**Solution**: Reinstall requirements: `pip install -r requirements.txt --upgrade`

## Citation

If you use Pro-Visualize in your research, please cite:

```
[Citation information to be added upon publication]
```

## Contributing

Contributions are welcome! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/new-feature`)
5. Create a Pull Request

### Development Guidelines

- Follow PEP 8 style guidelines
- Maintain the two-layer architecture (backend/frontend separation)
- Include docstrings for all functions and classes
- Add unit tests for new features
- Update documentation as needed

## License

[License information to be added]

## Contact

For questions, issues, or collaboration inquiries:

- GitHub Issues: [Project Issues Page]
- Email: [Contact email]

## Acknowledgments

- Built with [Streamlit](https://streamlit.io/)
- Visualization libraries: Plotly, Matplotlib, Seaborn
- Statistical analysis: SciPy, scikit-learn

## Roadmap

Planned features and improvements:

- [ ] Additional normalization methods
- [ ] Batch effect correction tools
- [ ] Enhanced pathway database integration
- [ ] Custom color palette builder
- [ ] Automated report generation
- [ ] Multi-omics integration support
- [ ] RESTful API for programmatic access

---

**Version**: 1.0.0  
**Last Updated**: November 2025  
**Status**: Active Development
