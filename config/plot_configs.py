"""
Centralized configuration for plot styling, validation thresholds, and
column-name fallbacks used across all Pro-Visualize modules.

Why a central config:
- Removes duplicated constants (TF list was defined identically in
  quant_visualizer.py and comparative_visualizer.py).
- Lets hardcoded thresholds (CV bands, deviation buckets, regex patterns)
  become user-configurable without code changes.
- Gives every module a single source of truth for theme names → Plotly
  templates and for column-detection fallbacks.

Importing convention: from config.plot_configs import THEMES, COLUMN_FALLBACKS, ...
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Themes — display label → Plotly template name. Used by every module's
# global-settings panel selectbox.
# ---------------------------------------------------------------------------

THEMES: dict[str, str] = {
    "Standard White": "plotly_white",
    "Dark Mode": "plotly_dark",
    "Minimal": "simple_white",
}

DEFAULT_THEME_LABEL: str = "Standard White"


# ---------------------------------------------------------------------------
# Color palettes — colorblind-safe by default. Okabe-Ito 8-color is the
# de-facto safe palette for scientific plots (works for all 3 common types
# of color-vision deficiency). Use this as the default when user hasn't
# set custom group colors.
# ---------------------------------------------------------------------------

OKABE_ITO_PALETTE: list[str] = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
]

COLORBLIND_PALETTE: list[str] = OKABE_ITO_PALETTE


# ---------------------------------------------------------------------------
# Human transcription factor gene set — used for highlighting in volcano
# and rank-order plots. Previously duplicated in quant_visualizer.py:29-37
# and comparative_visualizer.py:22-30; this is now the single source.
# ---------------------------------------------------------------------------

HUMAN_TRANSCRIPTION_FACTORS: frozenset[str] = frozenset({
    "ATF1", "ATF2", "ATF3", "ATF4", "ATF5", "ATF6", "ATF7",
    "CREB1", "CREB3", "CREB5",
    "FOS", "FOSB", "JUN", "JUNB", "JUND",
    "MYC", "MYCN", "MAX", "MAD",
    "CEBPA", "CEBPB", "CEBPD",
    "EGR1", "EGR2", "EGR3", "EGR4",
    "SP1", "SP2", "SP3", "SP4",
    "KLF1", "KLF2", "KLF4", "KLF5",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "NFKB1", "NFKB2", "RELA", "RELB",
    "TP53", "TP63", "TP73",
    "SOX2", "SOX9", "POU5F1", "NANOG",
    "GATA1", "GATA2", "GATA3", "GATA4", "GATA6",
    "RUNX1", "RUNX2", "RUNX3",
    "TCF3", "TCF4", "LEF1",
})


# ---------------------------------------------------------------------------
# Statistical-significance defaults. Surfaced as the initial values for
# Comparative-module FDR/FC sliders and SCP volcano cutoffs.
# ---------------------------------------------------------------------------

SIGNIFICANCE_DEFAULTS: dict[str, float] = {
    "fdr": 0.05,
    "log2fc": 1.0,
    "pval": 0.05,
}

SIGNIFICANCE_PRESETS: dict[str, dict[str, float]] = {
    "Most stringent": {"fdr": 0.01, "log2fc": 1.5},
    "Standard":       {"fdr": 0.05, "log2fc": 1.0},
    "Exploratory":    {"fdr": 0.10, "log2fc": 0.5},
}


# ---------------------------------------------------------------------------
# QC thresholds — defaults for Targeted/DIA QC control charts. Surfaced as
# user-configurable sliders in P5; kept here so visualizers can reference
# them without re-defining.
# ---------------------------------------------------------------------------

QC_THRESHOLDS: dict[str, float] = {
    "cv_low": 20.0,         # %CV — "good"
    "cv_high": 30.0,        # %CV — "warn above this"
    "control_sigma": 2.0,   # σ-bands on Levey-Jennings charts
}


# ---------------------------------------------------------------------------
# Dilution-series deviation buckets. Used to color expected-vs-observed
# fold-change comparisons. Currently hardcoded at dilution_series.py:114-119;
# P4 will surface these as sliders.
# ---------------------------------------------------------------------------

DEVIATION_BUCKETS: dict[str, float] = {
    "good": 0.2,   # |observed - expected| ≤ 0.2 → green
    "warn": 0.5,   # 0.2 < diff ≤ 0.5         → yellow; > 0.5 → red
}


# ---------------------------------------------------------------------------
# Column-name fallbacks. Maps a "logical role" (what the app needs) to
# common header strings used by various proteomics tools (DIA-NN, Spectronaut,
# MaxQuant, Skyline). Used by sanity.py to auto-detect columns and by the
# data-loading wizard to pre-fill text inputs.
# ---------------------------------------------------------------------------

COLUMN_FALLBACKS: dict[str, list[str]] = {
    "protein_id": [
        # Order matters — first match wins. Listed roughly by frequency
        # observed in the demo data sets (Pro-Vizualize-2.0-Demo).
        "Protein", "ProteinIds", "Protein.Group", "Protein.Ids",
        "Protein Group", "Protein IDs", "Accession", "Majority protein IDs",
    ],
    "gene": [
        "Gene Name", "Genes", "Gene", "Gene names",
        "Gene Symbol", "Gene.Symbol",
    ],
    "sample": [
        # SCP demo uses "Run"; Quant/Comp uses "Level3"; Targeted-QC uses
        # "Replicate Name"; DIA-NN long-format reports use "Run" (long) /
        # "File.Name" (stats). Filename comes last as it's path-prefixed.
        "Run", "Level3", "Replicate Name", "Sample",
        "File.Name", "Filename", "Raw file",
    ],
    "group": [
        # SCP demo: lowercase "condition". Quant demo:
        # "attribute_ExperimentalGroup". Generic: "Group" / "Condition".
        "condition", "attribute_ExperimentalGroup", "Group",
        "Condition", "Experimental Group", "biological_group",
        "Treatment",
    ],
    "concentration": [
        "Concentration", "concentration", "Conc", "Conc.",
        "Amount", "Dose", "Level", "Dilution", "ng/uL",
    ],
    "replicate": [
        "Replicate", "replicate", "Rep", "Replicate Number",
        "Bio Rep", "Bio Replicate",
    ],
    "fold_change": [
        "log2FC", "Log2FC", "Log2FoldChange", "logFC", "log2(FC)", "FC",
    ],
    "fdr": [
        "Imputed.FDR", "FDR", "adj.P.Val", "padj",
        "q.value", "qvalue", "Q.Value",
    ],
    "pvalue": [
        "p.value", "pvalue", "P.Value", "PValue", "p_val", "log10(p)",
    ],
}


# ---------------------------------------------------------------------------
# DIA-QC run-name regex patterns. Currently hardcoded in DiaQcVisualizer
# (lines 109-115, 128-131). P5 will let users add their own patterns.
# Patterns are tried in order; the first match wins.
# ---------------------------------------------------------------------------

DIA_RUN_NAME_PATTERNS: dict[str, list[str]] = {
    "date": [
        r"(\d{8})",   # YYYYMMDD
        r"(\d{6})",   # YYMMDD
    ],
    "amount_ng": [
        r"(\d+)\s*ng",
        r"(\d+)ng",
    ],
    "well": [
        r"([A-H]\d{1,2})",
    ],
    "project": [
        r"([A-Z]+\d+)",
    ],
}


# ---------------------------------------------------------------------------
# Plot-export defaults — used by PlotManager static-export buttons.
# ---------------------------------------------------------------------------

EXPORT_DEFAULTS: dict[str, int | float | str] = {
    "png_scale": 2,           # 2× DPI for retina-quality PNG
    "html_plotlyjs": "cdn",   # keeps HTML files small; offline mode = "include"
    "matplotlib_dpi": 150,    # for BytesIO PNGs from matplotlib visualizers
}
