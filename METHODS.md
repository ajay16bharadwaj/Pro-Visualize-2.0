# Methods — Pro-Visualize 2.0

Ready-to-adapt language for the Methods section of the Pro-Visualize paper.
It documents the analytical choices the application makes so they read as
deliberate design rather than oversights. Verified against the implementation
(branch `feature/p8-publication-hardening`); update if the code changes.

---

## Software and availability

Pro-Visualize 2.0 is a Streamlit application (Python 3.11). Analyses are
organized into modules for quality control (DIA and targeted), dilution-series
linearity, label-free quantification, comparative differential expression,
pathway enrichment, and single-cell proteomics (SCP). Each exported report
embeds a *Methods & Reproducibility* block recording the application version,
git commit, and the exact versions of all key dependencies
(`streamlit`, `pandas`, `numpy`, `scipy`, `scikit-learn`, `plotly`, `scanpy`,
`anndata`, `gseapy`, `gprofiler-official`, `umap-learn`), so every figure can
be traced to the environment that produced it. The application is distributed
as a pinned `requirements.txt` and a Docker image (non-root, Python 3.11.15).

## Statistical testing (Comparative module)

The Comparative module **does not perform differential-expression testing
in-app**. Fold-changes and adjusted p-values (FDR) are computed by the user's
upstream statistical workflow (e.g. MSstats, limma, or DESeq2) and imported as
pre-computed columns. The module visualizes these results (volcano, heatmap,
violin) and applies user-defined significance thresholds. Default volcano
thresholds are FDR ≤ 0.05 and |log₂ fold-change| ≥ 1, with guided presets
("most stringent" 0.01/1.5, "standard" 0.05/1.0, "exploratory" 0.10/0.5); all
are user-adjustable.

## Transformation, normalization, and missing values

Intensities of zero are treated as missing (`NaN`) prior to analysis. Log₂ is
used for fold-change, intensity distributions, and coefficient-of-variation
displays; log₁₀ is used for protein-abundance rank-order plots. Coefficient of
variation (CV %) is computed on the **linear** scale as 100 × (σ / μ).

For unsupervised analyses that require a complete matrix (PCA, hierarchical
clustering, sample-correlation), missing values are imputed by feature mean
(`sklearn.SimpleImputer`) and features standardized (`StandardScaler`) prior to
projection. Imputation is applied **only** for these exploratory visualizations
and never to reported statistics. Mean imputation can understate variance for
data missing not at random, which is common in bottom-up proteomics; results of
clustering/PCA should be interpreted with that caveat. Samples with no measured
values are dropped before projection, and analyses require ≥ 2 samples with
data.

## Dilution-series linearity (LOD/LOQ)

Per-protein response is fit by linear regression of log₂(intensity) on
log₂(concentration). The coefficient of determination (R²) quantifies linearity
in log–log space. Limits of detection and quantification are estimated from the
calibration fit as LOD = 2^(log₂C_min + 3.3·σ/slope) and
LOQ = 2^(log₂C_min + 10·σ/slope), where σ is the residual standard deviation of
the fit and slope is the regression slope. Proteins with a non-positive slope or
poor linearity (R² below threshold) are flagged and excluded from LOD/LOQ
estimation.

## Single-cell proteomics (SCP)

The SCP pipeline (scanpy/AnnData) applies log1p normalization, principal-
component analysis, neighbor-graph construction, UMAP embedding, and Leiden
community detection at a user-specified resolution. Optional batch correction
uses Harmony on the PCA coordinates. Differential expression between clusters
uses the Wilcoxon rank-sum test with Benjamini–Hochberg FDR correction
(`scanpy.tl.rank_genes_groups`), reporting per-group detection fractions.
Proteins below a configurable detection-rate threshold (default 10 % of cells)
are filtered prior to analysis. Pipeline parameters (PCs, neighbors, resolution,
batch key, DE test) are captured in the exported report.

## Quality control

DIA QC uses Levey–Jennings control charts with warning limits at mean ± Nσ
(N user-selectable, default 2) for ion-mobility, retention-time, mass-accuracy,
and peak-width metrics across runs, tracked on sentinel peptides. Targeted QC
reports per-peptide retention-time and peak-area CV (linear scale) and
peptide-level stability across injections.

## Pathway enrichment

Enrichment is computed against the GO (BP/CC/MF), KEGG, and Reactome libraries
via Enrichr. By default the **statistical background is the set of proteins
detected in the user's experiment** rather than the whole genome, queried
through Enrichr's Speedrichr background endpoints (and gseapy's custom
`background` for SCP). Using the detected proteome as the universe avoids
inflating terms composed of genes never observed in the experiment. A
whole-genome background (the classic Enrichr default) remains available as an
explicit option for comparison. Enrichr reports Benjamini–Hochberg adjusted
p-values per library.
