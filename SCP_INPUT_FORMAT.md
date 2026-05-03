# SCP Input Format Reference — Pro-Visualize

## Overview

The Single-Cell Proteomics (SCP) module expects three input files, matching the standard
DIA-NN output structure used in your timsTOF Pro workflows.

---

## File 1: PG Matrix (required)

**Source:** DIA-NN → `*_pg_matrix.tsv` (or any wide-format protein × sample matrix)

**Format:** Tab- or comma-separated, wide format

| Protein.Group | Run1 | Run2 | Run3 | … |
|---------------|------|------|------|---|
| GENE_A        | 1234 | 0    | 5678 | … |
| GENE_B        | 0    | 890  | 1100 | … |

- **First column** = protein group identifier (set "Protein column" in the UI).
  DIA-NN default is `Protein.Group`; you can also use `Genes` if your matrix was exported
  with gene-level rollup.
- **Remaining columns** = one column per sample/run (intensity values, e.g., `PG.MaxLFQ`).
  Missing values should be `NaN` or `0` — the class treats both as "not detected".
- Column headers must match (or be basenames of) the `Run` column in the Annotation file.

---

## File 2: Annotation (required)

**Format:** CSV or TSV

| Run          | cell_line       | condition | concentration | batch |
|--------------|-----------------|-----------|---------------|-------|
| Sample_P3_01 | Male_02i_CTR    | Control   | 0uM           | P3    |
| Sample_P3_02 | Female_083i_CTR | Treated   | 3uM           | P3    |
| Sample_P4_01 | Female_Line2    | Treated   | 10uM          | P4    |

- **Sample column** (set in UI, default `Run`) must match the column headers in the PG matrix.
  Full file paths in the PG matrix are matched by basename (without extension), so
  `/path/to/Sample_P3_01.d` matches annotation value `Sample_P3_01`.
- All remaining columns become grouping/metadata variables available in the UI
  for coloring, QC grouping, DE groupby, and activity-score violins.
- Columns with ≤30 unique values or string dtype appear in the "Group by" selectors.

---

## File 3: report.stats (optional but recommended)

**Source:** DIA-NN → `report.stats.tsv`

**Format:** Tab-separated (standard DIA-NN output)

| File.Name        | Proteins.Identified | FWHM.RT | Normalisation.Instability | MS1.Signal |
|------------------|---------------------|---------|---------------------------|------------|
| /path/Sample.raw | 1452                | 0.12    | 0.034                     | 8.3e7      |

- The `File.Name` column (set via "Run column" in UI, default `File.Name`) is matched
  to PG matrix columns by basename.
- The following columns are automatically pulled in if present:
  - `Proteins.Identified` — detection count (corroborates `n_proteins`)
  - `FWHM.RT` — chromatographic peak width (technical covariate candidate)
  - `Normalisation.Instability` — DIA-NN normalization stability (strong technical confound)
  - `MS1.Signal` — overall MS1 signal level
- These columns appear in QC violin plots and the covariate-correlation heatmap,
  and can be selected as regression covariates in the Preprocessing tab.

---

## Recommended Regression Covariates

From your Bruker AIP analyses, the most reliable technical confounders are:

1. `log_n_proteins` — log1p-transformed protein count per cell (auto-computed).
   Highly correlated with PC1 (~0.95). **Always include.**
2. `Normalisation.Instability` — DIA-NN normalization quality score (~-0.83 with PC1).
   **Include if present.**
3. `FWHM.RT` — only regress out if NOT biologically confounded with your treatment
   (check the covariate-correlation heatmap first).
4. `MS1.Signal` — usually biologically confounded; regress only after checking.

---

## Activity Score Gene Sets

Gene names in the "Activity Scoring" tab must match the **protein group identifiers**
in your PG matrix — typically HGNC gene symbols (e.g., `COX1`, `NDUFB8`, `MYH7`).

Example gene sets from your CM cardiomyocyte work:

```
OxPhos:    COX1, COX2, NDUFB8, NDUFB9, ATP5A1, ATP5B, UQCRC1
Myogenesis: MYH6, MYH7, TNNT2, ACTC1, MYL2, MYL3, TNNC1
EMT:        VIM, FN1, CDH2, SNAI1, ZEB1, TWIST1
UPS_Stress: PSMA1, PSMA2, PSMA3, PSMB5, PSMC1, PSMD1, UBB, UBA1
ISR:        EIF2AK4, ATF4, DDIT3, PPP1R15A, SLC7A5
```

---

## Typical Workflow

```
Upload PG Matrix + Annotation + report.stats
  ↓
Tab 1 — Upload & QC
  • Review n_proteins, pct_detected, FWHM.RT distributions
  • Set filters: min_proteins ≥ 100, protein detection ≥ 10%

Tab 2 — Preprocessing
  • Normalize (median total count) → log1p
  • Regress: [log_n_proteins, Normalisation.Instability]
  • Run PCA (50 PCs)
  • Check covariate correlation heatmap (confirm regression worked)
  • Optional: run Harmony if batch (plate) effect is present

Tab 3 — Embedding & Clustering
  • Run Neighbors (k=15, 30 PCs) → UMAP → Leiden (res=0.5)
  • Explore cluster composition by cell_line / condition / batch

Tab 4 — Differential Expression
  • Group by: concentration (or condition, cm_subtype_v2)
  • Reference: 0uM (or Control)
  • Method: Wilcoxon (min 3 cells, 25% detection filter)
  • View volcano, download DE table, generate heatmap

Tab 5 — Activity Scoring
  • Paste gene sets (OxPhos, Myogenesis, EMT, ISR, UPS_Stress)
  • View per-cell scores on UMAP
  • Compare score distributions across dose/subtype/cell-line groups
```
