# Pro-Visualize 2.0 — Production Readiness Plan

> **Status:** Active — execution started 2026-05-03 on `feature/p0-foundation`.
> **Branching baseline:** `develop` (after merge of `feature/single-cell-integration`).
> **Source:** Originally drafted via interactive planning session; this file in the repo is the source of truth going forward — update it as phases land.

---

## 🟢 Current Working State (last updated 2026-05-10 by Sonnet 4.6)

**Active branch:** `feature/p7-deploy-ready` — P7 implementation complete, PR pending → develop.

### Phase completion summary

| Phase | Status | Merge commit / tag | Branch |
|-------|--------|-------------------|--------|
| P0 — Foundation | ☑ Done | PR #1 merged → develop; `v2.0.0-p0` tagged | `feature/p0-foundation` |
| P1 — Report Builder | ☑ Done | PR #2 merged → develop; `v2.0.0-p1` tagged | `feature/p1-report-builder` |
| P2 — Comparative | ☑ Done | PR #3 merged → develop; `v2.0.0-p2` tagged | `feature/p2-comparative` |
| P3 — Quantification | ☑ Done | PR #4 merged → develop; `v2.0.0-p3` tagged | `feature/p3-quant` |
| P4 — Dilution Series | ☑ Done | PR #5 merged → develop; `v2.0.0-p4` tagged | `feature/p4-dilution` |
| P5 — QC (DIA) | ☑ Done | PR #6 merged → main (sync'd to develop); `v2.0.0-p5` tagged | `feature/p5-qc` |
| P6 — SCP Polish | ☑ Done | PR #7 merged → develop; `v2.0.0-p6` tagged | `feature/p6-scp-polish` |
| P7 — Deploy-Ready | ⏳ PR open | PR #8 → develop; `v2.0.0` tag on release | `feature/p7-deploy-ready` |
| P8 — Publication Hardening | ⏳ PR pending | umbrella → develop; `v2.0.1` on release | `feature/p8-publication-hardening` (subs: p8a–p8d) |

### What's done in P2 (Comparative module upgrade)

- ✅ **`@handle_plotting_errors`** added to all 5 plot methods in `comparative_visualizer.py`: `plot_volcano`, `plot_comparative_heatmap`, `plot_expression_violin`, `plot_enrichment_manhattan`, `plot_enrichment_dotplot`.
- ✅ **`safe_render` isolation** — all 6 inner tabs (overview/select/volcano/heat/expr/path) wrapped in `try/except` with scoped error cards; failure in one tab cannot crash the others.
- ✅ **Heatmap migrated to `MplPlotManager`** — replaces raw `st.image()`. Heatmap now has title/figsize/dpi editing, PNG download, and "Add to Report". `plot_comparative_heatmap` accepts `title`/`figsize`/`dpi` kwargs so edits actually re-render.
- ✅ **Enrichment threshold presets** — "Most stringent" (FDR<0.01, |FC|>1.5), "Standard" (0.05/1.0), "Exploratory" (0.10/0.5) surfaced in Selection tab via selectbox + Apply button. Uses `SIGNIFICANCE_PRESETS` from `config/plot_configs.py`.
- ✅ **Enrichr retry + timeout** — `HTTPAdapter(max_retries=Retry(total=3, backoff_factor=0.5, status_forcelist=[429,500,502,503,504]))` + `timeout=30` on both POST and GET. Actionable error message on failure: "Enrichr API may be down — try again or use the offline gene-list export."
- ✅ **`HUMAN_TRANSCRIPTION_FACTORS` moved to config** — removed class attribute from `comparative_visualizer.py`; now imported from `config/plot_configs.py` (single source of truth).
- ✅ **"Add to Report" wired** — `.module = "comparative"` set on all 5 PlotManager/MplPlotManager instances (volcano, violin, manhattan, dotplot ×N, heatmap). Figures now appear under "comparative" in the Report tab, not "unknown".

### What's done in P3 (Quantification module)

- ✅ **Venn/UpSet/Dendrogram/IntensityDist migrated to `MplPlotManager`** — title editing, DPI, PNG download, and "Add to Report" wired in.
- ✅ **Correlation Matrix tab implemented** — Pearson/Spearman with hierarchical clustering reorder, upper triangle masked, values annotated for ≤30 samples.
- ✅ **CV vs Intensity scatter** — new plot tab showing protein-level CV% vs log2 mean intensity.
- ✅ **`HUMAN_TRANSCRIPTION_FACTORS` removed from class** — imported from `config/plot_configs` (single source of truth).
- ✅ **`to_hex()` replaces fragile RGB parser** — `pc.unlabel_rgb` removed; `utils/helpers.to_hex()` used for group color picker.
- ✅ **All PlotManager instances have `.module = "quant"`** — "Add to Report" routes correctly.
- ✅ **Gene symbol coverage banner** — `sanity.gene_resolution_report` / `render_validation` wired at top of dashboard.

### What's done in P4 (Dilution module)

- ✅ **Configurable column names** — `sample_col`, `concentration_col`, `replicate_col`, `group_col` params added to `DilutionSeriesVisualizer.__init__`; metadata normalized to internal names at init so all downstream code is unchanged. UI text inputs in `_upload_data_section` under "Column Configuration" expander.
- ✅ **`DEVIATION_BUCKETS` configurable** — two sliders (`good_thresh`, `warn_thresh`) in "Global Plot Settings" expander; passed through to `plot_relative_abundance_ratios` and `_classify_deviation_color`.
- ✅ **CSV exports** — per-protein R² table (`get_r2_table()`), CV-by-concentration matrix (`get_cv_by_concentration_matrix()`), completeness summary (`get_completeness_summary()`) — all downloadable from their respective tabs.
- ✅ **R² histogram + ranked R² table** — new "📉 Linearity (R²)" tab with `plot_r2_histogram()`, median R² caption, sortable table, and CSV download.
- ✅ **LOD/LOQ plot + table** — new "🔬 LOD/LOQ" tab with `plot_lod_loq()` (scatter LOD vs R²), `get_lod_loq_table()` (slope-based: LOD = 2^(log₂C_min + 3.3σ/slope)), sortable table with in-range indicator, CSV download.
- ✅ **Sanity checks** — `run_sanity_checks()` checks for non-positive concentrations, duplicate replicate-concentration pairs, uneven geometric ratios, and negative-slope proteins; displayed in collapsible banner at top of dashboard.
- ✅ **ReportBuilder wired** — `.module = "dilution"` set on all 10 `PlotManager` instances.

### What's done in P5 (DIA QC module)

- ✅ **PlotManager wired into all 10 DIA QC plots** — IM control/drift, RT control/drift/pred-error/peak-width/elution, mass-accuracy dist/sentinel/trend. All have PNG/SVG/HTML export + "Add to Report".
- ✅ **`.module = "dia_qc"`** on every PlotManager — Report tab correctly categorises DIA QC figures.
- ✅ **Configurable σ-threshold slider** (1–3σ, default 2) in IM and RT control chart tabs; orange dotted ±Nσ warning lines drawn at the chosen threshold.
- ✅ **`QC_THRESHOLDS` and `DIA_RUN_NAME_PATTERNS` imported from `config/plot_configs`** — config-driven defaults, no hardcoded values in `DiaQcVisualizer`.
- ✅ **SCP bonus fix** — Windows backslash paths in `scp_visualizer.py` normalised before `os.path.basename()`, fixing annotation↔PG-matrix matching for DIA-NN files from Windows.

### What's done in P6 (SCP module polish)

- ✅ **`_pm()` / `_mpl_pm()` helpers** — all 14+ SCP PlotManagers pre-stamped with `.module = "scp"`; every SCP figure now routes to the Report tab.
- ✅ **DE heatmap → `MplPlotManager`** — migrated from bare `st.image()`. Gains title editing, DPI, PNG download, "Add to Report". `plot_de_heatmap()` accepts `title`/`figsize`/`dpi` kwargs.
- ✅ **Session save / load** — "💾 Save / Load Session" expander in Preprocessing tab. Serialises full `SCPVisualizer` (AnnData + `pp_state`) via pickle. Scientists can resume analysis across sessions.
- ✅ **Expression Overlay UMAP** — new "🎨 Expression Overlay" tab in Embedding. Select any protein from a dropdown; UMAP coloured by log-normalised expression (Viridis). `plot_expression_umap()` + `get_protein_names()` added to `SCPVisualizer`.
- ✅ **`scp_expr_umap` and `scp_de_heatmap`** added to `PLOT_KEYS` so cached figures clear on pipeline resets.

### What's done in P7 (Deploy-Ready)

- ✅ **`Dockerfile`** — Python 3.11-slim + `libgbm1`/`libasound2`/`libxshmfence1` for kaleido Chromium; `/_stcore/health` healthcheck; exposes 8501.
- ✅ **`docker-compose.yml`** — single service with `./data` + `./reports` volume mounts, `env_file: .env.docker`, `restart: unless-stopped`.
- ✅ **`.env.docker`** — committed with safe defaults (`STREAMLIT_SERVER_MAX_UPLOAD_SIZE=500`).
- ✅ **`.dockerignore`** — excludes `.git`, `.venv`, `__pycache__`, `.claude/`, `data/`, `reports/`, etc.
- ✅ **`Makefile`** — `build`, `run`, `stop`, `logs`, `test`, `lint`, `clean` targets.
- ✅ **`.streamlit/config.toml`** — `maxUploadSize=500`, `enableXsrfProtection=true`, `headless=true`.
- ✅ **`tests/`** — 50 pytest smoke tests across all 5 visualizers + ReportBuilder + sanity helpers:
  - `test_dilution_visualizer.py` (7 tests) — instantiation + 5 plot methods + LOD/LOQ
  - `test_quant_visualizer.py` (7 tests) — instantiation + PCA, correlation, CV vs intensity, rank-order
  - `test_comparative_visualizer.py` (4 tests) — instantiation + volcano, violin, heatmap
  - `test_dia_qc_visualizer.py` (5 tests) — parquet loading, metadata extraction, error cases
  - `test_scp_visualizer.py` (5 tests) — instantiation, QC metrics, preprocess, PCA, plot
  - `test_report_builder.py` (9 tests) — add/remove/reorder, HTML export, ZIP structure
  - `test_sanity.py` (13 tests) — validate_columns, check_sample_alignment, missingness, gene resolution
- ✅ **Welcome tab redesign** — quick-start cards per module with expected file-format table; replaced static markdown with `st.container(border=True)` cards.

**Next steps (release):**
1. Open PR #8: `feature/p7-deploy-ready` → `develop`
2. Smoke test: `pytest tests/ -q` passes (50/50)
3. Merge PR → develop, tag `v2.0.0-p7`
4. Merge `develop` → `main`, tag `v2.0.0`

### What's done in P8 (Publication Hardening — pre-paper audit)

Branch `feature/p8-publication-hardening`, cut from `feature/p7-deploy-ready`
(P7 was not yet on `develop`, so the Docker/tests work being hardened lived
only on the P7 branch). Four reviewable sub-checkpoints, each merged `--no-ff`:

- **P8a — Reproducible build** (`feature/p8a-repro-deploy`, commit `97035f7`):
  pinned all loose deps to exact versions (scanpy 1.11.5, anndata 0.12.11,
  gseapy 1.2.1, igraph 1.0.0, leidenalg 0.11.0, jinja2 3.1.6); enabled
  harmonypy 2.0.0 (guarded import). New `requirements-dev.txt`
  (pytest 9.0.3, ruff 0.15.16). Dockerfile pinned to `python:3.11.15-slim`,
  non-root `appuser`, full kaleido/Chromium libs. New `.github/workflows/ci.yml`
  (lint + pytest; docker build asserting non-root + kaleido PNG export).
  `pyproject.toml` ruff config (F + E9). Makefile test/lint targets made
  CI-safe. Lint hygiene + `.DS_Store` gitignore.
- **P8b — Enrichment background toggle** (`feature/p8b-enrichment-background`,
  commit `bed7f14`): the one scientific-correctness fix. Enrichment can now use
  the study's **detected proteins** as the statistical background (default)
  instead of the whole genome, via Enrichr Speedrichr endpoints (comparative)
  and gseapy `background=` (SCP); whole-genome remains an explicit option.
  UI selector + background-N caption in both modules. Network-mocked tests.
- **P8c — Robustness polish** (`feature/p8c-robustness-polish`, commit
  `e93dc94`): new `utils/logging_config.py` (root logging, called from
  `app.py`); replaced silent `except: pass` in `plot_manager.py` with logging;
  empty/degenerate-data guards on Venn/UpSet/clustering. Audit note: error
  isolation was already comprehensive (PlotManager try/except + QC tab-method
  decorators + tab-level `safe_render`), so blanket-decorating was unnecessary.
- **P8d — Reproducibility report** (`feature/p8d-repro-report`, commit
  `9a62e25`): report exports now embed a Methods & Reproducibility block
  (app/git/python/package versions); ZIP gains `provenance.json` and nests
  params under `parameters.json`; `add_figure` guards non-serializable params.
- **P8 (umbrella) — Methods doc**: `METHODS.md` with ready-to-paste language
  for the paper (statistics, transformation/imputation, LOD/LOQ, SCP, QC,
  enrichment background).

Test suite grew 50 → 60; `ruff check .` clean. Docker-level gates (build,
non-root, kaleido PNG) are exercised by CI (local docker daemon was down).

**P8 release steps:**
1. Open PR: `feature/p8-publication-hardening` → `develop` (after P7 merges).
2. Confirm CI green (lint + pytest 60/60 + docker build/import/PNG).
3. Merge → develop; merge `develop` → `main`; tag `v2.0.1`.

### Environment notes (critical for resumption)

- **Python venv** at `.venv/`: `source .venv/bin/activate`. `kaleido==0.2.1` installed.
- **`python` not on PATH unactivated** — system python3 is 3.9.6. Always activate venv first.
- **Demo data** at `/Volumes/VanEykJLab-Files/ByPerson/Ajay/Pro-Vizualize-2.0-Demo/` — see `reference_demo_data.md` memory.
- **`.claude/settings.local.json`** has harness config drift in working tree — leave it; do not stage.

### What P2 looks like (next phase)

Cut `feature/p2-comparative` from `develop`. P2 deliverables (Section 3.4):
- Add `@handle_plotting_errors` + `safe_render` to all Comparative tabs (currently no error handling at all).
- Migrate heatmap to `MplPlotManager`.
- Surface enrichment presets: "Most stringent" (FDR<0.01, |log2FC|>1.5), "Standard" (0.05/1.0), "Exploratory" (0.1/0.5).
- Wrap Enrichr/g:Profiler calls in `HTTPAdapter(max_retries=3)` + 30s timeout; actionable error on failure.
- Move `HUMAN_TRANSCRIPTION_FACTORS` import to `config/plot_configs.py` (remove duplicate in `comparative_visualizer.py:22-30`).
- Hook "Add to Report" into `ReportBuilder` for all Comparative plots.
- Critical files: `modules/comparative_module.py`, `visualizations/comparative_visualizer.py`.

---

## Context

Pro-Visualize is a Streamlit app helping proteomics scientists explore their data across 6 analysis types: QC (DIA + Targeted), Dilution Series, Quantification, Comparative, Pathway Enrichment, and Single-Cell Proteomics (SCP). The SCP module (most recently built) sets a clear quality bar — `PlotManager`-based per-figure editing, granular `try/except` per stage, `pp_state` workflow gating, AnnData layer-based data lineage, and rich user feedback. The older modules drift from this standard in five concrete ways:

1. **Inconsistent exception isolation** — an error in one tab can cascade or silently fail (Comparative has no decorator at all; Quant uses it on only one method).
2. **Incomplete figure customization** — `PlotManager` is used in some modules but skipped for matplotlib outputs (Venn, UpSet, dendrograms, heatmaps); QC modules don't use it at all.
3. **No report generation** — there's no way to bundle the figures a scientist ends a session with into a deliverable. Only DIA-QC has any export (raw CSV).
4. **No static image export** — `kaleido` isn't installed; no PNG/SVG download buttons on Plotly figures.
5. **Sparse sanity checks & validation drift** — early modules expect exact column names; the SCP module shows a much better pattern (parameterized columns, fallbacks, live previews, "X of Y matched" feedback).

End goal: a polished, deployable app where (a) every figure is high quality and editable, (b) any failure is contained to its tab and shows a useful message, (c) a scientist can leave with an HTML report capturing their edited figures + parameters + interpretation notes, and (d) sanity checks catch bad input before it produces nonsense plots.

## Goals

- **Consistent quality bar** across all modules matching the SCP reference (`scp_module.py:37-126`, `scp_module.py:1057-1093`).
- **Editable figures** — every plot supports title, axis labels, height, marker size, colors, theme. Matplotlib figures get an equivalent (or convert to Plotly where reasonable).
- **HTML report builder + ZIP bundle export** — per module, "Add to Report" queues a figure; the Report tab produces (a) a self-contained interactive HTML and (b) a ZIP bundle containing PNG/SVG/HTML per figure plus a `parameters.json` snapshot. Two outputs, one queue.
- **Robust isolation** — failure in one tab/figure never affects others. A standardized `safe_render()` wrapper catches, logs, and shows a recoverable error state.
- **Sanity checks** — input validation, post-computation diagnostics ("23 of 47 proteins resolved to gene symbols"), and warning banners for suspicious results.
- **Deployable** — pinned deps, clear setup instructions, secrets via `st.secrets`, no hardcoded paths.

## Non-Goals (this plan)

- Not rewriting visualizations from scratch — we upgrade in place.
- Not replacing `gprofiler-official` / `gseapy` / Enrichr — current pathway sources stay.
- Not adding authentication or multi-user workflows — out of scope for v2.0.
- Not building a PDF report path until HTML report is solid (HTML preserves Plotly interactivity, which is the bigger value to scientists). The ZIP bundle below covers users who want raw assets.
- Not adding multi-user auth or session sharing — single-user-per-instance for v2.0.

## Architecture Overview

Three layers, in build order:

```
┌─────────────────────────────────────────────────────────────────┐
│ Layer 3: Per-Module Upgrades                                     │
│   QC (DIA, Targeted) │ Dilution │ Quant │ Comparative │ SCP polish│
└─────────────────────────────────────────────────────────────────┘
                              ▲
┌─────────────────────────────────────────────────────────────────┐
│ Layer 2: Cross-Cutting Features                                  │
│   Report Builder │ Static Export │ Sanity Checks │ Tab Isolation │
└─────────────────────────────────────────────────────────────────┘
                              ▲
┌─────────────────────────────────────────────────────────────────┐
│ Layer 1: Foundation (refactor existing utilities)                │
│   utils/plot_manager.py │ utils/helpers.py │ config/plot_configs │
└─────────────────────────────────────────────────────────────────┘
```

Build Layer 1 → Layer 2 → Layer 3. Each layer is independently testable and ships value before the next starts.

---

## Layer 1: Foundation — Refactor Existing Utilities

### 1.1 Extend `utils/plot_manager.py` (currently 108 lines)

Current state: handles Plotly only; supports title, height, marker size, axis labels (`utils/plot_manager.py:49-108`). Loses custom edits when re-generated.

Changes:
- **Persist edits across regeneration**: store edit dict in `st.session_state[f"{key}_edits"]`, reapply after `_generate_plot`.
- **Add color/theme override per-plot** (in addition to global): font family, font size, legend position, x/y axis ranges (auto/manual), gridlines toggle.
- **Add static export buttons**: PNG, SVG, HTML — using `kaleido` for static and `fig.to_html()` for interactive.
- **Add "Add to Report" button** — calls `report_builder.add_figure(key, fig, title, notes)`.
- **Matplotlib support**: a parallel `MplPlotManager` (or unified subclass) for figures returned as `BytesIO` — supports title, dpi, figsize via re-rendering with stored parameters; static-only export.

Critical files to modify:
- `utils/plot_manager.py` — extend
- `utils/helpers.py` — add `safe_render(label, fn, *args, **kwargs)` wrapper that catches, logs, displays an inline error card with "Retry" button. This is the new isolation primitive.

### 1.2 Centralize config in `config/plot_configs.py` (currently empty)

Move scattered hardcoded values:

```python
# config/plot_configs.py
THEMES = {"Standard White": "plotly_white", "Dark Mode": "plotly_dark", "Minimal": "simple_white"}
COLORBLIND_PALETTE = [...]  # Okabe-Ito 8-color
HUMAN_TRANSCRIPTION_FACTORS = {...}  # currently duplicated in quant_visualizer.py:29 and comparative_visualizer.py:22
SIGNIFICANCE_DEFAULTS = {"fdr": 0.05, "log2fc": 1.0, "pval": 0.05}
QC_THRESHOLDS = {"cv_low": 20, "cv_high": 30, "control_sigma": 2}
DEVIATION_BUCKETS = {"good": 0.2, "warn": 0.5}  # dilution_series.py:114-119
COLUMN_FALLBACKS = {
    "protein_id": ["Protein", "ProteinIds", "Accession"],
    "gene": ["Gene Name", "Gene", "Gene names", "Genes"],
    ...
}
```

This single file becomes the source of truth — referenced from every module/visualizer.

### 1.3 Add `utils/sanity.py`

Helper functions for repeated patterns:
- `validate_columns(df, required, optional=None) -> ValidationResult` — returns missing list and warnings.
- `check_sample_alignment(protein_df, annotation_df, sample_col)` — reports unmatched samples.
- `check_value_ranges(series, expected_range, name)` — flags suspicious distributions (e.g., negative intensities, all-zero proteins).
- `summarize_missingness(df) -> dict` — used by every module's data overview.
- `gene_resolution_report(df, gene_col)` — "X/Y proteins have a gene symbol" caption like SCP does (`scp_module.py:209-226`).

### 1.4 Pin new dependencies

Add to `requirements.txt`:
- `kaleido==0.2.1` — Plotly static image export (PNG/SVG/PDF from `fig.to_image()`).
- `jinja2>=3.1.0` — already a Streamlit transitive dep, but explicit-pin for report templates.

(Defer `weasyprint` until PDF report path is requested — HTML covers the immediate need and avoids native dep complications on deploy.)

---

## Layer 2: Cross-Cutting Features

### 2.1 Report Builder — `utils/report_builder.py` (new) + `templates/report.html.j2` (new)

Two output formats from one queued state: **interactive HTML** and **per-figure ZIP bundle**. The user picks one or both at export time.

Design:

```python
class ReportBuilder:
    """Per-session report state, lives in st.session_state.report"""

    def add_figure(self, module: str, key: str, fig, title: str, params: dict, notes: str = "") -> None
    def add_table(self, module: str, key: str, df, caption: str) -> None
    def add_section(self, module: str, heading: str, body_md: str) -> None
    def remove(self, module: str, key: str) -> None
    def reorder(self, module: str, keys: list[str]) -> None
    def render_preview(self) -> None       # streamlit preview pane

    # Two export paths:
    def export_html(self) -> bytes         # Jinja2 → self-contained HTML w/ interactive Plotly
    def export_zip(self) -> bytes          # ZIP: figs/<module>_<key>.{png,svg,html} + parameters.json + notes.md
```

UI surfaces (every module gets):
1. **Per-figure** "Add to Report" button (in `PlotManager` editor).
2. **Per-module** "Report Items" expander listing what's queued, with reorder/remove.
3. **Top-level Report tab** (new tab in `app.py`) — full preview, edit notes/headings, **two download buttons**: "Download Interactive HTML" and "Download ZIP Bundle".

**Interactive HTML template** (Jinja2): single-file HTML with `include_plotlyjs='cdn'`, one section per module, captions, parameters table, timestamp, app version. Plotly figures use `fig.to_html(full_html=False)` for interactivity; matplotlib figures embed as inline base64-PNG.

**ZIP bundle structure**:
```
report.zip
├── parameters.json          # all plot params keyed by module/key
├── notes.md                 # user-entered narrative text
├── manifest.json            # ordered list of figures with metadata
└── figures/
    ├── comparative_volcano.png
    ├── comparative_volcano.svg
    ├── comparative_volcano.html
    ├── quant_pca.png
    └── ...
```

Both paths share the same queue and `add_figure` API — no duplication.

### 2.2 Tab Isolation Primitive

Wrap each top-level tab body in `safe_render`:

```python
# app.py
with tab_quant:
    safe_render("Quantification", render_quant_module)
```

Internally, `safe_render` catches all exceptions, logs full traceback, and renders a friendly error card with the exception type, message, traceback in expander, "Retry" button, and "Reset module" button (clears that module's session_state keys).

Plus: every individual figure rendered through `PlotManager` gets the same isolation — one broken plot doesn't break the tab.

### 2.3 Static Image Export

Once `kaleido` is installed, add to `PlotManager.render_plot_and_editor()`:

```python
col1, col2, col3 = st.columns(3)
col1.download_button("⬇ PNG", fig.to_image(format="png", scale=2), file_name=f"{key}.png")
col2.download_button("⬇ SVG", fig.to_image(format="svg"), file_name=f"{key}.svg")
col3.download_button("⬇ HTML", fig.to_html(include_plotlyjs="cdn"), file_name=f"{key}.html")
```

Matplotlib figures get the existing `BytesIO` PNG download.

### 2.4 Sanity Check Framework

Each module's data load path runs `sanity.validate_*` and surfaces results in a top-level "Data Quality" expander:

```
✓ All required columns present (5/5)
⚠ 12 of 847 proteins lack gene symbols
⚠ 3 samples in protein file not found in annotation: SAMPLE_X, SAMPLE_Y, SAMPLE_Z
✗ Negative values detected in intensity matrix (n=47) — log transform will fail
```

Severity drives behavior: `✗` blocks computation; `⚠` warns but allows; `✓` confirms.

---

## Layer 3: Per-Module Upgrades

### 3.1 QC Module (DIA + Targeted) — Highest Risk, Brittle Today

**Current state:**
- `dia_qc_tab.py` — uses `@handle_plotting_errors` but no PlotManager. Hardcoded date regex (`DiaQcVisualizer.py:109-115, 128-131`) breaks on non-standard run names.
- `targeted_qc_tab.py` — same issues plus single regex pattern for `Replicate Name` (`targettedQCVisualization.py:115-123`), hardcoded peptide-count defaults, no metadata review step.

**Changes:**
1. Wire PlotManager into all DIA + Targeted plots. Convert plot methods to `(self, **kwargs) -> go.Figure` signature.
2. Replace hardcoded regex with configurable patterns from `config/plot_configs.py` + a UI fallback "I can't parse — let me map columns manually."
3. Add metadata review step to Targeted module to match DIA's pattern (`dia_qc_tab.py:125-126` workflow gate).
4. Convert sentinel peptide selection into a saved/named workflow (preset peptides per instrument type).
5. Add CV thresholds and σ-bands to control charts as configurable inputs (currently `±1σ, ±2σ, ±3σ` is hardcoded).
6. Add CSV download buttons for QC summary tables (parity with DIA's existing one).
7. Hook into ReportBuilder.

**Critical files:**
- `modules/qc_tabs/dia_qc_tab.py`
- `modules/qc_tabs/targeted_qc_tab.py`
- `visualizations/DiaQcVisualizer.py`
- `visualizations/targettedQCVisualization.py`

### 3.2 Dilution Series Module — Closest to SCP Pattern Today

**Current state:** Already uses PlotManager for all 8 plots (`dilution_module.py:142-325`); has theme + colorblind toggle (`dilution_module.py:45-125`). Missing: deviation thresholds hardcoded, no exports, "Concentration" column name not flexible.

**Changes:**
1. Surface `DEVIATION_BUCKETS` as user-configurable sliders.
2. Make `Concentration`/`Replicate`/`Group` column names configurable (currently hardcoded — `dilution_series.py:30-32`).
3. Add CSV export for: per-protein R² fit table, CV-by-concentration matrix, completeness summary.
4. Add R² histogram + ranked R² table to highlight proteins with poor linearity.
5. **New plot**: LOD/LOQ estimation per protein using common slope-based methods (CCβ-style or 3.3·σ/slope) with a sortable summary table.
6. **New sanity check**: warn if `Concentration` is non-monotonic per replicate, or if log-fit coefficient is negative when slope should be positive.
7. Hook into ReportBuilder.

**Critical files:**
- `modules/dilution_module.py`
- `visualizations/dilution_series.py`

### 3.3 Quantification Module

**Current state:** Mixed Plotly + matplotlib outputs (`quant_visualizer.py`); only Plotly plots use PlotManager. TF list duplicated (`quant_visualizer.py:29-37` vs `comparative_visualizer.py:22-30`). Correlation matrix is a stub (`quant_module.py:321`). Color picker logic fragile (`quant_module.py:114-118`).

**Changes:**
1. Migrate Venn/UpSet/Dendrogram to `MplPlotManager` so they get title editing + PNG download + Add to Report.
2. **Implement Correlation Matrix tab**: sample×sample Pearson/Spearman with hierarchical clustering, mask upper triangle, annotate values for ≤30 samples.
3. Move TF list to `config/plot_configs.py` (single source of truth).
4. Replace fragile RGB-to-hex parsing with `utils/helpers.py:to_hex(color)`.
5. Use `sanity.gene_resolution_report` to show "Gene Name" column status with same UX as SCP.
6. **New plot**: protein-level CV vs intensity scatter (helps users decide on filters before downstream).
7. Hook into ReportBuilder.

**Critical files:**
- `modules/quant_module.py`
- `visualizations/quant_visualizer.py`

### 3.4 Comparative Module

**Current state:** No `@handle_plotting_errors` decorator usage at all. Heatmap method bypasses PlotManager (`comparative_module.py:193-197`). TF list duplicated. Enrichr API has no timeout/retry (`comparative_module.py:252-254`). Default thresholds are lenient and unguided.

**Changes:**
1. Add `@handle_plotting_errors` + `safe_render` to all tabs.
2. Migrate heatmap to `MplPlotManager`.
3. Surface enrichment defaults with guidance: "Most stringent" (FDR<0.01, |log2FC|>1.5), "Standard" (0.05/1.0), "Exploratory" (0.1/0.5).
4. Wrap Enrichr / g:Profiler calls in `requests.adapters.HTTPAdapter(max_retries=3)` + 30s timeout; on failure show actionable error ("API may be down — try again or use offline gene-list export").
5. **Cross-comparison support**: if user uploads multiple comparative files, allow side-by-side volcano + shared/unique protein UpSet.
6. **New plot**: rank-rank hypergeometric overlap (RRHO) for two comparisons (commonly requested by reviewers).
7. Move TF list to config.
8. Hook into ReportBuilder.

**Critical files:**
- `modules/comparative_module.py`
- `visualizations/comparative_visualizer.py`

### 3.5 SCP Module — Polish

**Current state:** Already gold-standard architecturally. Improvements are incremental.

**Changes:**
1. Persist PlotManager edits across regeneration (currently lost — line `utils/plot_manager.py:25-33`).
2. Add static exports (PNG/SVG/HTML) to all SCP plots.
3. Add "Save AnnData state" / "Load AnnData state" — pickle the `adata` after preprocessing for resume-on-next-session.
4. Add per-cell DE inspection: click a UMAP cluster → show top markers in a side pane.
5. Add expression-overlay UMAP (continuous color by single-protein abundance, multi-select).
6. Hook into ReportBuilder.
7. Surface Harmony / Leiden parameters in the report so methods are reproducible.

**Critical files:**
- `modules/scp_module.py`
- `visualizations/scp_visualizer.py`

---

## Suggested Improvements (Beyond Original Request)

### A. Welcome / Landing Page

Currently the Welcome tab is static (`app.py:28-36`). Replace with:
- Quick-start cards per module ("Have a DIA-NN report.tsv? Click here.")
- Sample data downloads (one tiny fixture per module — ~50 proteins, ~10 samples).
- Recent session state ("Continue last analysis?")

### B. Data Loading Wizard

Common pain point: users don't know which column maps to what. Add a one-time "Detect columns" button per module:
- Reads first uploaded file.
- Uses `COLUMN_FALLBACKS` from config to auto-fill text inputs.
- Shows a preview table with detected column → role mapping.
- Lets user override before "Confirm."

### C. Logging & Telemetry

Add `utils/logging_config.py` setting up a rotating file handler in `~/.pro_visualize/logs/`. On error cards, show "Log file: ..." for debug. No PII; just analysis events + tracebacks.

### D. Lightweight Tests

Add `tests/` with smoke tests:
- `test_visualizers.py` — load fixture data, instantiate each visualizer, run each plot method, assert non-empty output.
- `test_sanity.py` — unit test the validation helpers.
- `test_report_builder.py` — assert HTML output contains expected sections.

Aim for ~30 minutes runtime, run via `pytest tests/`.

### E. Deployment Hygiene (Self-hosted Docker target)

Confirmed deployment target: **self-hosted Docker** (internal server / cloud VM). This shapes the work as follows:
- **`Dockerfile`** — Python 3.11 slim base + system deps for `kaleido` (Chromium prerequisites: `libgbm1`, `libasound2`, `libxshmfence1`), copy app, install requirements, expose 8501, run `streamlit run app.py --server.address=0.0.0.0`.
- **`docker-compose.yml`** — single service mounting `./data` for sample fixtures and `./reports` for any persisted outputs; env file for tunables (max upload size, log level).
- **`.streamlit/config.toml`** — `[server]` `maxUploadSize=500` (MB) since Docker has no Streamlit Cloud cap; `[theme]` defaults; `enableXsrfProtection=true`.
- **`.dockerignore`** — exclude `.git`, `__pycache__`, `.claude/`, `tests/fixtures/`.
- **`st.secrets`** — for any future API keys (Enrichr/gprofiler are open today; LLM keys for the "Pro-Viz Chat" tab would go here).
- **`make` targets**: `make build`, `make run`, `make test`, `make lint`, `make logs`.
- **Health check endpoint** — Streamlit's `/_stcore/health` is enough; document it in README for reverse-proxy setup.
- **Resource notes** — call out in README that SCP AnnData on large datasets needs ≥4 GB RAM; Docker compose can set memory limits.

### F. Optional: Plotly→Streamlit Renderer Upgrade

Streamlit 1.36 supports `st.plotly_chart(use_container_width=True, theme=None)`. Confirm we're using `theme=None` everywhere so user themes win over Streamlit's default override.

---

## Phasing & Estimated Effort

The plan is sized so each phase ships value independently — you can stop after any phase and the app is still better than today.

| Phase | Status | Scope | Effort | Ship Value |
|-------|--------|-------|--------|-----------|
| **P0** | ☑ Done (`v2.0.0-p0`) | Foundation: extend `PlotManager`, `safe_render`, config, sanity, kaleido | 1-2 days | Tab isolation + static exports work everywhere |
| **P1** | ☑ Done (pending merge) | Report Builder (HTML + ZIP bundle) + template + Report tab | 1-2 days | Scientists leave with interactive HTML or raw asset ZIP |
| **P2** | ☐ Pending | Comparative module upgrade (most user-facing, most fragile) | 1 day | Critical path stabilized |
| **P3** | ☐ Pending | Quantification module upgrade + correlation matrix | 1 day | Feature-complete |
| **P4** | ☐ Pending | Dilution module upgrade + LOD/LOQ feature | 0.5-1 day | Adds genuine new science |
| **P5** | ☐ Pending | QC module upgrade (DIA + Targeted) | 1.5-2 days | Brittle code path hardened |
| **P6** | ☐ Pending | SCP polish + state persistence | 1 day | Reproducibility |
| **P7** | ⏳ PR open | Welcome page, wizard, Dockerfile/compose, tests | 1-2 days | Deploy-ready (self-hosted Docker) |

Total: ~8-12 working days of focused effort. P0+P1 first because everything else inherits from them.

> **How to update phase status as work lands:** mark each phase ☑ Done with the merge commit / tag (`v2.0.0-p<N>`) when its PR merges to `develop`. Move the ⏳ marker to the next phase being worked.

---

## Branching & Checkpoint Strategy

Every phase ships as its own branch with a working, testable checkpoint. No phase merges to `develop` until its smoke test passes. This means rollback is one branch deletion, and at any point `develop` is a runnable app.

### Step 0: Pre-Plan Cleanup (DONE on 2026-05-03)

Original branch was `feature/single-cell-integration` with uncommitted SCP work.

1. ☑ Committed pending SCP work (`f26d708 — feat(scp): add enrichment filters and multi-gene mapping`).
2. ☑ Pushed `feature/single-cell-integration` to origin.
3. ☑ Merged `feature/single-cell-integration` → `develop` with `--no-ff` (`0d27126`); `develop` pushed to origin.
4. ☑ Cut `feature/p0-foundation` from `develop`. All subsequent phase branches will be cut from `develop`.

### Branch Naming Convention

```
feature/p0-foundation
feature/p1-report-builder
feature/p2-comparative
feature/p3-quant
feature/p4-dilution
feature/p5-qc
feature/p6-scp-polish
feature/p7-deploy-ready
```

### Per-Phase Workflow (every phase follows this)

```
1. git checkout develop && git pull
2. git checkout -b feature/p<N>-<name>
3. <do the work — see "Handover Criteria" per phase below>
4. Run smoke test for phase
5. Open PR to develop
6. Self-review or co-review
7. Merge PR (squash) → develop
8. Tag the merge: git tag v2.0.0-p<N>
9. Repeat for P<N+1>
```

### Definition of Done — applies to EVERY phase

A phase is done when ALL of these are true:
- ✅ Code compiles and `streamlit run app.py` launches without error.
- ✅ All previously-working tabs still load (no regressions).
- ✅ Phase-specific smoke test (listed below) passes manually.
- ✅ At least one new feature from this phase is exercised end-to-end in the running app.
- ✅ No `print()` debug statements; all logs go through `logger`.
- ✅ Branch is rebased on latest `develop` before PR.
- ✅ PR description lists the files changed and the smoke-test steps performed.

### Per-Phase Handover Criteria

**P0 → P1 handover:**
- New: `utils/sanity.py`, populated `config/plot_configs.py`, extended `utils/plot_manager.py`, `utils/helpers.py:safe_render()`, `requirements.txt` includes `kaleido`.
- Existing modules untouched but should still launch.
- **Smoke test**: Open SCP module, generate a plot, click PNG/SVG/HTML download — all three files open and render. Wrap one tab in `safe_render`, raise an exception inside, confirm the error card appears and the other tabs are unaffected.
- **Hand-off artifact**: P1 begins by importing `safe_render`, `PlotManager` (extended), and `config.plot_configs` — these are the contracts P1 depends on.

**P1 → P2 handover:**
- New: `utils/report_builder.py`, `templates/report.html.j2`, Report tab added to `app.py`.
- `PlotManager.render_plot_and_editor()` now exposes "Add to Report" button.
- `st.session_state.report` is the canonical queue.
- **Smoke test**: Add 2 SCP figures to report, open Report tab, edit notes, download interactive HTML and ZIP bundle. Both open correctly; HTML figures are interactive; ZIP contains PNG+SVG+HTML+parameters.json.
- **Hand-off artifact**: P2 onwards uses the now-stable `report_builder.add_figure()` API; subsequent module upgrades just wire it in.

**P2 → P3 handover:**
- Comparative module fully migrated to `PlotManager` (heatmap on `MplPlotManager`), `safe_render` wrapping, sanity checks on load, Enrichr retry/timeout, TF list moved to config.
- **Smoke test**: Upload comparative fixture → generate volcano, heatmap, enrichment dotplot → edit each → add all to report → download both report formats. Intentionally use bad input file and confirm error stays scoped to Comparative tab.
- **Hand-off artifact**: TF list now lives in `config/plot_configs.py`; P3 (Quant) imports from there instead of redefining.

**P3 → P4 handover:**
- Quant module migrated; correlation matrix tab implemented (was a stub at `quant_module.py:321`); Venn/UpSet/Dendrogram on `MplPlotManager`.
- **Smoke test**: Upload quant fixture → generate every plot tab including new correlation matrix → confirm CV-vs-intensity plot appears → add to report → export.
- **Hand-off artifact**: `MplPlotManager` is now battle-tested across two modules; P4 uses it freely for any matplotlib output.

**P4 → P5 handover:**
- Dilution module: configurable column names, R²/CV/LOD-LOQ tables exportable as CSV, deviation buckets configurable, sanity checks for non-monotonic concentration.
- **Smoke test**: Upload dilution fixture with non-default column names → use UI to map them → generate all plots including new LOD/LOQ panel → export CSVs and full report.
- **Hand-off artifact**: Column-mapping UI pattern (used to handle non-default `Concentration` etc.) becomes the template for QC's metadata-mapping flow in P5.

**P5 → P6 handover:**
- DIA QC and Targeted QC migrated to `PlotManager`; configurable regex / column mapping; Targeted now has a metadata-review step matching DIA; sentinel peptide presets saved.
- **Smoke test**: Upload non-standard DIA-NN report (different column casing or regex format) → use new mapping UI → all QC plots render → control-chart sigma bands are configurable → CSV export works for both QC types.
- **Hand-off artifact**: All five user-facing modules now use the same patterns. P6 only refines the SCP module on top of stable foundations.

**P6 → P7 handover:**
- SCP: persistent `PlotManager` edits, AnnData save/load, expression-overlay UMAP, methods captured in report.
- **Smoke test**: Run full SCP pipeline → save state → restart app → load state → resume from clustering tab → continue to DE → all parameters appear in final report.
- **Hand-off artifact**: All five modules + SCP now feature-complete; P7 only adds deployment, fixtures, tests, and welcome experience.

**P7 → Production:**
- `Dockerfile`, `docker-compose.yml`, `.streamlit/config.toml`, `.dockerignore`, `Makefile` added.
- `tests/` with at least one smoke test per visualizer; CI runs `pytest`.
- Welcome page rebuilt with quick-start cards and sample fixtures in `tests/fixtures/`.
- **Smoke test**: `docker compose up` from a clean clone → access `http://localhost:8501` → walk through every module with sample data → produce a report → download both formats. `pytest tests/` passes.
- **Tag**: `v2.0.0` on `develop`, then merge `develop → main` for the deploy-ready release.

### Rollback Strategy per Phase

If a phase reveals a fundamental flaw at smoke-test time:
1. Don't merge that phase's PR.
2. Either fix on the same branch (preferred) or close the PR and re-cut from `develop`.
3. `develop` is always the last known-good state.

Tags `v2.0.0-p0` through `v2.0.0-p7` give you exact restore points.

---

## Verification Plan

### Per-module smoke test

For each module after upgrade:
1. Launch app: `streamlit run app.py`.
2. Upload provided fixture data for that module.
3. Walk through every tab, generating each plot.
4. Edit one plot's title, height, color → verify edit persists after re-generate.
5. Click "Add to Report" on at least 2 plots.
6. Download PNG, SVG, HTML for one figure each — open and confirm they render.
7. Intentionally break input (delete a required column, submit empty file) → confirm error stays in tab and other tabs still work.

### Cross-module test

1. Run analyses in 3 modules in one session, add 5+ figures to report.
2. Open Report tab, reorder figures, edit notes, generate both outputs.
3. Download interactive HTML, open in fresh browser tab — every figure interactive, parameters table populated, no broken links.
4. Download ZIP bundle, unzip — verify `figures/` has PNG+SVG+HTML per figure, `parameters.json` is valid, `notes.md` matches what was entered.

### Docker test (P7)

1. `docker compose build && docker compose up`.
2. Hit `http://localhost:8501` → run a full smoke test inside the container.
3. Confirm `kaleido` static export works (Chromium deps installed correctly).
4. Confirm uploads up to ~500 MB succeed.

### Automated tests (Phase 7)

`pytest tests/ -q` — under 60 seconds, runs in CI.

---

## Risk Register

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| `kaleido` breaks on deploy (Chromium dep) | Medium | Test in clean Docker early in P0 with the full Chromium dep list (`libgbm1`, `libasound2`, `libxshmfence1`); fall back to client-side download via `plotly.js` if blocked |
| Persisting PlotManager edits introduces stale-state bugs | Medium | Tag edits with figure-data hash; invalidate when underlying data changes |
| HTML report file size large with many embedded Plotly figs | Low | Use `include_plotlyjs='cdn'`; offer "lite" mode that exports static PNG instead of interactive |
| Refactoring TF list / config breaks existing plots | Low | Add unit test asserting same TF set before/after move |
| Migrating matplotlib plots to MplPlotManager regresses look | Medium | Visual-diff fixture screenshots in test suite |
| Enrichr API instability | Already happens | Wrapping with retries + clear UX for failure (P4) |

---

## Critical Files Reference

Foundation:
- `utils/plot_manager.py` (extend)
- `utils/helpers.py` (extend with `safe_render`, `to_hex`)
- `utils/sanity.py` (new)
- `utils/report_builder.py` (new)
- `config/plot_configs.py` (populate from empty)
- `templates/report.html.j2` (new)
- `requirements.txt` (add kaleido, jinja2 pin)

Modules (in upgrade order):
- `modules/comparative_module.py`, `visualizations/comparative_visualizer.py`
- `modules/quant_module.py`, `visualizations/quant_visualizer.py`
- `modules/dilution_module.py`, `visualizations/dilution_series.py`
- `modules/qc_module.py`, `modules/qc_tabs/dia_qc_tab.py`, `modules/qc_tabs/targeted_qc_tab.py`
- `visualizations/DiaQcVisualizer.py`, `visualizations/targettedQCVisualization.py`
- `modules/scp_module.py`, `visualizations/scp_visualizer.py`

Reference patterns to follow (already in repo):
- Session state helpers: `modules/scp_module.py:37-51`
- Global settings + theme/color UI: `modules/scp_module.py:82-126`
- Tab structure pattern: `modules/scp_module.py:1057-1093`
- Per-stage try/except + logging: `modules/scp_module.py:177-199, 670-672, 855-867`
- PlotManager usage: `modules/scp_module.py:301-305`
- Decorator definition: `utils/helpers.py:8-21`
- PlotManager class: `utils/plot_manager.py:4-108`
