"""
visualizations/scp_visualizer.py

Single-Cell Proteomics (SCP) Visualizer for Pro-Visualize.

Input format:
  - PG Matrix  : DIA-NN protein group matrix (wide TSV/CSV).
                 Rows = protein groups, Columns = sample run names + one protein ID column.
  - report.stats: DIA-NN per-run stats (File.Name, FWHM.RT, Normalisation.Instability, …).
  - Annotation  : CSV/TSV with a sample column matching PG-matrix columns + grouping columns.

Preprocessing pipeline (in order):
  compute_qc_metrics → filter_samples → filter_proteins →
  preprocess (normalize → log1p → regress → scale) →
  run_pca → [run_harmony] → run_neighbors → run_umap → run_leiden →
  run_de → compute_activity_scores
"""

import os
import logging
from io import BytesIO

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from scipy.sparse import issparse

try:
    import scanpy as sc
    import anndata as ad
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False

try:
    import harmonypy as hm
    HARMONY_AVAILABLE = True
except ImportError:
    HARMONY_AVAILABLE = False

try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Helper
# ─────────────────────────────────────────────────────────────────────────────

def _to_dense(X):
    """Return a dense numpy array regardless of sparse/dense input."""
    if issparse(X):
        return X.toarray()
    return np.asarray(X)


# ─────────────────────────────────────────────────────────────────────────────
# Main class
# ─────────────────────────────────────────────────────────────────────────────

class SCPVisualizer:
    """
    End-to-end single-cell proteomics analysis backed by AnnData / scanpy.

    Typical usage (matches the Streamlit module workflow):
        viz = SCPVisualizer(pg_df, stats_df, annotation_df, sample_col='Run')
        viz.compute_qc_metrics()
        viz.filter_samples(min_proteins=100)
        viz.filter_proteins(min_detection_pct=10)
        viz.preprocess(regress_covariates=['log_n_proteins','Normalisation.Instability'])
        viz.run_pca()
        viz.run_neighbors()
        viz.run_umap()
        viz.run_leiden(resolution=0.5)
        de = viz.run_de(groupby='condition', groups=['Treated'], reference='Control')
        scores = viz.compute_activity_scores({'OxPhos': ['...'], 'EMT': ['...']})
    """

    # DIA-NN report.stats columns we try to pull in
    DIANN_STATS_COLS = [
        "Proteins.Identified",
        "Precursors.Identified",
        "FWHM.RT",
        "Normalisation.Instability",
        "MS1.Signal",
    ]

    def __init__(
        self,
        pg_matrix_df: pd.DataFrame,
        stats_df: pd.DataFrame,
        annotation_df: pd.DataFrame,
        sample_col: str,
        run_col: str = "File.Name",
        protein_col: str = "Protein.Group",
    ):
        if not SCANPY_AVAILABLE:
            raise ImportError(
                "scanpy is required for SCP analysis. "
                "Install with:  pip install scanpy"
            )

        self.sample_col = sample_col
        self.run_col = run_col
        self.protein_col = protein_col

        self._pg_raw = pg_matrix_df.copy()
        self._stats_raw = stats_df.copy() if stats_df is not None else pd.DataFrame()
        self._annotation_raw = annotation_df.copy()

        # Track which preprocessing steps have been applied
        self.pp_state = {
            "qc_computed": False,
            "filtered": False,
            "normalized": False,
            "log_transformed": False,
            "regressed": False,
            "scaled": False,
            "pca_computed": False,
            "batch_corrected": False,
            "neighbors_computed": False,
            "umap_computed": False,
            "clustered": False,
            "de_computed": False,
        }

        self.adata = self._build_anndata()
        logger.info(
            f"SCPVisualizer ready: {self.adata.n_obs} cells × {self.adata.n_vars} proteins"
        )

    # ─────────────────────────────────────────────────────────────────────────
    # CONSTRUCTION
    # ─────────────────────────────────────────────────────────────────────────

    def _build_anndata(self) -> "sc.AnnData":
        pg_df = self._pg_raw.copy()

        # Identify protein column (first column if not found)
        if self.protein_col not in pg_df.columns:
            self.protein_col = pg_df.columns[0]
            logger.warning(f"Protein column not found; using first column: {self.protein_col}")

        pg_df = pg_df.set_index(self.protein_col)
        matrix_samples = pg_df.columns.tolist()
        anno_samples = self._annotation_raw[self.sample_col].astype(str).tolist()

        # Try direct match first, then basename match
        direct_match = [s for s in anno_samples if s in matrix_samples]
        if direct_match:
            used_samples = direct_match
            anno_index_col = self.sample_col
        else:
            def _stem(p):
                return os.path.splitext(os.path.basename(str(p).replace("\\", "/")))[0]

            base_map = {_stem(c): c for c in matrix_samples}
            matched_matrix_cols = []
            matched_anno_rows = []
            for s in anno_samples:
                s_base = _stem(s)
                if s_base in base_map:
                    matched_matrix_cols.append(base_map[s_base])
                    matched_anno_rows.append(s)
            if not matched_matrix_cols:
                raise ValueError(
                    "No annotation samples matched PG-matrix columns. "
                    "Verify that the 'Sample Column' values match column names in the PG matrix."
                )
            # Build a mapping from original annotation value → matrix column name
            self._annotation_raw["_matrix_col"] = self._annotation_raw[self.sample_col].apply(
                lambda s: base_map.get(_stem(s), np.nan)
            )
            used_samples = self._annotation_raw["_matrix_col"].dropna().tolist()
            anno_index_col = "_matrix_col"

        # Save non-sample metadata columns (e.g. Genes, Gene.Names) before subsetting
        meta_cols = [c for c in pg_df.columns if c not in set(used_samples)]
        var_meta = pg_df[meta_cols].copy() if meta_cols else pd.DataFrame(index=pg_df.index)

        # Subset PG matrix to matched samples only
        pg_df = pg_df[[s for s in used_samples if s in pg_df.columns]]

        # Build obs from annotation (indexed by matrix column name)
        obs_df = (
            self._annotation_raw.set_index(anno_index_col)
            .loc[pg_df.columns.tolist()]
            .copy()
        )
        obs_df.index = obs_df.index.astype(str)

        # X: cells × proteins  (NaN → 0 deferred to preprocess step)
        X = pg_df.T.values.astype(float)
        var_df = var_meta.reindex(pg_df.index)

        adata = sc.AnnData(X=X, obs=obs_df, var=var_df)
        adata.layers["input"] = X.copy()

        # Merge DIA-NN stats
        self._merge_stats(adata)
        return adata

    def _merge_stats(self, adata):
        """Merge report.stats QC metrics into adata.obs."""
        if self._stats_raw.empty or self.run_col not in self._stats_raw.columns:
            logger.warning("report.stats unavailable or missing run column — skipping stats merge.")
            return

        stats = self._stats_raw.copy()
        stats["_base"] = stats[self.run_col].apply(
            lambda x: os.path.splitext(os.path.basename(str(x).replace("\\", "/")))[0]
        )
        adata_bases = pd.Series(
            [os.path.splitext(os.path.basename(str(s).replace("\\", "/")))[0] for s in adata.obs.index],
            index=adata.obs.index,
        )
        stats_idx = stats.set_index("_base")
        available = [c for c in self.DIANN_STATS_COLS if c in stats_idx.columns]
        if not available:
            logger.warning("No standard DIA-NN stats columns found in report.stats.")
            return

        for col in available:
            adata.obs[col] = [
                stats_idx.loc[b, col] if b in stats_idx.index else np.nan
                for b in adata_bases.values
            ]
            adata.obs[col] = pd.to_numeric(adata.obs[col], errors="coerce")

        logger.info(f"Merged stats: {available}")

    # ─────────────────────────────────────────────────────────────────────────
    # QC & FILTERING
    # ─────────────────────────────────────────────────────────────────────────

    def compute_qc_metrics(self):
        """Compute per-cell and per-protein QC metrics from the input data."""
        X = _to_dense(self.adata.layers["input"]).astype(float)
        X_nan = np.where(X == 0, np.nan, X)

        # Per-cell
        self.adata.obs["n_proteins"] = np.sum(~np.isnan(X_nan), axis=1).astype(int)
        self.adata.obs["pct_detected"] = (
            self.adata.obs["n_proteins"] / self.adata.n_vars * 100
        )
        self.adata.obs["total_intensity"] = np.nansum(X_nan, axis=1)
        self.adata.obs["median_intensity"] = np.nanmedian(X_nan, axis=1)
        self.adata.obs["log_n_proteins"] = np.log1p(self.adata.obs["n_proteins"])

        # Per-protein
        self.adata.var["n_cells_detected"] = np.sum(~np.isnan(X_nan), axis=0).astype(int)
        self.adata.var["detection_rate"] = (
            self.adata.var["n_cells_detected"] / self.adata.n_obs * 100
        )
        self.adata.var["mean_intensity"] = np.nanmean(X_nan, axis=0)
        with np.errstate(invalid="ignore"):
            self.adata.var["cv_pct"] = (
                np.nanstd(X_nan, axis=0) / np.nanmean(X_nan, axis=0) * 100
            )

        self.pp_state["qc_computed"] = True
        logger.info(
            f"QC metrics computed. n_proteins: "
            f"{self.adata.obs['n_proteins'].min():.0f}–{self.adata.obs['n_proteins'].max():.0f}"
        )

    def filter_samples(
        self, min_proteins: int = 100, max_proteins: int = None
    ) -> int:
        """Remove samples with fewer than min_proteins detected proteins."""
        if not self.pp_state["qc_computed"]:
            self.compute_qc_metrics()
        n_before = self.adata.n_obs
        mask = self.adata.obs["n_proteins"] >= min_proteins
        if max_proteins is not None:
            mask &= self.adata.obs["n_proteins"] <= max_proteins
        self.adata = self.adata[mask].copy()
        removed = n_before - self.adata.n_obs
        logger.info(f"Sample filter: {n_before} → {self.adata.n_obs} (removed {removed})")
        self.pp_state["filtered"] = True
        return removed

    def filter_proteins(self, min_detection_pct: float = 10.0) -> int:
        """Remove proteins detected in fewer than min_detection_pct% of cells."""
        if not self.pp_state["qc_computed"]:
            self.compute_qc_metrics()
        n_before = self.adata.n_vars
        mask = self.adata.var["detection_rate"] >= min_detection_pct
        self.adata = self.adata[:, mask].copy()
        removed = n_before - self.adata.n_vars
        logger.info(f"Protein filter: {n_before} → {self.adata.n_vars} (removed {removed})")
        return removed

    def filter_by_qc_metric(
        self, metric: str, threshold: float, direction: str = "above"
    ) -> int:
        """Filter samples by a QC metric (direction: 'above' or 'below')."""
        if metric not in self.adata.obs.columns:
            raise ValueError(f"Metric '{metric}' not in obs columns.")
        n_before = self.adata.n_obs
        mask = (
            self.adata.obs[metric] >= threshold
            if direction == "above"
            else self.adata.obs[metric] <= threshold
        )
        self.adata = self.adata[mask].copy()
        return n_before - self.adata.n_obs

    # ─────────────────────────────────────────────────────────────────────────
    # PREPROCESSING
    # ─────────────────────────────────────────────────────────────────────────

    def preprocess(
        self,
        normalize: bool = True,
        log_transform: bool = True,
        regress_covariates: list = None,
        scale: bool = False,
    ):
        """
        Run the core preprocessing pipeline.
        Modifies adata.X in-place and stores each stage in a layer.
        """
        adata = self.adata

        # Replace NaN with 0 (undetected = 0 in SCP)
        X = _to_dense(adata.layers["input"]).astype(float)
        X = np.nan_to_num(X, nan=0.0)
        adata.X = X
        adata.layers["raw"] = X.copy()

        if normalize:
            sc.pp.normalize_total(adata, target_sum=None)
            adata.layers["normalized"] = _to_dense(adata.X).copy()
            self.pp_state["normalized"] = True
            logger.info("✓ Normalized (median library size)")

        if log_transform:
            sc.pp.log1p(adata)
            adata.layers["log1p"] = _to_dense(adata.X).copy()
            self.pp_state["log_transformed"] = True
            logger.info("✓ log1p transform applied")

        if regress_covariates:
            valid_covs = []
            for cov in regress_covariates:
                if cov not in adata.obs.columns:
                    logger.warning(f"Regression covariate '{cov}' missing — skipped")
                    continue
                if adata.obs[cov].isna().any():
                    adata.obs[cov] = adata.obs[cov].fillna(adata.obs[cov].median())
                    logger.warning(f"Filled NaN in '{cov}' with median before regression")
                valid_covs.append(cov)
            if valid_covs:
                adata.layers["before_regression"] = _to_dense(adata.X).copy()
                sc.pp.regress_out(adata, valid_covs)
                adata.layers["after_regression"] = _to_dense(adata.X).copy()
                self.pp_state["regressed"] = True
                logger.info(f"✓ Regressed out: {', '.join(valid_covs)}")

        if scale:
            sc.pp.scale(adata, max_value=None)
            adata.layers["scaled"] = _to_dense(adata.X).copy()
            self.pp_state["scaled"] = True
            logger.info("✓ Standard scaling applied")

        return self

    # ─────────────────────────────────────────────────────────────────────────
    # DIMENSIONALITY REDUCTION & CLUSTERING
    # ─────────────────────────────────────────────────────────────────────────

    def run_pca(self, n_comps: int = 50):
        """Run PCA. n_comps capped by data dimensions automatically."""
        n_comps = min(n_comps, self.adata.n_obs - 1, self.adata.n_vars - 1)
        sc.tl.pca(self.adata, n_comps=n_comps)
        self.pp_state["pca_computed"] = True
        logger.info(f"✓ PCA ({n_comps} components)")

    def run_harmony(self, batch_key: str):
        """Apply Harmony batch correction on the PCA embedding."""
        if not HARMONY_AVAILABLE:
            raise ImportError(
                "harmonypy is required for batch correction. "
                "Install with:  pip install harmonypy"
            )
        if not self.pp_state["pca_computed"]:
            self.run_pca()
        if batch_key not in self.adata.obs.columns:
            raise ValueError(f"Batch key '{batch_key}' not found in obs.")

        ho = hm.run_harmony(
            self.adata.obsm["X_pca"],
            self.adata.obs,
            batch_key,
            random_state=42,
        )
        self.adata.obsm["X_pca_harmony"] = ho.Z_corr.T
        self.adata.obsm["X_pca"] = self.adata.obsm["X_pca_harmony"].copy()
        self.pp_state["batch_corrected"] = True
        logger.info(f"✓ Harmony batch correction on '{batch_key}'")

    def run_neighbors(
        self, n_neighbors: int = 15, n_pcs: int = 30, use_harmony: bool = False
    ):
        """Compute the neighborhood graph (prerequisite for UMAP/Leiden)."""
        if not self.pp_state["pca_computed"]:
            self.run_pca()
        rep = (
            "X_pca_harmony"
            if (use_harmony and self.pp_state["batch_corrected"])
            else "X_pca"
        )
        sc.pp.neighbors(
            self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep
        )
        self.pp_state["neighbors_computed"] = True
        logger.info(f"✓ Neighbors (k={n_neighbors}, n_pcs={n_pcs})")

    def run_umap(self, min_dist: float = 0.3, spread: float = 1.0):
        """Compute UMAP embedding."""
        if not self.pp_state["neighbors_computed"]:
            self.run_neighbors()
        sc.tl.umap(self.adata, min_dist=min_dist, spread=spread, random_state=42)
        self.pp_state["umap_computed"] = True
        logger.info("✓ UMAP computed")

    def run_leiden(self, resolution: float = 0.5) -> int:
        """Run Leiden clustering. Returns number of clusters found."""
        if not self.pp_state["neighbors_computed"]:
            self.run_neighbors()
        sc.tl.leiden(self.adata, resolution=resolution, random_state=42)
        n_clusters = self.adata.obs["leiden"].nunique()
        self.pp_state["clustered"] = True
        logger.info(f"✓ Leiden (res={resolution}) → {n_clusters} clusters")
        return n_clusters

    # ─────────────────────────────────────────────────────────────────────────
    # DIFFERENTIAL EXPRESSION
    # ─────────────────────────────────────────────────────────────────────────

    def run_de(
        self,
        groupby: str,
        groups: list = None,
        reference: str = "rest",
        method: str = "wilcoxon",
        min_cells: int = 3,
        min_detection_pct: float = 0.25,
    ) -> dict:
        """
        Run differential expression using Wilcoxon rank-sum (or t-test).

        Returns dict {group_name: DataFrame} with columns:
            protein, log2FC, pval, pval_adj, score, pct_group, pct_ref,
            protein_class (quantitative / switch_gained / switch_lost / low_detect)
        """
        if groupby not in self.adata.obs.columns:
            raise ValueError(f"Column '{groupby}' not found in obs.")

        # Use log1p layer for DE
        adata_de = self.adata.copy()
        if "log1p" in adata_de.layers:
            adata_de.X = adata_de.layers["log1p"].copy()

        all_groups = adata_de.obs[groupby].astype(str).unique().tolist()
        groups_to_test = [str(g) for g in groups] if groups else all_groups
        de_results = {}

        for group in groups_to_test:
            n_g = (adata_de.obs[groupby].astype(str) == group).sum()
            if reference == "rest":
                n_r = (adata_de.obs[groupby].astype(str) != group).sum()
            else:
                n_r = (adata_de.obs[groupby].astype(str) == str(reference)).sum()

            if n_g < min_cells or n_r < min_cells:
                logger.warning(f"Skipping '{group}': n={n_g} vs n_ref={n_r} (min={min_cells})")
                continue

            # Subset to relevant cells
            if reference != "rest":
                mask = adata_de.obs[groupby].astype(str).isin([group, str(reference)])
                adata_sub = adata_de[mask].copy()
            else:
                adata_sub = adata_de.copy()

            # Filter proteins by detection
            X_sub = _to_dense(adata_sub.X)
            det_rate = (X_sub != 0).mean(axis=0)
            prot_mask = det_rate >= min_detection_pct
            adata_sub = adata_sub[:, prot_mask].copy()

            sc.tl.rank_genes_groups(
                adata_sub,
                groupby=groupby,
                groups=[group],
                reference=str(reference) if reference != "rest" else "rest",
                method=method,
                pts=True,
                tie_correct=True,
                use_raw=False,
            )

            result = sc.get.rank_genes_groups_df(adata_sub, group=group)
            result = result.rename(
                columns={
                    "names": "protein",
                    "logfoldchanges": "log2FC",
                    "pvals": "pval",
                    "pvals_adj": "pval_adj",
                    "scores": "score",
                }
            )
            result["group"] = group
            result["reference"] = reference

            # Detection fractions
            X_full = _to_dense(adata_sub.X)
            grp_idx = (adata_sub.obs[groupby].astype(str) == group).values
            ref_idx = ~grp_idx if reference == "rest" else (
                adata_sub.obs[groupby].astype(str) == str(reference)
            ).values

            gene_names = adata_sub.var_names.tolist()
            pct_g = pd.Series(
                (X_full[grp_idx] != 0).sum(axis=0) / grp_idx.sum(), index=gene_names
            )
            pct_r = pd.Series(
                (X_full[ref_idx] != 0).sum(axis=0) / ref_idx.sum(), index=gene_names
            )
            result["pct_group"] = result["protein"].map(pct_g).values
            result["pct_ref"] = result["protein"].map(pct_r).values

            # Classify proteins (quantitative vs switch-like)
            result = self._classify_de_proteins(result)

            de_results[group] = result
            n_sig = (result["pval_adj"] < 0.05).sum()
            logger.info(
                f"DE '{group}' vs '{reference}': "
                f"{n_sig} significant (adj p<0.05)"
            )

        self.adata.uns["de_results"] = de_results
        self.pp_state["de_computed"] = True
        return de_results

    @staticmethod
    def _classify_de_proteins(
        df: pd.DataFrame,
        high: float = 0.25,
        low: float = 0.05,
    ) -> pd.DataFrame:
        """Classify DE proteins as quantitative, switch-gained, switch-lost, or low-detect."""
        df = df.copy()
        gained = (df["pct_group"] >= high) & (df["pct_ref"] < low)
        lost = (df["pct_ref"] >= high) & (df["pct_group"] < low)
        quant = (df["pct_group"] >= high) & (df["pct_ref"] >= high)
        df["protein_class"] = "low_detect"
        df.loc[gained, "protein_class"] = "switch_gained"
        df.loc[lost, "protein_class"] = "switch_lost"
        df.loc[quant, "protein_class"] = "quantitative"
        return df

    def get_de_results(self, group: str = None) -> pd.DataFrame:
        """Return DE results for one group or all groups concatenated."""
        if "de_results" not in self.adata.uns:
            raise ValueError("Run run_de() first.")
        de = self.adata.uns["de_results"]
        if group:
            return de.get(group, pd.DataFrame())
        return pd.concat(de.values(), ignore_index=True) if de else pd.DataFrame()

    def get_de_groups(self) -> list:
        """Return list of groups for which DE results exist."""
        return list(self.adata.uns.get("de_results", {}).keys())

    # ─────────────────────────────────────────────────────────────────────────
    # ACTIVITY SCORING
    # ─────────────────────────────────────────────────────────────────────────

    def _build_gene_symbol_map(self) -> dict:
        """Return {gene_symbol: var_name} from adata.var gene-name columns.

        Used so that enrichment results (gene symbols) can be resolved to
        var_names (which may be UniProt accessions) before scoring.
        """
        _GENE_COL_CANDIDATES = ["Genes", "Gene.Names", "Gene names", "Gene", "gene_names"]
        gene_col = next((c for c in _GENE_COL_CANDIDATES if c in self.adata.var.columns), None)
        if not gene_col:
            return {}
        mapping = {}
        for var_name, gene_sym in self.adata.var[gene_col].items():
            for sym in str(gene_sym).split(";"):
                sym = sym.strip()
                if sym.lower() not in ("", "nan", "none", "na"):
                    mapping[sym] = var_name
        return mapping

    def compute_activity_scores(self, gene_sets: dict) -> list:
        """Compute per-cell activity scores using scanpy's score_genes.

        gene_sets : dict  {score_name: [gene1, gene2, ...]}

        Genes can be var_names (protein group IDs) OR gene symbols — both are
        resolved via a symbol→var_name lookup so enrichment pathway genes
        (always gene symbols) work even when var_names are UniProt accessions.

        Returns list of successfully computed score column names.
        """
        if "log1p" in self.adata.layers:
            self.adata.X = self.adata.layers["log1p"].copy()

        available = set(self.adata.var_names)
        sym_to_var = self._build_gene_symbol_map()
        computed = []

        for name, genes in gene_sets.items():
            resolved = []
            for g in genes:
                if g in available:
                    resolved.append(g)
                elif g in sym_to_var:
                    resolved.append(sym_to_var[g])
            resolved = list(dict.fromkeys(resolved))  # deduplicate, preserve order

            if len(resolved) < 3:
                logger.warning(
                    f"Gene set '{name}': only {len(resolved)} genes matched in data (need ≥3). "
                    f"Tried direct var_name match and gene-symbol lookup."
                )
                continue
            score_col = f"{name}_score"
            sc.tl.score_genes(
                self.adata, gene_list=resolved, score_name=score_col, use_raw=False
            )
            logger.info(
                f"Scored '{name}' ({len(resolved)} genes): "
                f"mean={self.adata.obs[score_col].mean():.3f}"
            )
            computed.append(score_col)

        return computed

    def run_gsea_enrichment(
        self,
        de_group: str,
        pval_thresh: float = 0.05,
        fc_thresh: float = 1.0,
        gene_sets: list = None,
        direction: str = "both",
    ) -> "pd.DataFrame":
        import gseapy as gp

        if gene_sets is None:
            gene_sets = ["KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023"]

        de_df = self.get_de_results(de_group)
        if de_df.empty:
            raise ValueError(f"No DE results for group '{de_group}'.")

        sig = de_df[(de_df["pval_adj"] < pval_thresh) & (de_df["log2FC"].abs() > fc_thresh)]
        if direction == "up":
            sig = sig[sig["log2FC"] > 0]
        elif direction == "down":
            sig = sig[sig["log2FC"] < 0]

        dep_list = sig["protein"].dropna().tolist()
        if len(dep_list) < 3:
            raise ValueError(
                f"Only {len(dep_list)} DEPs pass the current filters (need ≥3). "
                "Try relaxing the p-value or log₂FC thresholds."
            )

        # Resolve gene symbols for Enrichr.
        # Priority 1: gene-name column already in adata.var (preserved from PG matrix)
        # Priority 2: fetch from UniProt for just the significant DEPs (lazy / on demand)
        # Priority 3: split the protein ID on ";" as a last resort
        _GENE_COL_CANDIDATES = ["Genes", "Gene.Names", "Gene names", "Gene", "gene_names"]
        gene_col = next(
            (c for c in _GENE_COL_CANDIDATES if c in self.adata.var.columns), None
        )

        if gene_col is None:
            logger.info("No gene-name column in adata.var — querying UniProt for significant DEPs.")
            try:
                acc_map = {p: self._parse_uniprot_accession(p) for p in dep_list}
                fetched = self._fetch_uniprot_gene_names(list(dict.fromkeys(acc_map.values())))
                if fetched:
                    if "Genes" not in self.adata.var.columns:
                        self.adata.var["Genes"] = ""
                    for protein, acc in acc_map.items():
                        if acc in fetched and protein in self.adata.var.index:
                            self.adata.var.at[protein, "Genes"] = fetched[acc]
                    gene_col = "Genes"
                    logger.info(f"UniProt resolved {len(fetched)} accessions to gene symbols.")
            except Exception as exc:
                logger.warning(f"UniProt fallback failed: {exc}")

        gene_list: list[str] = []
        for protein in dep_list:
            if gene_col and protein in self.adata.var.index:
                raw = str(self.adata.var.at[protein, gene_col])
            else:
                raw = protein
            for sym in raw.split(";"):
                sym = sym.strip()
                if sym and sym.lower() != "nan":
                    gene_list.append(sym)
        gene_list = list(dict.fromkeys(gene_list))  # deduplicate, preserve order
        logger.info(
            f"Enrichment input: {len(dep_list)} protein groups → {len(gene_list)} gene symbols"
            + (f" (via adata.var['{gene_col}'])" if gene_col else " (split from protein ID)")
        )

        if len(gene_list) < 3:
            raise ValueError(
                f"Only {len(gene_list)} unique gene symbols after parsing (need ≥3). "
                + (
                    f"Check that the '{gene_col}' column contains valid gene symbols."
                    if gene_col
                    else "Consider using 'Fetch Gene Names from UniProt' in the Upload tab."
                )
            )

        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism="human",
            outdir=None,
            verbose=False,
        )
        results = enr.results
        if results.empty:
            return pd.DataFrame()
        # Normalise column names across gseapy versions
        results.columns = [c.strip() for c in results.columns]
        adj_col = next(
            (c for c in results.columns if "adjusted" in c.lower() and "p" in c.lower()), None
        )
        if adj_col and adj_col != "Adjusted P-value":
            results = results.rename(columns={adj_col: "Adjusted P-value"})
        return results.sort_values("Adjusted P-value")

    def get_available_score_cols(self) -> list:
        return [c for c in self.adata.obs.columns if c.endswith("_score")]

    # ─────────────────────────────────────────────────────────────────────────
    # UNIPROT GENE NAME ANNOTATION
    # ─────────────────────────────────────────────────────────────────────────

    @staticmethod
    def _parse_uniprot_accession(protein_id: str) -> str:
        """Extract the primary UniProt accession from a protein group ID.

        Handles:
        - Semicolon-separated groups:  'P12345;Q67890'              → 'P12345'
        - Isoform suffixes:            'P12345-2'                   → 'P12345'
        - FASTA-style IDs:             'sp|P12345|GENE_HUMAN'       → 'P12345'
        - Contaminant prefixed IDs:    'contam_sp|O43790|KRT86_HUMAN' → 'O43790'
        """
        primary = str(protein_id).split(";")[0].strip()
        # FASTA-style: [prefix]|ACCESSION|NAME  — accession is always the middle field
        if "|" in primary:
            parts = primary.split("|")
            if len(parts) >= 2:
                primary = parts[1]
        return primary.split("-")[0].strip()

    @staticmethod
    def _fetch_uniprot_gene_names(
        accessions: list, batch_size: int = 50
    ) -> dict:
        """Batch-query the UniProt REST API (GET, ≤50 accessions per request).

        UniProt enforces a ~75-accession URL length limit; 50 keeps us safely under it.
        Returns {accession: first_gene_symbol}.  Failed batches are skipped so
        partial results are always returned.
        """
        import requests

        result = {}
        unique = list(dict.fromkeys(a for a in accessions if a))
        n_batches = (len(unique) + batch_size - 1) // batch_size
        for i in range(0, len(unique), batch_size):
            batch = unique[i : i + batch_size]
            query = " OR ".join(f"accession:{acc}" for acc in batch)
            batch_num = i // batch_size + 1
            try:
                resp = requests.get(
                    "https://rest.uniprot.org/uniprotkb/search",
                    params={
                        "query": query,
                        "fields": "accession,gene_names",
                        "format": "tsv",
                        "size": batch_size,
                    },
                    timeout=30,
                )
                resp.raise_for_status()
                for line in resp.text.strip().split("\n")[1:]:
                    parts = line.split("\t")
                    if len(parts) >= 2 and parts[1].strip():
                        result[parts[0].strip()] = parts[1].strip().split(" ")[0]
                logger.info(
                    f"UniProt batch {batch_num}/{n_batches}: "
                    f"{len(batch)} queried, {len(result)} total resolved."
                )
            except Exception as exc:
                logger.warning(f"UniProt batch {batch_num}/{n_batches} failed: {exc}")
        return result

    def annotate_gene_names_from_uniprot(self) -> int:
        """Query UniProt for gene symbols for all proteins in adata.var.

        Stores results in adata.var['Genes'].
        Returns the number of proteins that received a gene symbol.
        """
        var_names = self.adata.var_names.tolist()
        acc_map = {vn: self._parse_uniprot_accession(vn) for vn in var_names}
        gene_map = self._fetch_uniprot_gene_names(list(dict.fromkeys(acc_map.values())))

        self.adata.var["Genes"] = [
            gene_map.get(acc_map[vn], "") for vn in var_names
        ]
        n = int((self.adata.var["Genes"] != "").sum())
        logger.info(
            f"UniProt annotation: {n}/{len(var_names)} proteins mapped to gene symbols"
        )
        return n

    # ─────────────────────────────────────────────────────────────────────────
    # COVARIATE CORRELATION ANALYSIS
    # ─────────────────────────────────────────────────────────────────────────

    def get_covariate_pc_correlations(self, n_pcs: int = 5) -> pd.DataFrame:
        """Return DataFrame of correlations between obs numeric vars and top PCs."""
        if not self.pp_state["pca_computed"]:
            return pd.DataFrame()
        pca_mat = self.adata.obsm["X_pca"][:, :n_pcs]
        pc_labels = [f"PC{i+1}" for i in range(n_pcs)]
        numeric_cols = self.adata.obs.select_dtypes(include=[np.number]).columns.tolist()
        rows = []
        for col in numeric_cols:
            vals = self.adata.obs[col].fillna(0).values
            row = {"Variable": col}
            for i, pc in enumerate(pc_labels):
                row[pc] = (
                    float(np.corrcoef(vals, pca_mat[:, i])[0, 1])
                    if np.std(vals) > 0
                    else 0.0
                )
            rows.append(row)
        return pd.DataFrame(rows).set_index("Variable").round(4)

    # ─────────────────────────────────────────────────────────────────────────
    # QC PLOTS
    # ─────────────────────────────────────────────────────────────────────────

    def plot_qc_overview(self, color_by: str = None, **kwargs) -> go.Figure:
        """Violin panel: n_proteins, pct_detected, total_intensity."""
        if not self.pp_state["qc_computed"]:
            self.compute_qc_metrics()

        obs = self.adata.obs.reset_index()
        template = kwargs.get("template", "plotly_white")
        color_map = kwargs.get("color_discrete_map")
        metrics = ["n_proteins", "pct_detected", "total_intensity"]
        titles = ["Proteins Detected", "Detection Rate (%)", "Total Intensity (log scale)"]
        fig = make_subplots(rows=1, cols=3, subplot_titles=titles)
        safe_colors = px.colors.qualitative.Safe
        group_col = color_by

        if group_col and group_col in obs.columns:
            groups = obs[group_col].unique().tolist()
            for mi, metric in enumerate(metrics):
                for gi, grp in enumerate(groups):
                    vals = obs[obs[group_col] == grp][metric].dropna()
                    fill_color = (
                        color_map.get(grp) if color_map else safe_colors[gi % len(safe_colors)]
                    )
                    fig.add_trace(
                        go.Violin(
                            y=vals, name=grp, x0=grp,
                            box_visible=True, meanline_visible=True,
                            fillcolor=fill_color, line_color="black",
                            showlegend=(mi == 0), legendgroup=grp,
                        ),
                        row=1, col=mi + 1,
                    )
        else:
            for mi, metric in enumerate(metrics):
                fig.add_trace(
                    go.Violin(
                        y=obs[metric].dropna(), name=metric,
                        box_visible=True, meanline_visible=True,
                        fillcolor=safe_colors[0], line_color="black",
                        showlegend=False,
                    ),
                    row=1, col=mi + 1,
                )

        if "total_intensity" in metrics:
            fig.update_yaxes(type="log", row=1, col=3)

        fig.update_layout(template=template, height=500, title="Single-Cell QC Overview",
                          violinmode="overlay")
        return fig

    def plot_protein_detection_histogram(self, **kwargs) -> go.Figure:
        """Histogram of per-protein detection rates."""
        if not self.pp_state["qc_computed"]:
            self.compute_qc_metrics()
        fig = px.histogram(
            self.adata.var,
            x="detection_rate",
            nbins=50,
            title="Protein Detection Rate Distribution",
            labels={"detection_rate": "Detection Rate (%)", "count": "# Proteins"},
            template=kwargs.get("template", "plotly_white"),
        )
        fig.update_layout(height=400)
        return fig

    def plot_n_proteins_scatter(
        self, x_col: str = "n_proteins", y_col: str = "total_intensity",
        color_col: str = None, **kwargs
    ) -> go.Figure:
        """Scatter plot of two per-cell QC metrics."""
        if not self.pp_state["qc_computed"]:
            self.compute_qc_metrics()
        obs = self.adata.obs.reset_index()
        fig = px.scatter(
            obs, x=x_col, y=y_col, color=color_col,
            opacity=0.7,
            title=f"{y_col} vs {x_col}",
            template=kwargs.get("template", "plotly_white"),
            color_discrete_map=kwargs.get("color_discrete_map"),
        )
        fig.update_traces(marker=dict(size=7))
        fig.update_layout(height=500)
        return fig

    def plot_stats_metric(self, metric: str = "FWHM.RT", groupby: str = None, **kwargs) -> go.Figure:
        """Violin of a DIA-NN stats QC metric per group."""
        if metric not in self.adata.obs.columns:
            raise ValueError(f"Metric '{metric}' not found. Available: {list(self.adata.obs.columns)}")
        obs = self.adata.obs.reset_index()
        fig = px.violin(
            obs, y=metric, color=groupby, box=True, points="all",
            title=f"{metric} per Cell",
            template=kwargs.get("template", "plotly_white"),
            color_discrete_map=kwargs.get("color_discrete_map"),
        )
        fig.update_layout(height=500)
        return fig

    def plot_library_size_comparison(self, **kwargs) -> go.Figure:
        """Compare library sizes before / after normalization."""
        if "raw" not in self.adata.layers:
            raise ValueError("Run preprocess() first.")
        raw_sizes = _to_dense(self.adata.layers["raw"]).sum(axis=1)
        fig = make_subplots(rows=1, cols=2,
                            subplot_titles=["Before Normalization", "After Normalization"])
        fig.add_trace(go.Histogram(x=raw_sizes, nbinsx=40, name="Raw",
                                   marker_color="coral", showlegend=False), row=1, col=1)
        if "normalized" in self.adata.layers:
            norm_sizes = _to_dense(self.adata.layers["normalized"]).sum(axis=1)
            fig.add_trace(go.Histogram(x=norm_sizes, nbinsx=40, name="Normalized",
                                       marker_color="mediumseagreen", showlegend=False), row=1, col=2)
        fig.update_layout(template=kwargs.get("template", "plotly_white"),
                          height=400, title="Library Size Distribution")
        return fig

    def plot_covariate_correlation_heatmap(self, **kwargs) -> go.Figure:
        """Heatmap of numeric obs variables vs top PCs."""
        if not self.pp_state["pca_computed"]:
            raise ValueError("Run run_pca() first.")
        corr_df = self.get_covariate_pc_correlations()
        if corr_df.empty:
            raise ValueError("No numeric variables for correlation analysis.")
        fig = px.imshow(
            corr_df,
            color_continuous_scale="RdBu_r",
            zmin=-1, zmax=1,
            text_auto=".2f",
            title="Technical Covariate Correlations with PCs",
            template=kwargs.get("template", "plotly_white"),
        )
        fig.update_layout(height=max(350, len(corr_df) * 35))
        return fig

    # ─────────────────────────────────────────────────────────────────────────
    # EMBEDDING PLOTS
    # ─────────────────────────────────────────────────────────────────────────

    def _embedding_df(self, kind: str = "umap") -> tuple:
        """Return (DataFrame with embedding coords + obs, [col_x, col_y])."""
        if kind == "umap":
            if "X_umap" not in self.adata.obsm:
                raise ValueError("Run run_umap() first.")
            coords = self.adata.obsm["X_umap"]
            col_names = ["UMAP 1", "UMAP 2"]
        else:  # pca
            if "X_pca" not in self.adata.obsm:
                raise ValueError("Run run_pca() first.")
            coords = self.adata.obsm["X_pca"][:, :2]
            vr = self.adata.uns.get("pca", {}).get("variance_ratio", [0, 0])
            col_names = [f"PC1 ({vr[0]*100:.1f}%)", f"PC2 ({vr[1]*100:.1f}%)"]
        df = pd.DataFrame(coords, columns=col_names, index=self.adata.obs.index)
        df = pd.concat([df, self.adata.obs], axis=1).reset_index(names="Cell")
        return df, col_names

    def plot_pca(self, color_by: str = None, **kwargs) -> go.Figure:
        """Interactive PCA scatter plot."""
        if not self.pp_state["pca_computed"]:
            raise ValueError("Run run_pca() first.")
        df, cols = self._embedding_df("pca")
        color_col = color_by if color_by and color_by in df.columns else None
        is_cont = color_col and pd.api.types.is_numeric_dtype(df[color_col]) if color_col else False
        extra = dict(color_continuous_scale="RdYlBu_r") if is_cont else dict(
            color_discrete_map=kwargs.get("color_discrete_map")
        )
        fig = px.scatter(
            df, x=cols[0], y=cols[1], color=color_col,
            hover_name="Cell", opacity=0.8,
            title=f"PCA — {color_col}" if color_col else "PCA",
            template=kwargs.get("template", "plotly_white"),
            **extra,
        )
        fig.update_traces(marker=dict(size=8))
        fig.update_layout(height=600)
        return fig

    def plot_umap(self, color_by: str = None, **kwargs) -> go.Figure:
        """Interactive UMAP scatter plot."""
        if not self.pp_state["umap_computed"]:
            raise ValueError("Run run_umap() first.")
        df, cols = self._embedding_df("umap")
        color_col = color_by if color_by and color_by in df.columns else None
        is_cont = color_col and pd.api.types.is_numeric_dtype(df[color_col]) if color_col else False
        extra = dict(color_continuous_scale="RdYlBu_r") if is_cont else dict(
            color_discrete_map=kwargs.get("color_discrete_map")
        )
        fig = px.scatter(
            df, x=cols[0], y=cols[1], color=color_col,
            hover_name="Cell", opacity=0.7,
            title=f"UMAP — {color_col}" if color_col else "UMAP",
            template=kwargs.get("template", "plotly_white"),
            **extra,
        )
        fig.update_traces(marker=dict(size=6))
        fig.update_layout(height=600)
        return fig

    def plot_elbow(self, **kwargs) -> go.Figure:
        """PCA scree / elbow plot."""
        if not self.pp_state["pca_computed"]:
            raise ValueError("Run run_pca() first.")
        vr = self.adata.uns["pca"]["variance_ratio"]
        cumulative = np.cumsum(vr) * 100
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig.add_trace(
            go.Bar(x=list(range(1, len(vr) + 1)), y=vr * 100,
                   name="% Variance", marker_color="steelblue"),
        )
        fig.add_trace(
            go.Scatter(x=list(range(1, len(cumulative) + 1)), y=cumulative,
                       name="Cumulative %", line=dict(color="crimson", width=2),
                       mode="lines+markers", marker=dict(size=4)),
            secondary_y=True,
        )
        fig.update_layout(
            title="PCA Variance Explained", height=400,
            xaxis_title="PC", template=kwargs.get("template", "plotly_white"),
        )
        fig.update_yaxes(title_text="% Variance per PC", secondary_y=False)
        fig.update_yaxes(title_text="Cumulative %", secondary_y=True)
        return fig

    def plot_cluster_composition(
        self, groupby: str = "leiden", splitby: str = None, **kwargs
    ) -> go.Figure:
        """Stacked bar: cluster composition split by a metadata variable."""
        obs = self.adata.obs.copy()
        if groupby not in obs.columns:
            raise ValueError(f"Column '{groupby}' not in obs.")
        colors = px.colors.qualitative.Safe
        if splitby and splitby in obs.columns:
            ct = pd.crosstab(obs[groupby], obs[splitby], normalize="index") * 100
            fig = go.Figure()
            for i, col in enumerate(ct.columns):
                fig.add_trace(go.Bar(
                    x=ct.index.astype(str), y=ct[col], name=str(col),
                    marker_color=colors[i % len(colors)]
                ))
            fig.update_layout(barmode="stack",
                              title=f"Cluster Composition by {splitby}",
                              yaxis_title="% Cells")
        else:
            ct = obs[groupby].value_counts().sort_index()
            fig = px.bar(x=ct.index.astype(str), y=ct.values,
                         title=f"Cell Count per {groupby}",
                         labels={"x": groupby, "y": "# Cells"},
                         template=kwargs.get("template", "plotly_white"))
        fig.update_layout(template=kwargs.get("template", "plotly_white"),
                          height=500, xaxis_title=groupby)
        return fig

    # ─────────────────────────────────────────────────────────────────────────
    # DIFFERENTIAL EXPRESSION PLOTS
    # ─────────────────────────────────────────────────────────────────────────

    def plot_volcano_sc(
        self,
        de_df: pd.DataFrame,
        title: str = "Volcano Plot",
        fc_thresh: float = 1.0,
        pval_thresh: float = 0.05,
        n_label: int = 10,
        **kwargs,
    ) -> go.Figure:
        """
        Interactive volcano plot for SCP DE results.
        Distinguishes quantitative vs switch-like proteins with different shapes.
        """
        df = de_df.copy()
        df["-log10(pval_adj)"] = -np.log10(df["pval_adj"].clip(lower=1e-50))
        fc_cap = min(df["log2FC"].abs().quantile(0.98) * 1.3, 20)
        df["log2FC_disp"] = df["log2FC"].clip(-fc_cap, fc_cap)

        # Significance + direction
        df["sig"] = "Not significant"
        df.loc[(df["pval_adj"] < pval_thresh) & (df["log2FC"] > fc_thresh), "sig"] = "Up-regulated"
        df.loc[(df["pval_adj"] < pval_thresh) & (df["log2FC"] < -fc_thresh), "sig"] = "Down-regulated"
        df.loc[(df["pval_adj"] < pval_thresh) & (df["log2FC"].abs() <= fc_thresh), "sig"] = "Significant (low FC)"

        color_map = {
            "Up-regulated": "#c0392b",
            "Down-regulated": "#2980b9",
            "Significant (low FC)": "#f5b041",
            "Not significant": "#aaaaaa",
        }

        # Switch-like marker override
        symbol_map = {}
        if "protein_class" in df.columns:
            symbol_map = df.apply(
                lambda r: "triangle-up" if r.get("protein_class") == "switch_gained"
                else ("triangle-down" if r.get("protein_class") == "switch_lost"
                      else "circle"),
                axis=1,
            )

        fig = px.scatter(
            df, x="log2FC_disp", y="-log10(pval_adj)",
            color="sig", color_discrete_map=color_map,
            hover_data={c: True for c in ["protein", "log2FC", "pval_adj",
                                           "pct_group", "pct_ref", "protein_class"]
                        if c in df.columns},
            title=title,
            template=kwargs.get("template", "plotly_white"),
        )
        if symbol_map is not None and len(symbol_map) > 0:
            for trace in fig.data:
                pass  # custom symbols via update_traces is complex; using colour distinction

        fig.add_hline(y=-np.log10(pval_thresh), line_dash="dash", line_color="grey", opacity=0.5)
        fig.add_vline(x=fc_thresh, line_dash="dash", line_color="grey", opacity=0.5)
        fig.add_vline(x=-fc_thresh, line_dash="dash", line_color="grey", opacity=0.5)

        # Label top hits
        sig_df = df[df["sig"] != "Not significant"].copy()
        if len(sig_df) > 0:
            sig_df["rank"] = sig_df["-log10(pval_adj)"] + sig_df["log2FC"].abs() * 0.3
            top = sig_df.nlargest(n_label, "rank")
            for _, row in top.iterrows():
                fig.add_annotation(
                    x=row["log2FC_disp"], y=row["-log10(pval_adj)"],
                    text=row["protein"], showarrow=True,
                    arrowhead=2, arrowsize=1, arrowwidth=1,
                    ax=20, ay=-30, font=dict(size=9, color="black"),
                )

        n_up = (df["sig"] == "Up-regulated").sum()
        n_dn = (df["sig"] == "Down-regulated").sum()
        fig.add_annotation(x=0.97, y=0.97, xref="paper", yref="paper",
                           text=f"↑ {n_up}", showarrow=False,
                           font=dict(color="#c0392b", size=13, family="Arial Black"))
        fig.add_annotation(x=0.03, y=0.97, xref="paper", yref="paper",
                           text=f"↓ {n_dn}", showarrow=False,
                           font=dict(color="#2980b9", size=13, family="Arial Black"))
        fig.update_layout(height=700, xaxis_title="log₂ Fold Change",
                          yaxis_title="−log₁₀(adj. p-value)")
        return fig

    def plot_de_heatmap(
        self,
        de_df: pd.DataFrame,
        groupby: str,
        n_top: int = 25,
        pval_thresh: float = 0.05,
        fc_thresh: float = 0.5,
        title: str = "",
        figsize: tuple = None,
        dpi: int = 150,
        **kwargs,
    ) -> BytesIO:
        """Static seaborn clustermap of top DE proteins."""
        sig = de_df[
            (de_df["pval_adj"] < pval_thresh) & (de_df["log2FC"].abs() > fc_thresh)
        ].copy()
        top_proteins = sig.nsmallest(n_top, "pval_adj")["protein"].tolist()
        if not top_proteins:
            raise ValueError("No significant proteins found for heatmap.")

        prot_mask = self.adata.var_names.isin(top_proteins)
        adata_sub = self.adata[:, prot_mask].copy()
        layer_key = "log1p" if "log1p" in adata_sub.layers else None
        X = _to_dense(adata_sub.layers[layer_key] if layer_key else adata_sub.X)

        heat_df = pd.DataFrame(X.T, index=adata_sub.var_names, columns=adata_sub.obs.index)
        col_colors = None
        if groupby in adata_sub.obs.columns:
            groups = adata_sub.obs[groupby].astype(str)
            uniq = groups.unique()
            palette = sns.color_palette("tab10", len(uniq))
            col_colors = groups.map(dict(zip(uniq, palette)))

        plot_figsize = figsize if figsize else (12, max(8, len(heat_df) * 0.22))
        g = sns.clustermap(
            heat_df, method="ward", cmap="RdBu_r", z_score=0, center=0,
            col_colors=col_colors.to_frame() if col_colors is not None else None,
            yticklabels=(len(heat_df) <= 50),
            figsize=plot_figsize,
        )
        heading = title if title else f"Top DE Proteins (n={len(top_proteins)})"
        g.fig.suptitle(heading, y=1.02)
        buf = BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=dpi)
        buf.seek(0)
        plt.close(g.fig)
        return buf

    # ─────────────────────────────────────────────────────────────────────────
    # ACTIVITY SCORE PLOTS
    # ─────────────────────────────────────────────────────────────────────────

    def plot_expression_umap(self, protein: str, **kwargs) -> go.Figure:
        """UMAP coloured by a single protein's log-normalised expression."""
        if not self.pp_state["umap_computed"]:
            raise ValueError("Run run_umap() first.")
        if protein not in self.adata.var_names:
            raise ValueError(f"Protein '{protein}' not found in the dataset.")

        layer_key = "log1p" if "log1p" in self.adata.layers else None
        X = _to_dense(self.adata.layers[layer_key] if layer_key else self.adata.X)
        prot_idx = list(self.adata.var_names).index(protein)
        expr_vals = X[:, prot_idx]

        df, cols = self._embedding_df("umap")
        df["_expr"] = expr_vals

        fig = px.scatter(
            df, x=cols[0], y=cols[1], color="_expr",
            hover_name="Cell", opacity=0.8,
            color_continuous_scale="Viridis",
            title=f"UMAP — {protein} expression",
            labels={"_expr": "log(norm. intensity)"},
            template=kwargs.get("template", "plotly_white"),
        )
        fig.update_traces(marker=dict(size=6))
        fig.update_layout(height=600)
        return fig

    def get_protein_names(self) -> list:
        """Return sorted list of protein names for selection UI."""
        return sorted(self.adata.var_names.tolist())

    def plot_activity_umap(self, score_col: str, **kwargs) -> go.Figure:
        """UMAP coloured by a per-cell activity score (continuous colourscale)."""
        if score_col not in self.adata.obs.columns:
            raise ValueError(f"Score '{score_col}' not computed. Run compute_activity_scores() first.")
        return self.plot_umap(color_by=score_col, **kwargs)

    def plot_activity_violin(
        self, score_col: str, groupby: str, splitby: str = None, **kwargs
    ) -> go.Figure:
        """Interactive violin of an activity score across groups."""
        if score_col not in self.adata.obs.columns:
            raise ValueError(f"Score '{score_col}' not computed.")
        if groupby not in self.adata.obs.columns:
            raise ValueError(f"'{groupby}' not in obs.")
        obs = self.adata.obs.reset_index(names="Cell")
        color_col = splitby if splitby and splitby in obs.columns else groupby
        fig = px.violin(
            obs, x=groupby, y=score_col, color=color_col,
            box=True, points="all",
            title=f"{score_col} by {groupby}",
            template=kwargs.get("template", "plotly_white"),
            color_discrete_map=kwargs.get("color_discrete_map"),
        )
        fig.update_traces(jitter=0.05, pointpos=0, marker=dict(size=3, opacity=0.4))
        fig.update_layout(height=500, xaxis_title=groupby, yaxis_title=score_col)
        return fig

    def plot_activity_scores_panel(
        self, score_cols: list, groupby: str
    ) -> BytesIO:
        """
        Static multi-panel figure: UMAP + violin for each activity score.
        Returns a BytesIO PNG.
        """
        n = len(score_cols)
        if n == 0:
            raise ValueError("No score columns provided.")

        has_umap = self.pp_state["umap_computed"]
        n_rows = 2 if has_umap else 1
        fig, axes = plt.subplots(
            n_rows, n, figsize=(5 * n, 5 * n_rows), squeeze=False
        )

        umap_coords = self.adata.obsm.get("X_umap") if has_umap else None
        obs = self.adata.obs.copy()
        groups = obs[groupby].values if groupby in obs.columns else None
        unique_g = list(pd.unique(groups)) if groups is not None else []
        palette = dict(zip(unique_g, plt.cm.tab10.colors[: len(unique_g)]))

        for j, score_col in enumerate(score_cols):
            score_vals = obs[score_col].values
            name = score_col.replace("_score", "")

            # --- Row 0: UMAP coloured by score ---
            if has_umap and umap_coords is not None:
                ax = axes[0, j]
                vmin, vmax = np.percentile(score_vals, [2, 98])
                sc0 = ax.scatter(
                    umap_coords[:, 0], umap_coords[:, 1],
                    c=score_vals, cmap="RdYlBu_r",
                    s=10, alpha=0.7, vmin=vmin, vmax=vmax,
                    linewidths=0, rasterized=True,
                )
                plt.colorbar(sc0, ax=ax, shrink=0.6, label="Score")
                ax.set_title(f"{name}\n(per-cell activity)", fontsize=11, fontweight="bold")
                ax.set_xlabel("UMAP 1", fontsize=9)
                ax.set_ylabel("UMAP 2", fontsize=9)
                ax.tick_params(labelsize=7)

            # --- Row 1 (or 0): Violin ---
            vio_row = 1 if has_umap else 0
            ax_v = axes[vio_row, j]

            if groups is not None and len(unique_g) > 0:
                data_by_g = [
                    obs.loc[obs[groupby] == g, score_col].dropna().values
                    for g in unique_g
                ]
                positions = list(range(len(unique_g)))
                valid = [(d, p) for d, p in zip(data_by_g, positions) if len(d) > 0]
                if valid:
                    vp = ax_v.violinplot(
                        [d for d, _ in valid],
                        positions=[p for _, p in valid],
                        widths=0.6, showmeans=True, showmedians=False, showextrema=False,
                    )
                    for k, body in enumerate(vp["bodies"]):
                        body.set_facecolor(palette.get(unique_g[k], "grey"))
                        body.set_alpha(0.55)
                    vp["cmeans"].set_color("black")
                ax_v.set_xticks(positions)
                ax_v.set_xticklabels(
                    [str(g) for g in unique_g], rotation=40, ha="right", fontsize=8
                )

            ax_v.set_ylabel(f"{name} Score", fontsize=10)
            ax_v.set_title(f"{name} by {groupby}", fontsize=11, fontweight="bold", loc="left")
            ax_v.axhline(0, color="grey", lw=0.5, ls="--", alpha=0.5)
            ax_v.spines[["top", "right"]].set_visible(False)

        plt.tight_layout()
        buf = BytesIO()
        plt.savefig(buf, format="png", bbox_inches="tight", dpi=150)
        buf.seek(0)
        plt.close(fig)
        return buf

    # ─────────────────────────────────────────────────────────────────────────
    # UTILITY / ACCESSORS
    # ─────────────────────────────────────────────────────────────────────────

    def get_obs_df(self) -> pd.DataFrame:
        return self.adata.obs.copy()

    def get_var_df(self) -> pd.DataFrame:
        return self.adata.var.copy()

    def get_available_groupby_cols(self) -> list:
        obs = self.adata.obs
        return [c for c in obs.columns if obs[c].dtype == object or obs[c].nunique() <= 30]

    def get_available_numeric_cols(self) -> list:
        return self.adata.obs.select_dtypes(include=[np.number]).columns.tolist()

    def get_preprocessing_summary(self) -> dict:
        return {
            "n_cells": self.adata.n_obs,
            "n_proteins": self.adata.n_vars,
            "layers": list(self.adata.layers.keys()),
            "embeddings": list(self.adata.obsm.keys()),
            "uns_keys": list(self.adata.uns.keys()),
            "pp_state": self.pp_state.copy(),
        }

    def get_available_diann_stats(self) -> list:
        """Return which DIA-NN stats columns are present in obs."""
        return [c for c in self.DIANN_STATS_COLS if c in self.adata.obs.columns]
