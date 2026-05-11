import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scanpy", reason="scanpy not installed")

from visualizations.scp_visualizer import SCPVisualizer


@pytest.fixture
def scp_viz():
    np.random.seed(99)
    n_proteins, n_cells = 60, 20
    samples = [f"cell_{i}" for i in range(n_cells)]
    proteins = [f"PROT_{i}" for i in range(n_proteins)]

    # PG matrix: proteins × cells (plus a Protein.Group column)
    matrix = np.random.lognormal(2, 0.8, (n_proteins, n_cells))
    # Introduce ~20 % missingness
    mask = np.random.random((n_proteins, n_cells)) < 0.2
    matrix[mask] = 0.0

    pg_df = pd.DataFrame(matrix, columns=samples)
    pg_df.insert(0, "Protein.Group", proteins)

    annotation_df = pd.DataFrame({
        "Run": samples,
        "condition": (["Ctrl"] * (n_cells // 2)) + (["Treat"] * (n_cells // 2)),
    })

    return SCPVisualizer(
        pg_matrix_df=pg_df,
        stats_df=None,
        annotation_df=annotation_df,
        sample_col="Run",
    )


def test_instantiation(scp_viz):
    assert scp_viz is not None
    assert scp_viz.adata.n_obs == 20
    assert scp_viz.adata.n_vars == 60


def test_compute_qc_metrics(scp_viz):
    scp_viz.compute_qc_metrics()
    assert "n_proteins" in scp_viz.adata.obs.columns
    assert "detection_rate" in scp_viz.adata.var.columns


def test_plot_protein_detection_histogram(scp_viz):
    scp_viz.compute_qc_metrics()
    fig = scp_viz.plot_protein_detection_histogram()
    assert fig is not None


def test_filter_and_preprocess(scp_viz):
    scp_viz.compute_qc_metrics()
    scp_viz.filter_samples(min_proteins=1)
    scp_viz.filter_proteins(min_detection_pct=5.0)
    scp_viz.preprocess()
    assert scp_viz.pp_state["log_transformed"]


def test_pca_and_plot(scp_viz):
    scp_viz.compute_qc_metrics()
    scp_viz.filter_samples(min_proteins=1)
    scp_viz.filter_proteins(min_detection_pct=5.0)
    scp_viz.preprocess()
    scp_viz.run_pca()
    fig = scp_viz.plot_pca(color_by="condition")
    assert fig is not None
