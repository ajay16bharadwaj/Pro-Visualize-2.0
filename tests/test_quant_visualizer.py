import numpy as np
import pandas as pd
import pytest

from visualizations.quant_visualizer import QuantificationVisualizer


@pytest.fixture
def viz():
    groups = {"GroupA": ["S1", "S2", "S3"], "GroupB": ["S4", "S5", "S6"]}
    samples = [s for g in groups.values() for s in g]

    annotation_df = pd.DataFrame({
        "Level3": samples,
        "attribute_ExperimentalGroup": [g for g, ss in groups.items() for _ in ss],
    })

    np.random.seed(0)
    n_proteins = 40
    proteins = [f"P{i}" for i in range(n_proteins)]
    data = {"Protein": proteins}
    for s in samples:
        vals = np.random.lognormal(2, 0.5, n_proteins)
        # introduce a few NaN to exercise missingness handling
        vals[np.random.choice(n_proteins, 5, replace=False)] = np.nan
        data[s] = vals
    protein_df = pd.DataFrame(data)

    return QuantificationVisualizer(protein_df, annotation_df)


def test_instantiation(viz):
    assert viz is not None
    assert len(viz.sample_cols) == 6


def test_plot_protein_counts(viz):
    fig = viz.plot_protein_counts()
    assert fig is not None


def test_plot_pca_by_annotation(viz):
    fig = viz.plot_pca_by_annotation(color_by="attribute_ExperimentalGroup")
    assert fig is not None


def test_plot_correlation_matrix_pearson(viz):
    fig = viz.plot_correlation_matrix(method="pearson")
    assert fig is not None


def test_plot_correlation_matrix_spearman(viz):
    fig = viz.plot_correlation_matrix(method="spearman")
    assert fig is not None


def test_plot_cv_vs_intensity(viz):
    fig = viz.plot_cv_vs_intensity()
    assert fig is not None


def test_plot_protein_rank_order(viz):
    fig = viz.plot_protein_rank_order()
    assert fig is not None


# --- Empty / degenerate-data guards (C3 robustness) --------------------------

@pytest.fixture
def empty_viz():
    """A visualizer whose intensity matrix is entirely missing."""
    samples = ["S1", "S2", "S3", "S4"]
    annotation_df = pd.DataFrame({
        "Level3": samples,
        "attribute_ExperimentalGroup": ["A", "A", "B", "B"],
    })
    proteins = [f"P{i}" for i in range(10)]
    data = {"Protein": proteins}
    for s in samples:
        data[s] = [np.nan] * len(proteins)
    return QuantificationVisualizer(pd.DataFrame(data), annotation_df)


def test_clustering_guard_on_empty_data(empty_viz):
    with pytest.raises(ValueError):
        empty_viz._prepare_data_for_clustering()


def test_upset_guard_on_empty_data(empty_viz):
    with pytest.raises(ValueError, match="no intersections|No proteins"):
        empty_viz.plot_upset()


def test_venn_guard_on_empty_data(empty_viz):
    with pytest.raises(ValueError, match="No proteins"):
        empty_viz.plot_venn_diagram(selected_groups=["A", "B"])
