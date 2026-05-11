import numpy as np
import pandas as pd
import pytest

from visualizations.comparative_visualizer import ComparativeVisualizer


@pytest.fixture
def comp_data():
    samples = [f"S{i}" for i in range(1, 9)]
    groups = ["Ctrl"] * 4 + ["Treat"] * 4

    annotation_df = pd.DataFrame({
        "Level3": samples,
        "attribute_ExperimentalGroup": groups,
    })

    np.random.seed(7)
    n_proteins = 50
    proteins = [f"P{i}" for i in range(n_proteins)]
    data = {"Protein": proteins, "Gene Name": proteins}
    for s in samples:
        data[s] = np.random.lognormal(2, 0.5, n_proteins)
    protein_df = pd.DataFrame(data)

    comparative_df = pd.DataFrame({
        "Protein": proteins,
        "log2FC": np.random.normal(0, 1.5, n_proteins),
        "Imputed.FDR": np.random.uniform(0, 1, n_proteins),
        "Comparison": ["Treat_vs_Ctrl"] * n_proteins,
    })

    column_config = {
        "protein_id": "Protein",
        "sample_id": "Level3",
        "grouping": "attribute_ExperimentalGroup",
        "comp_protein_id": "Protein",
        "fold_change": "log2FC",
        "fdr": "Imputed.FDR",
        "comparison_label": "Comparison",
    }

    return protein_df, annotation_df, comparative_df, column_config


@pytest.fixture
def viz(comp_data):
    protein_df, annotation_df, comparative_df, column_config = comp_data
    return ComparativeVisualizer(protein_df, annotation_df, comparative_df, column_config)


def test_instantiation(viz):
    assert viz is not None


def test_plot_volcano(viz):
    fig = viz.plot_volcano(fdr_cutoff=0.05, fc_cutoff=1.0)
    assert fig is not None


def test_plot_expression_violin(comp_data, viz):
    protein_df, *_ = comp_data
    proteins = protein_df["Protein"].tolist()[:5]
    fig = viz.plot_expression_violin(protein_list=proteins)
    assert fig is not None


def test_plot_comparative_heatmap(comp_data, viz):
    protein_df, *_ = comp_data
    proteins = protein_df["Protein"].tolist()[:10]
    result = viz.plot_comparative_heatmap(protein_list=proteins)
    # returns BytesIO
    assert result is not None
    assert hasattr(result, "read")
