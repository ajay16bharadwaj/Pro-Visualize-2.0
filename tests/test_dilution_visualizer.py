import numpy as np
import pandas as pd
import pytest

from visualizations.dilution_series import DilutionSeriesVisualizer


@pytest.fixture
def viz():
    concentrations = [1, 2, 4, 8, 16]
    reps = ["R1", "R2"]
    samples = [f"S_{c}_{r}" for c in concentrations for r in reps]

    metadata = pd.DataFrame({
        "Sample": samples,
        "Concentration": [c for c in concentrations for _ in reps],
        "Replicate": [r for _ in concentrations for r in reps],
        "Group": ["GroupA"] * len(samples),
    })

    np.random.seed(42)
    n_proteins = 30
    proteins = [f"PROT_{i}" for i in range(n_proteins)]
    data = {"Protein.Group": proteins, "Genes": proteins}
    for s in samples:
        data[s] = np.random.lognormal(mean=2.0, sigma=0.5, size=n_proteins)
    protein_df = pd.DataFrame(data)

    return DilutionSeriesVisualizer(protein_df, metadata)


def test_instantiation(viz):
    assert viz is not None
    assert len(viz.sample_cols) == 10


def test_plot_intensity_distribution(viz):
    fig = viz.plot_intensity_distribution()
    assert fig is not None


def test_plot_protein_counts_per_sample(viz):
    fig = viz.plot_protein_counts_per_sample()
    assert fig is not None


def test_plot_cv_distribution(viz):
    fig = viz.plot_cv_distribution()
    assert fig is not None


def test_plot_pca(viz):
    fig = viz.plot_pca()
    assert fig is not None


def test_plot_r2_histogram(viz):
    fig = viz.plot_r2_histogram()
    assert fig is not None


def test_plot_lod_loq(viz):
    fig = viz.plot_lod_loq()
    assert fig is not None
