import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from visualizations.DiaQcVisualizer import DiaQcVisualizer


@pytest.fixture
def parquet_path(tmp_path):
    """Write a minimal DIA-NN-style parquet file."""
    n = 200
    runs = [f"20240101_ProjectX_L_{ng}ng_A01_{i}" for ng in [50, 100] for i in range(1, n // 4 + 1)]
    run_col = np.resize(runs, n)

    df = pd.DataFrame({
        "Run": run_col,
        "Precursor.Id": [f"PEPTIDE{i % 20}_1" for i in range(n)],
        "Precursor.Quantity": np.random.lognormal(2, 0.5, n),
        "Decoy": [0] * n,
        "Q.Value": np.random.uniform(0, 0.01, n),
        "Predicted.RT": np.random.uniform(10, 90, n),
        "RT": np.random.uniform(10, 90, n),
        "FWHM": np.random.uniform(0.1, 0.5, n),
        "IM": np.random.uniform(0.8, 1.2, n),
        "Predicted.IM": np.random.uniform(0.8, 1.2, n),
        "Mass.Evidence": np.random.normal(0, 5, n),
    })
    path = tmp_path / "report.parquet"
    df.to_parquet(path, index=False)
    return str(path)


def test_instantiation(parquet_path):
    viz = DiaQcVisualizer(parquet_path)
    assert viz is not None


def test_get_metadata(parquet_path):
    viz = DiaQcVisualizer(parquet_path)
    meta = viz.get_metadata()
    assert isinstance(meta, pd.DataFrame)
    assert "Run" in meta.columns


def test_get_processed_data(parquet_path):
    viz = DiaQcVisualizer(parquet_path)
    df = viz.get_processed_data()
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_missing_run_column_raises(tmp_path):
    bad = pd.DataFrame({"NotRun": ["x"], "Q.Value": [0.001]})
    path = tmp_path / "bad.parquet"
    bad.to_parquet(path, index=False)
    with pytest.raises(KeyError):
        DiaQcVisualizer(str(path))


def test_non_parquet_raises(tmp_path):
    path = tmp_path / "file.csv"
    path.write_text("Run,Q.Value\nx,0.001\n")
    with pytest.raises(ValueError):
        DiaQcVisualizer(str(path))
