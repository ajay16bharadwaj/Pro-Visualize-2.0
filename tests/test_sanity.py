import numpy as np
import pandas as pd
import pytest

from utils.sanity import (
    ValidationResult,
    check_sample_alignment,
    check_value_ranges,
    detect_column,
    gene_resolution_report,
    has_blocking,
    summarize_missingness,
    validate_columns,
)


def test_validate_columns_all_present():
    df = pd.DataFrame({"A": [1], "B": [2], "C": [3]})
    result = validate_columns(df, required=["A", "B"])
    assert not has_blocking([result])


def test_validate_columns_missing_required():
    df = pd.DataFrame({"A": [1]})
    result = validate_columns(df, required=["A", "Missing"])
    assert has_blocking([result])


def test_validate_columns_missing_optional_is_warning():
    df = pd.DataFrame({"A": [1]})
    result = validate_columns(df, required=["A"], optional=["OptCol"])
    # optional missing → warning, not blocking
    assert not result.is_blocking
    assert result.severity == "warn"


def test_check_sample_alignment_all_match():
    protein_df = pd.DataFrame({"Protein": ["P1"], "S1": [1.0], "S2": [2.0]})
    annotation_df = pd.DataFrame({"Sample": ["S1", "S2"]})
    result = check_sample_alignment(protein_df, annotation_df, "Sample")
    assert not result.is_blocking


def test_check_sample_alignment_missing_samples():
    protein_df = pd.DataFrame({"Protein": ["P1"], "S1": [1.0]})
    annotation_df = pd.DataFrame({"Sample": ["S1", "S2", "S3"]})
    result = check_sample_alignment(protein_df, annotation_df, "Sample")
    assert isinstance(result, ValidationResult)


def test_check_value_ranges_no_issue():
    s = pd.Series([1.0, 2.0, 3.0])
    result = check_value_ranges(s, name="intensity")
    assert not result.is_blocking


def test_check_value_ranges_negative_values():
    s = pd.Series([-5.0, 2.0, 3.0])
    result = check_value_ranges(s, name="intensity", expected_non_negative=True)
    assert result.is_blocking


def test_summarize_missingness_full():
    df = pd.DataFrame({"A": [1.0, 2.0], "B": [3.0, 4.0]})
    summary = summarize_missingness(df)
    assert summary["missing_cells"] == 0
    assert summary["pct_missing"] == 0.0


def test_summarize_missingness_partial():
    df = pd.DataFrame({"A": [1.0, np.nan], "B": [np.nan, 4.0]})
    summary = summarize_missingness(df)
    assert summary["missing_cells"] == 2
    assert 0 < summary["pct_missing"] <= 100


def test_gene_resolution_report_all_resolved():
    df = pd.DataFrame({"Gene": ["GAPDH", "TP53", "EGFR"]})
    result = gene_resolution_report(df, gene_col="Gene")
    assert not result.is_blocking


def test_gene_resolution_report_partial():
    df = pd.DataFrame({"Gene": ["GAPDH", None, ""]})
    result = gene_resolution_report(df, gene_col="Gene")
    assert isinstance(result, ValidationResult)


def test_detect_column_finds_first_match():
    df = pd.DataFrame({"Gene names": [1], "Protein": [2]})
    found = detect_column(df, candidates=["Gene Name", "Gene names", "Gene"])
    assert found == "Gene names"


def test_detect_column_returns_none_when_absent():
    df = pd.DataFrame({"A": [1]})
    found = detect_column(df, candidates=["Gene Name", "Gene names"])
    assert found is None
