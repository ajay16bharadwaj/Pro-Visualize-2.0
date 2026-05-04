"""
Input-validation and post-load sanity-check helpers shared across modules.

Severity model (matches the "Data Quality" expander UI in modules):
    "ok"    → ✓  green check, informational
    "warn"  → ⚠  yellow warning, allow user to proceed
    "error" → ✗  red error, blocks downstream computation

Each helper returns a `ValidationResult` so the caller can render a
consistent UI block via `render_validation()`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Result dataclass + UI renderer
# ---------------------------------------------------------------------------


@dataclass
class ValidationResult:
    """One validation outcome.

    `severity` drives UI behaviour: "error" blocks the pipeline, "warn"
    surfaces a non-blocking notice, "ok" is informational. `details` is
    free-form context (e.g. list of missing columns) shown in an expander.
    """

    severity: str  # "ok" | "warn" | "error"
    message: str
    details: list[str] = field(default_factory=list)

    @property
    def is_blocking(self) -> bool:
        return self.severity == "error"

    @property
    def icon(self) -> str:
        return {"ok": "✓", "warn": "⚠", "error": "✗"}.get(self.severity, "•")


def has_blocking(results: Iterable[ValidationResult]) -> bool:
    """True if any result is severity=error — caller should refuse to run."""
    return any(r.is_blocking for r in results)


def render_validation(results: Sequence[ValidationResult], header: str = "Data Quality") -> bool:
    """Render a list of results inside a Streamlit expander.

    Returns True if the data is safe to use (no blocking errors), False
    otherwise. Streamlit is imported lazily so this module can be unit-tested
    without a Streamlit runtime.
    """
    import streamlit as st  # noqa: WPS433 — intentional lazy import

    if not results:
        return True

    blocking = has_blocking(results)
    expand_default = blocking or any(r.severity == "warn" for r in results)
    label = f"{header} — {sum(1 for r in results if r.severity == 'ok')}/{len(results)} OK"

    with st.expander(label, expanded=expand_default):
        for r in results:
            line = f"{r.icon} {r.message}"
            if r.severity == "error":
                st.error(line)
            elif r.severity == "warn":
                st.warning(line)
            else:
                st.success(line)
            if r.details:
                with st.container():
                    for d in r.details[:25]:
                        st.caption(f"  • {d}")
                    if len(r.details) > 25:
                        st.caption(f"  …and {len(r.details) - 25} more")

    return not blocking


# ---------------------------------------------------------------------------
# Column validation
# ---------------------------------------------------------------------------


def validate_columns(
    df: pd.DataFrame,
    required: Sequence[str],
    optional: Sequence[str] | None = None,
    *,
    label: str = "input file",
) -> ValidationResult:
    """Confirm `df` has every column in `required`.

    Optional columns are reported as a warning only when missing — useful for
    things like an annotation file's `Group` column being nice-to-have.
    """
    optional = optional or []
    missing_required = [c for c in required if c not in df.columns]
    missing_optional = [c for c in optional if c not in df.columns]

    if missing_required:
        return ValidationResult(
            severity="error",
            message=f"{label}: missing {len(missing_required)} required column(s)",
            details=[f"Missing: {c}" for c in missing_required]
            + ([f"Available: {', '.join(df.columns[:20])}"] if len(df.columns) <= 20 else []),
        )

    if missing_optional:
        return ValidationResult(
            severity="warn",
            message=f"{label}: {len(missing_optional)} optional column(s) missing",
            details=[f"Missing: {c}" for c in missing_optional],
        )

    return ValidationResult(
        severity="ok",
        message=f"{label}: all {len(required)} required columns present",
    )


def detect_column(df: pd.DataFrame, candidates: Sequence[str]) -> str | None:
    """First candidate header present in `df.columns`, or None.

    Case-sensitive — proteomics tools are inconsistent enough that exact
    matches are safer than fuzzy. For wizard-style auto-detection, callers
    should fall back to a manual selectbox if this returns None.
    """
    cols = set(df.columns)
    for c in candidates:
        if c in cols:
            return c
    return None


# ---------------------------------------------------------------------------
# Sample alignment between protein matrix and annotation file
# ---------------------------------------------------------------------------


def check_sample_alignment(
    protein_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    sample_col: str,
) -> ValidationResult:
    """Compare annotation sample IDs against protein-matrix column headers.

    Common failure mode: scientist exported protein matrix with sample IDs
    that have a `.raw` suffix or a path prefix, but the annotation file uses
    bare names. Surface mismatches up front.
    """
    if sample_col not in annotation_df.columns:
        return ValidationResult(
            severity="error",
            message=f"Annotation file missing sample-link column: {sample_col!r}",
        )

    annotation_samples = set(annotation_df[sample_col].astype(str))
    protein_columns = set(protein_df.columns.astype(str))

    matched = annotation_samples & protein_columns
    only_annotation = sorted(annotation_samples - protein_columns)
    only_protein = sorted(protein_columns - annotation_samples)

    if not matched:
        return ValidationResult(
            severity="error",
            message="No samples match between annotation file and protein matrix",
            details=[
                f"Annotation samples (first 10): {', '.join(sorted(annotation_samples)[:10])}",
                f"Protein columns (first 10): {', '.join(sorted(protein_columns)[:10])}",
            ],
        )

    if only_annotation:
        return ValidationResult(
            severity="warn",
            message=f"{len(matched)}/{len(annotation_samples)} annotation samples found in protein matrix",
            details=[f"In annotation but missing from protein matrix: {s}" for s in only_annotation],
        )

    return ValidationResult(
        severity="ok",
        message=f"All {len(matched)} annotation samples matched in protein matrix",
        details=([f"Extra protein columns (likely metadata): {', '.join(only_protein[:10])}"]
                 if only_protein else []),
    )


# ---------------------------------------------------------------------------
# Numeric range / distribution sanity
# ---------------------------------------------------------------------------


def check_value_ranges(
    series: pd.Series,
    *,
    name: str,
    expected_non_negative: bool = True,
    flag_all_zero: bool = True,
) -> ValidationResult:
    """Flag suspicious distributions in a numeric column.

    Default rules tuned for proteomics intensity data:
    - Negative intensities are usually a bug (log-transform downstream will fail).
    - All-zero columns indicate empty samples that will break PCA/clustering.
    """
    numeric = pd.to_numeric(series, errors="coerce")
    n_total = len(numeric)
    n_nan = int(numeric.isna().sum())
    n_neg = int((numeric < 0).sum()) if expected_non_negative else 0
    n_zero = int((numeric == 0).sum())

    if expected_non_negative and n_neg > 0:
        return ValidationResult(
            severity="error",
            message=f"{name}: {n_neg} negative values detected (expected non-negative)",
        )

    if flag_all_zero and n_total > 0 and (n_zero + n_nan) == n_total:
        return ValidationResult(
            severity="error",
            message=f"{name}: all values are zero or missing",
        )

    if n_nan / max(n_total, 1) > 0.5:
        return ValidationResult(
            severity="warn",
            message=f"{name}: {n_nan}/{n_total} values are missing ({n_nan / n_total:.0%})",
        )

    return ValidationResult(
        severity="ok",
        message=f"{name}: distribution looks healthy ({n_total - n_nan} non-null, {n_nan} NaN)",
    )


# ---------------------------------------------------------------------------
# Missingness summary — used by every module's data overview
# ---------------------------------------------------------------------------


def summarize_missingness(df: pd.DataFrame) -> dict:
    """Return overall + per-column missingness metrics."""
    if df.empty:
        return {"rows": 0, "cols": 0, "total_cells": 0, "missing_cells": 0, "pct_missing": 0.0,
                "per_column_pct": {}}

    total = int(df.size)
    missing = int(df.isna().sum().sum())
    per_col = (df.isna().mean() * 100).round(2).to_dict()
    return {
        "rows": int(df.shape[0]),
        "cols": int(df.shape[1]),
        "total_cells": total,
        "missing_cells": missing,
        "pct_missing": round(100 * missing / total, 2) if total else 0.0,
        "per_column_pct": per_col,
    }


# ---------------------------------------------------------------------------
# Gene-symbol resolution report — mirrors SCP's "X/Y proteins have a gene
# symbol" caption (scp_module.py:209-226).
# ---------------------------------------------------------------------------


def gene_resolution_report(df: pd.DataFrame, gene_col: str) -> ValidationResult:
    """Report what fraction of rows have a usable gene symbol.

    "Usable" = non-empty, not literal "nan"/"none"/"na". Matches the
    cleaning logic SCP uses when building gene-symbol → var_name maps.
    """
    if gene_col not in df.columns:
        return ValidationResult(
            severity="warn",
            message=f"No gene-symbol column ({gene_col!r}) found — TF highlighting and "
                    f"enrichment will fall back to protein IDs",
        )

    total = len(df)
    if total == 0:
        return ValidationResult(severity="warn", message="Empty DataFrame")

    series = df[gene_col].astype(str).str.strip()
    invalid = series.str.lower().isin({"", "nan", "none", "na"})
    resolved = int((~invalid).sum())
    pct = 100 * resolved / total

    if resolved == 0:
        return ValidationResult(
            severity="error",
            message=f"No proteins have a gene symbol in column {gene_col!r}",
        )

    if pct < 50:
        return ValidationResult(
            severity="warn",
            message=f"Only {resolved}/{total} ({pct:.0f}%) proteins have a gene symbol",
            details=[
                "Downstream pathway enrichment and TF highlighting may miss many proteins.",
                "Consider mapping protein IDs to gene symbols via UniProt before upload.",
            ],
        )

    return ValidationResult(
        severity="ok",
        message=f"{resolved}/{total} ({pct:.0f}%) proteins have a gene symbol",
    )
