import json
import zipfile
from io import BytesIO

import pandas as pd
import plotly.graph_objects as go
import pytest

from utils.report_builder import ReportBuilder


@pytest.fixture
def builder():
    return ReportBuilder()


def _make_fig():
    fig = go.Figure(go.Scatter(x=[1, 2, 3], y=[4, 5, 6]))
    fig.update_layout(title="Test Figure")
    return fig


def test_add_figure_and_len(builder):
    builder.add_figure("quant", "pca", _make_fig(), "PCA Plot", {})
    assert len(builder.items) == 1


def test_add_table(builder):
    df = pd.DataFrame({"A": [1, 2], "B": [3, 4]})
    builder.add_table("quant", "summary", df, "Summary Table")
    assert len(builder.items) == 1


def test_add_section(builder):
    builder.add_section("quant", "Methods", "We used PCA.")
    assert len(builder.items) == 1


def test_remove(builder):
    builder.add_figure("quant", "pca", _make_fig(), "PCA", {})
    builder.add_figure("quant", "heatmap", _make_fig(), "Heatmap", {})
    builder.remove("quant", "pca")
    keys = [item["key"] for item in builder.items]
    assert "pca" not in keys
    assert "heatmap" in keys


def test_reorder(builder):
    builder.add_figure("quant", "a", _make_fig(), "A", {})
    builder.add_figure("quant", "b", _make_fig(), "B", {})
    builder.add_figure("quant", "c", _make_fig(), "C", {})
    builder.reorder("quant", ["c", "a", "b"])
    quant_items = [i for i in builder.items if i["module"] == "quant"]
    assert [i["key"] for i in quant_items] == ["c", "a", "b"]


def test_export_html_contains_title(builder):
    builder.add_figure("quant", "pca", _make_fig(), "PCA Plot", {"n_components": 2})
    html_bytes = builder.export_html()
    assert isinstance(html_bytes, bytes)
    assert b"PCA Plot" in html_bytes


def test_export_html_valid_utf8(builder):
    builder.add_figure("comp", "volcano", _make_fig(), "Volcano", {})
    html_bytes = builder.export_html()
    html_bytes.decode("utf-8")  # raises if not valid


def test_export_zip_structure(builder):
    builder.add_figure("quant", "pca", _make_fig(), "PCA", {"param": 1})
    builder.add_section("quant", "Notes", "Analysis notes.")
    zip_bytes = builder.export_zip()
    assert isinstance(zip_bytes, bytes)

    with zipfile.ZipFile(BytesIO(zip_bytes)) as zf:
        names = zf.namelist()
        assert any("parameters.json" in n for n in names)
        assert any("figures/" in n for n in names)

    # parameters.json must be valid JSON
    with zipfile.ZipFile(BytesIO(zip_bytes)) as zf:
        params_name = next(n for n in zf.namelist() if "parameters.json" in n)
        params = json.loads(zf.read(params_name))
        assert isinstance(params, (dict, list))


def test_empty_builder_export(builder):
    html = builder.export_html()
    assert isinstance(html, bytes)
    z = builder.export_zip()
    assert isinstance(z, bytes)


# --- Reproducibility provenance (C4) -----------------------------------------

def test_provenance_has_versions(builder):
    prov = builder.provenance()
    assert prov["app_version"] == ReportBuilder.APP_VERSION
    assert "python_version" in prov
    # at least the core packages are resolvable in the test env
    assert "pandas" in prov["packages"]
    assert "streamlit" in prov["packages"]


def test_html_contains_methods_section(builder):
    builder.add_figure("quant", "pca", _make_fig(), "PCA Plot", {})
    html = builder.export_html().decode("utf-8")
    assert "Methods &amp; Reproducibility" in html or "Methods & Reproducibility" in html
    assert "pandas" in html  # a package version row is rendered


def test_zip_contains_provenance(builder):
    builder.add_figure("quant", "pca", _make_fig(), "PCA", {"param": 1})
    with zipfile.ZipFile(BytesIO(builder.export_zip())) as zf:
        assert "provenance.json" in zf.namelist()
        prov = json.loads(zf.read("provenance.json"))
        assert "packages" in prov and "python_version" in prov
        params = json.loads(zf.read("parameters.json"))
        assert "_provenance" in params and "figures" in params


def test_add_figure_drops_non_serializable_params(builder):
    # A circular reference cannot be serialized even with default=str, so it is
    # the realistic "unserializable" case the guard protects export from.
    circular: dict = {}
    circular["self"] = circular

    builder.add_figure("quant", "pca", _make_fig(), "PCA",
                       {"good": 2, "bad": circular})
    stored = builder.items[0]["params"]
    assert stored == {"good": 2}  # the bad param was dropped, not crashed on
    # ...and export still succeeds.
    assert isinstance(builder.export_zip(), bytes)


def test_add_figure_non_dict_params(builder):
    builder.add_figure("quant", "pca", _make_fig(), "PCA", params=["not", "a", "dict"])
    assert builder.items[0]["params"] == {}
