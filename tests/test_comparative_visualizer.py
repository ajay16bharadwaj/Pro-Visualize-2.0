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


# --- Enrichment background toggle (network mocked) ---------------------------

class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload


def _enrichr_row(term, pval, adj):
    # [rank, term, p_value, z, combined, [genes], adj_p_value, ...]
    return [1, term, pval, 0.0, 0.0, ["GENEA", "GENEB"], adj]


class _FakeSession:
    """Routes Enrichr/Speedrichr calls to canned data. Classic GETs return a
    distinct p-value from the background POSTs so the two modes are
    distinguishable. Records the background list it was given."""

    def __init__(self):
        self.background_seen = None

    def _lib_from(self, url, data):
        if data and "backgroundType" in data:
            return data["backgroundType"]
        return url.split("backgroundType=")[-1]

    def post(self, url, files=None, data=None, timeout=None):
        if url.endswith("/addList"):
            return _FakeResponse({"userListId": 42})
        if url.endswith("/addbackground"):
            self.background_seen = data["background"]
            return _FakeResponse({"backgroundid": "bg-1"})
        if url.endswith("/backgroundenrich"):
            lib = self._lib_from(url, data)
            return _FakeResponse({lib: [_enrichr_row("Pathway X", 0.04, 0.20)]})
        raise AssertionError(f"unexpected POST {url}")

    def get(self, url, timeout=None):
        lib = self._lib_from(url, None)
        return _FakeResponse({lib: [_enrichr_row("Pathway X", 0.01, 0.05)]})


def test_enrichment_whole_genome_mode(viz, monkeypatch):
    session = _FakeSession()
    monkeypatch.setattr(viz, "_enrichr_session", lambda: session)
    df = viz.run_enrichment_analysis(["GENEA", "GENEB"], organism="human")
    assert not df.empty
    assert session.background_seen is None  # classic path: no background posted
    assert (df["p_value"] == 0.01).all()


def test_enrichment_detected_background_mode(comp_data, viz, monkeypatch):
    session = _FakeSession()
    monkeypatch.setattr(viz, "_enrichr_session", lambda: session)
    background = comp_data[0]["Gene Name"].tolist()
    df = viz.run_enrichment_analysis(["GENEA", "GENEB"], organism="human",
                                     background_genes=background)
    assert not df.empty
    # The custom background was actually posted to Speedrichr...
    assert session.background_seen is not None
    assert set(session.background_seen.split("\n")) == set(background)
    # ...and the background-corrected p-values differ from whole-genome mode.
    assert (df["p_value"] == 0.04).all()
