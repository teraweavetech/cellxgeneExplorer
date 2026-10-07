"""Core tests (no network, no R). Run with: pytest"""
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from cxg_explorer import agent as A
from cxg_explorer.example import make_example

PKG = Path(A.__file__).parent


@pytest.fixture(scope="session")
def example(tmp_path_factory):
    return make_example(tmp_path_factory.mktemp("data") / "example.h5ad", n_cells=2500, n_genes=120)


def test_bh_matches_definition():
    p = np.array([0.01, 0.04, 0.03, 0.2])
    q = A.bh(p)
    assert np.allclose(q, [0.04, 0.16 / 3, 0.16 / 3, 0.2])  # computed by hand


def test_read_count_rows_matches_anndata(example):
    import anndata as ad
    a = ad.read_h5ad(example)
    rows = np.sort(np.random.default_rng(0).choice(a.n_obs, 300, replace=False))
    with h5py.File(example) as f:
        m = A.read_count_rows(f, "raw/X", rows, max_nnz=4000, gap=20)
    assert abs(m - a.raw.X[rows]).sum() == 0


def test_fisher_tables():
    df = pd.DataFrame({"ct": ["a"] * 100 + ["b"] * 100, "g": (["x"] * 50 + ["y"] * 50) * 2,
                       "detected": [True] * 80 + [False] * 20 + [True] * 10 + [False] * 90})
    vr, bt = A.fisher_tables(df, "ct", "g")
    row = vr[(vr.scope == "all groups") & (vr.ct == "a")].iloc[0]
    assert row.pos_in_celltype == 80 and row.neg_in_other_cells == 90 and row.odds_ratio > 1
    assert len(bt) == 2 and bt.fdr.between(0, 1).all()


def test_full_analysis_and_summary(example, tmp_path):
    from cxg_explorer import notebook_agent as NA
    from cxg_explorer import summary as SM
    NA.OUT_DIR = tmp_path / "res"
    D = A.CxgDataset(example)
    vid, sym, _ = D.resolve_gene("example1")
    assert sym == "EXAMPLE1"
    res = NA.run_analysis(D, vid, sym, ["classical monocyte"], study="Test")
    out = res["folder"]
    ov = pd.read_csv(out / "1_cellular_level/overall_summary.csv").iloc[0]
    assert ov.n_cells_total == 2500 and 0 < ov.n_cells_expressing < 2500
    assert all(ok for _, ok, _ in A.qc_checks(out)), A.qc_checks(out)
    for f in ["report.html", "HOW_IT_IS_CALCULATED.html", "METHODS.md", "study_info.json", "samples_used.csv"]:
        assert (out / f).exists(), f
    runs = SM.discover_runs(tmp_path / "res")
    keys, labels = SM.study_keys(runs, "EXAMPLE1")
    S = SM.build_summary({k: v[0] for k, v in keys.items()}, labels, "EXAMPLE1",
                         mapping_path=tmp_path / "res/_cross_study/map.csv", online_citations=False)
    assert len(S["inventory"]) == 1 and "classical monocyte" in set(S["celltype_long"].cell_type)
    C = SM.build_celltype_summary(S, ("fine", "classical monocyte"))
    assert C["overview"].present.all()


def test_no_pictograms_in_source():
    pat = re.compile("[\U0001F300-\U0001FAFF\u2600-\u27BF\u2713\u2717\u2460-\u2464]")
    for f in PKG.glob("*.py"):
        assert not pat.search(f.read_text()), f"pictogram found in {f.name}"
