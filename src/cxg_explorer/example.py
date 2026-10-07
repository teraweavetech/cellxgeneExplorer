# cellxgeneExplorer
# Copyright (c) 2026 TeraWeave Technologies. Licensed under the PolyForm Noncommercial License 1.0.0 (see LICENSE.md).
# Required Notice: Copyright (c) 2026 TeraWeave Technologies
"""Small synthetic CELLxGENE-style dataset for trying the tools offline (no download needed)."""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def make_example(path="example_data/example_pbmc.h5ad", seed=7, n_cells=6000, n_genes=300):
    """Write a CELLxGENE-schema-like H5AD: 4 cell types, 8 donors, normal vs COVID-19, raw counts in raw.X.

    The gene EXAMPLE1 is enriched in classical monocytes (more so in COVID-19), and genes G1-G5 are
    co-expressed with it (G6-G7 anti-correlated), so every part of the analysis, including DE, has
    something to find.
    """
    import anndata as ad
    rng = np.random.default_rng(seed)
    cts = ["classical monocyte", "CD4-positive, alpha-beta T cell", "B cell", "natural killer cell"]
    ct = rng.choice(cts, n_cells, p=[.35, .35, .15, .15])
    donors = np.array([f"P{i}" for i in range(8)])
    donor = rng.choice(donors, n_cells)
    dis = np.where(np.isin(donor, donors[:4]), "normal", "COVID-19")
    base = rng.gamma(0.6, 0.5, n_genes)
    lib = rng.lognormal(0, 0.35, n_cells)
    deff = dict(zip(donors, rng.normal(0, 0.3, 8)))
    z = (ct == "classical monocyte") & (rng.random(n_cells) < np.where(dis == "COVID-19", 0.6, 0.35))
    lam = np.outer(lib * np.exp([deff[d] for d in donor]), base)
    lam[:, 0] = lib * np.where(z, 2.5, np.where(ct == "classical monocyte", 0.15, 0.05))
    lam[np.ix_(np.where(z)[0], [1, 2, 3, 4, 5])] *= 3.0
    lam[np.ix_(np.where(z)[0], [6, 7])] *= 0.3
    counts = sp.csr_matrix(rng.poisson(lam * rng.gamma(5, 1 / 5, lam.shape)).astype(np.float32))
    libs = np.asarray(counts.sum(1)).ravel()
    libs[libs == 0] = 1
    X = counts.multiply(1 / libs[:, None]).tocsr()
    X.data = np.log1p(X.data * 1e4)
    var = pd.DataFrame({"feature_name": ["EXAMPLE1"] + [f"G{i}" for i in range(1, n_genes)]},
                       index=[f"ENSG{i:011d}" for i in range(n_genes)])
    obs = pd.DataFrame({"cell_type": ct, "donor_id": donor, "disease": dis, "tissue": "blood",
                        "tissue_type": "tissue", "suspension_type": "cell", "sex": "female",
                        "assay": "10x 3' v3", "organism": "Homo sapiens"}, index=[f"c{i}" for i in range(n_cells)])
    obs["sample_id"] = obs["donor_id"] + "_" + rng.choice(["A", "B"], n_cells)
    obs = obs.astype("category")
    cen = {c: rng.normal(0, 10, 2) for c in cts}
    a = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var,
                   obsm={"X_umap": np.vstack([cen[c] for c in ct]) + rng.normal(0, 1.5, (n_cells, 2))})
    a.raw = ad.AnnData(X=counts, var=var)
    a.uns["title"] = "Example PBMC dataset (synthetic)"
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path)
    return path
