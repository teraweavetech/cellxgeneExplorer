# cellxgeneExplorer
# Copyright (c) 2026 TeraWeave Technologies. Licensed under the PolyForm Noncommercial License 1.0.0 (see LICENSE.md).
# Required Notice: Copyright (c) 2026 TeraWeave Technologies
"""Analysis engine: CELLxGENE discovery and download, three-level gene analysis, DE (NEBULA + Fisher),
methods text, calculation guide and reports."""
import argparse
import base64
import difflib
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import h5py
import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

try:
    from tqdm.auto import tqdm
except ImportError:  # tqdm is optional
    tqdm = None

CURATION_API = "https://api.cellxgene.cziscience.com/curation/v1"
CROSSREF_API = "https://api.crossref.org/works"
COLLECTION_URL = "https://cellxgene.cziscience.com/collections/"
UUID_RE = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", re.I)
DOI_RE = re.compile(r"10\.\d{4,9}/\S+", re.I)
HEADERS = {"User-Agent": "cellxgene-gene-agent/1.0 (research tool)"}

INTERACTIVE = False
FIGS = []     # (section, caption, path)
TABLES = []   # (section, caption, csv path, preview DataFrame)
NOTES = []    # (section, html text)


# ----------------------------------------------------------------------------- utils
def say(msg=""):
    print(msg, flush=True)


def header(msg):
    say("\n" + "=" * 78 + f"\n  {msg}\n" + "=" * 78)


def ask(prompt, default=None):
    if not INTERACTIVE:
        return default
    sfx = f" [{default}]" if default not in (None, "") else ""
    try:
        ans = input(f"\n» {prompt}{sfx}: ").strip()
    except EOFError:
        ans = ""
    return ans if ans else default


def safe_name(s, maxlen=60):
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_")
    return s[:maxlen] or "unnamed"


def labels_of(lst):
    if not lst:
        return ""
    if isinstance(lst, dict):
        lst = [lst]
    return ", ".join(sorted({str(x.get("label", "")) for x in lst if isinstance(x, dict)}))


def human_size(n):
    if not n:
        return "?"
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def get_json(url, params=None, timeout=120, retries=3):
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            if attempt == retries - 1:
                raise
            say(f"   (request failed: {e}; retrying)")
            time.sleep(2 * (attempt + 1))


# ----------------------------------------------------------------------------- discovery (CELLxGENE only)
def norm_words(s):
    return re.sub(r"[^a-z0-9 ]+", " ", (s or "").lower()).split()


def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    d = re.sub(r"^(https?://)?(dx\.)?doi\.org/", "", d)
    return re.sub(r"^doi:\s*", "", d).rstrip(".")


def title_similarity(query, target):
    q, t = norm_words(query), norm_words(target)
    if not q or not t:
        return 0.0
    seq = difflib.SequenceMatcher(None, " ".join(q), " ".join(t)).ratio()
    coverage = len(set(q) & set(t)) / len(set(q))
    return 0.5 * seq + 0.5 * coverage


def collection_dois(c):
    dois = {norm_doi(c.get("doi"))}
    for link in c.get("links") or []:
        if str(link.get("link_type", "")).upper() == "DOI" or "doi.org" in str(link.get("link_url", "")):
            dois.add(norm_doi(link.get("link_url")))
    return {d for d in dois if d}


def fetch_all_collections(cache_dir):
    cache = cache_dir / "cellxgene_collections.json"
    if cache.exists() and time.time() - cache.stat().st_mtime < 86400:
        return json.loads(cache.read_text())
    say("   Fetching the CELLxGENE Discover collection index ...")
    data = get_json(f"{CURATION_API}/collections")
    cache.write_text(json.dumps(data))
    return data


def crossref_dois(title):
    """Resolve a paper title to DOI(s) via Crossref, used only to match CELLxGENE collections."""
    try:
        js = get_json(CROSSREF_API, params={"query.bibliographic": title, "rows": 5}, timeout=30, retries=2)
    except Exception:
        return []
    out = []
    for it in js.get("message", {}).get("items", []):
        t = (it.get("title") or [""])[0]
        if title_similarity(title, t) >= 0.8:
            out.append(norm_doi(it.get("DOI")))
    return out


def find_collections(query, cache_dir, top=8):
    m = UUID_RE.search(query)
    if m:
        return [(1.0, get_json(f"{CURATION_API}/collections/{m.group(0)}"), "collection ID")]

    collections = fetch_all_collections(cache_dir)
    doi_m = DOI_RE.search(query)
    target_dois = {norm_doi(doi_m.group(0))} if doi_m else set(crossref_dois(query))

    scored = []
    for c in collections:
        why = ""
        if target_dois & collection_dois(c):
            score, why = 1.5, "DOI match"
        else:
            name_s = title_similarity(query, c.get("name", ""))
            ds_s = max([title_similarity(query, d.get("title", "")) for d in (c.get("datasets") or []) if d] or [0])
            desc_s = 0.7 * title_similarity(query, c.get("description", ""))
            score = max(name_s, ds_s, desc_s)
            why = "collection title" if score == name_s else ("dataset title" if score == ds_s else "description")
        scored.append((score, c, why))
    scored.sort(key=lambda x: -x[0])
    return scored[:top]


def choose_collection(query, cache_dir):
    hits = find_collections(query, cache_dir)
    if not hits:
        sys.exit("No CELLxGENE collections found.")
    say("\nCandidate CELLxGENE collections:")
    for i, (score, c, why) in enumerate(hits, 1):
        pm = c.get("publisher_metadata") or {}
        journal = f"{pm.get('journal', '')} {pm.get('published_year', '')}".strip()
        n_ds = len(c.get("datasets") or [])
        say(f"  [{i}] ({why}, score {min(score, 1.0):.2f}) {c.get('name')}")
        say(f"       {journal} | DOI: {c.get('doi') or '-'} | {n_ds} dataset(s)")
    if hits[0][0] < 0.45:
        say("\n   No strong match. The paper may not be on CELLxGENE, or try its DOI instead of the title.")
    pick = int(ask("Pick a collection number", "1")) - 1
    c = hits[pick][1]
    if not c.get("datasets") or "assets" not in (c["datasets"][0] or {}):
        c = get_json(f"{CURATION_API}/collections/{c['collection_id']}")
    return c


def choose_dataset(collection, preset=None):
    ds = [d for d in collection.get("datasets", []) if d]
    if not ds:
        sys.exit("This collection has no datasets.")
    ds.sort(key=lambda d: -(d.get("cell_count") or 0))
    say(f"\nDatasets in: {collection.get('name')}")
    for i, d in enumerate(ds, 1):
        h5 = next((a for a in d.get("assets", []) if a.get("filetype") == "H5AD"), {})
        say(f"  [{i}] {d.get('title')}")
        say(f"       {d.get('cell_count', '?'):,} cells | {labels_of(d.get('organism'))} | "
            f"{labels_of(d.get('assay'))} | H5AD {human_size(h5.get('filesize'))}")
        tissues = labels_of(d.get("tissue"))
        if tissues:
            say(f"       tissue: {tissues[:110]}")
    pick = int(preset or ask("Pick a dataset number", "1")) - 1
    return ds[pick]


def download_h5ad(dataset, data_dir):
    asset = next((a for a in dataset.get("assets", []) if a.get("filetype") == "H5AD"), None)
    if not asset:
        sys.exit("No H5AD asset for this dataset.")
    url, size = asset["url"], asset.get("filesize") or 0
    dest = data_dir / f"{safe_name(dataset.get('title'), 50)}__{dataset['dataset_id'][:8]}.h5ad"
    if dest.exists() and (not size or dest.stat().st_size == size):
        say(f"   Using cached file {dest}")
        return dest
    if size > 5 * 1024 ** 3:
        if str(ask(f"File is {human_size(size)}. Download anyway? (y/n)", "y")).lower() != "y":
            sys.exit("Aborted.")
    say(f"   Downloading {human_size(size)} from CELLxGENE ...")
    tmp = dest.with_suffix(".part")
    with requests.get(url, stream=True, headers=HEADERS, timeout=600) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", size) or 0)
        bar = tqdm(total=total, unit="B", unit_scale=True) if tqdm else None
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                fh.write(chunk)
                if bar:
                    bar.update(len(chunk))
        if bar:
            bar.close()
    tmp.rename(dest)
    return dest


# ----------------------------------------------------------------------------- reading the H5AD
def _decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def read_h5_index(group):
    """Index of an encoded dataframe group (works across anndata on-disk versions)."""
    try:
        from anndata.io import read_elem
    except ImportError:
        from anndata.experimental import read_elem
    return np.asarray(read_elem(group).index.astype(str))


def read_gene_column(f, path, gidx, n_obs, want_rowsums=False, max_nnz=20_000_000):
    """Return (column vector, row sums or None) for one gene, reading the matrix in chunks."""
    node = f[path]
    col = np.zeros(n_obs, dtype=np.float32)
    rs = np.zeros(n_obs, dtype=np.float64) if want_rowsums else None

    if isinstance(node, h5py.Dataset):  # dense
        step = max(1, max_nnz // max(1, node.shape[1]))
        for s in range(0, n_obs, step):
            block = node[s:s + step]
            col[s:s + step] = block[:, gidx]
            if want_rowsums:
                rs[s:s + step] = block.sum(1)
        return col, rs

    enc = node.attrs.get("encoding-type", "csr_matrix")
    enc = enc.decode() if isinstance(enc, bytes) else enc
    data, indices, indptr = node["data"], node["indices"], node["indptr"][:]

    if enc == "csc_matrix":
        a, b = indptr[gidx], indptr[gidx + 1]
        col[indices[a:b]] = data[a:b]
        if want_rowsums:
            for s in range(0, len(data), max_nnz):
                rs += np.bincount(indices[s:s + max_nnz], weights=data[s:s + max_nnz], minlength=n_obs)
        return col, rs

    # CSR: walk rows in chunks bounded by number of non-zeros
    bar = tqdm(total=int(indptr[-1]), unit="nnz", unit_scale=True, desc="   reading") if tqdm else None
    s = 0
    while s < n_obs:
        e = int(np.searchsorted(indptr, indptr[s] + max_nnz, side="right")) - 1
        e = min(max(e, s + 1), n_obs)
        a, b = indptr[s], indptr[e]
        idx, dat = indices[a:b], data[a:b]
        rows = np.repeat(np.arange(e - s), np.diff(indptr[s:e + 1]))
        hit = idx == gidx
        col[s + rows[hit]] = dat[hit]
        if want_rowsums:
            rs[s:e] = np.bincount(rows, weights=dat, minlength=e - s)
        if bar:
            bar.update(int(b - a))
        s = e
    if bar:
        bar.close()
    return col, rs


def looks_like_counts(f, path, n=200_000):
    node = f[path]
    sample = node[: min(2000, node.shape[0])] if isinstance(node, h5py.Dataset) else node["data"][:n]
    sample = np.asarray(sample).ravel()
    sample = sample[sample != 0]
    return sample.size > 0 and np.allclose(sample, np.round(sample)) and sample.max() > 10


class CxgDataset:
    def __init__(self, path):
        self.path = Path(path)
        say(f"   Opening {self.path.name} (backed, memory-light) ...")
        self.adata = ad.read_h5ad(self.path, backed="r")
        self.obs = self.adata.obs.copy()
        self.var = self.adata.var.copy()
        self.n_obs = self.adata.n_obs
        self.h5 = h5py.File(self.path, "r")
        self.has_raw = "raw" in self.h5 and "X" in self.h5["raw"]
        self.raw_var_names = read_h5_index(self.h5["raw"]["var"]) if self.has_raw else None
        self.embedding_key = self._pick_embedding()
        self.title = str(self.adata.uns.get("title", self.path.stem)) if "title" in self.adata.uns else self.path.stem
        say(f"   {self.n_obs:,} cells x {self.adata.n_vars:,} genes | embedding: {self.embedding_key} | "
            f"raw counts layer: {'yes' if self.has_raw else 'no'}")

    def _pick_embedding(self):
        keys = list(self.adata.obsm.keys())
        for pref in ["X_umap", "umap"]:
            for k in keys:
                if k.lower() == pref.lower():
                    return k
        for k in keys:
            if "umap" in k.lower():
                return k
        for k in keys:
            if "tsne" in k.lower():
                return k
        return keys[0] if keys else None

    def embedding(self):
        if not self.embedding_key:
            return None
        return np.asarray(self.adata.obsm[self.embedding_key])[:, :2].astype(float)

    def symbols(self):
        if "feature_name" in self.var.columns:
            return self.var["feature_name"].astype(str)
        return pd.Series(self.var.index.astype(str), index=self.var.index)

    def resolve_gene(self, query):
        q = query.strip()
        sym = self.symbols()
        ids_nover = pd.Series(self.var.index.astype(str).str.replace(r"\.\d+$", "", regex=True), index=self.var.index)
        for series in (sym, ids_nover):
            hits = series.index[series.str.upper() == q.upper()]
            if len(hits):
                if len(hits) > 1:
                    say(f"   Note: {len(hits)} features match '{q}'; using {hits[0]}.")
                return hits[0], sym[hits[0]], []
        # Some datasets append the Ensembl ID to the symbol, e.g. 'GENE_ENSG...'
        hits = sym.index[sym.str.upper().str.split("_").str[0] == q.upper()]
        if len(hits):
            return hits[0], sym[hits[0]], []
        sugg = difflib.get_close_matches(q.upper(), sym.str.upper().unique().tolist(), n=8, cutoff=0.6)
        return None, None, sugg

    def gene_expression(self, var_id):
        """Return (expression used for plots, detected boolean, description of the values)."""
        gidx = int(np.where(self.var.index == var_id)[0][0])
        x_counts = looks_like_counts(self.h5, "X")
        say("   Extracting expression from X ...")
        expr, rs = read_gene_column(self.h5, "X", gidx, self.n_obs, want_rowsums=x_counts)
        if x_counts:
            lib = np.where(rs > 0, rs, 1.0)
            detected = expr > 0
            expr = np.log1p(expr / lib * 1e4).astype(np.float32)
            desc = "log1p(CP10k), computed from the raw counts in X"
            PARAMS["x_def"] = "x_i = ln(1 + 10^4 × c_i / L_i), L_i = total counts of cell i (X held raw counts)"
        else:
            desc = "normalized values in X as deposited by the authors (CELLxGENE schema; usually log-normalized)"
            PARAMS["x_def"] = "x_i = normalized expression of the gene in cell i (CELLxGENE X matrix as deposited)"
            detected = expr > 0
            if self.has_raw and var_id in set(self.raw_var_names):
                ridx = int(np.where(self.raw_var_names == var_id)[0][0])
                say("   Extracting raw counts from raw.X for detection calls ...")
                raw_col, _ = read_gene_column(self.h5, "raw/X", ridx, self.n_obs)
                detected = raw_col > 0
        return expr, detected, desc



# ----------------------------------------------------------------------------- methods & AI interpretation
GROUP_COL = "disease"   # grouping/condition column; changeable in the notebook panel
PARAMS = {"min_cells": 3, "min_pct": 1.0, "min_group": 10, "x_def": "normalized values in X as deposited"}
STATS = []      # subject-level statistical tests run during the analysis
RUN_INFO = {}   # dataset / run context used by the Methods text and the interpretation

LLM_MODEL = "gpt-oss:20b"
LLM_URL = "http://localhost:11434"   # Ollama default; an OpenAI-compatible URL ending in /v1 also works


def de_methods():
    d = RUN_INFO.get("de")
    if not d:
        return ""
    g = RUN_INFO.get("gene", "the gene")
    scopes = []
    if "pooled" in d["scopes"]:
        scopes.append("with all selected groups pooled (group as a fixed-effect covariate)")
    if "within" in d["scopes"]:
        scopes.append("separately within each group")
    return f"""**Differential expression.** Within each cell type, cells were classified as {g}-positive (at least one \
count) or {g}-negative. Genes differentially expressed between positive and negative cells were identified with \
NEBULA (v{d.get('nebula_version', '?')}; negative binomial mixed model, NEBULA-LN), using raw counts, a random \
intercept per donor, each cell's total counts as the offset, and the standardized log number of detected genes \
per cell as a covariate. Analyses were run {' and '.join(scopes)}. A cell type was analysed when it had at least \
{d['min_cells']} positive and {d['min_cells']} negative cells and at least {d['min_donors']} donors contributing at least \
{d['min_cells_per_donor']} cells of each class. At most {d['max_cells']:,} cells per class were used, randomly \
downsampled with a fixed seed. Genes detected in at least {d['min_det'] * 100:g}% of positive or negative cells were \
tested, excluding {g} itself. P-values were adjusted with the Benjamini-Hochberg procedure within each analysis, \
using genes whose fit converged according to NEBULA's convergence code (code > -20); significance was set at FDR \
< {d['fdr']:g}.

**Cell proportions.** Within each cell type, the numbers of {g}-positive and {g}-negative cells were compared with \
two-sided Fisher's exact tests: (i) the cell type against all other cells (overall and within each group), and \
(ii) each pair of groups. Odds ratios with Woolf 95% confidence intervals (0.5 continuity correction) are reported, \
and p-values were Benjamini-Hochberg adjusted within each table. These tests treat cells as independent \
observations and therefore do not account for donor-to-donor variation.

"""


def methods_markdown():
    """Methods section written in prose, filled in with this run's details."""
    r, p = RUN_INFO, PARAMS
    norm = ("Expression values were the normalized matrix (X) deposited on CELLxGENE, used without "
            "further transformation." if "deposited" in p["x_def"] else
            "Because X contained raw counts, values were normalized to 10,000 counts per cell and log-transformed "
            "(log1p).")
    det = ("A cell was considered to express the gene when at least one UMI/read for the gene was detected in the raw "
           "counts (raw.X)." if r.get("has_raw") else
           "A cell was considered to express the gene when its value in X was greater than zero.")
    pops = r.get("populations") or []
    pop_txt = (f"For the cell-specific analysis, the following population(s) were selected: {', '.join(pops)}"
               + ("; they were analysed individually and as one combined group." if len(pops) > 1 else ".")
               if pops else "No population was selected for the cell-specific analysis.")
    grp = r.get("groups")
    grp_txt = (f" The analysis was restricted to cells whose `{GROUP_COL}` was one of: {', '.join(grp)}."
               if r.get("groups_subset") else
               (f" All `{GROUP_COL}` groups were included ({', '.join(grp)})." if grp else ""))
    return f"""**Data source.** The dataset "{r.get('dataset', '?')}" was obtained from CELLxGENE Discover \
(collection: {r.get('collection', 'n/a')}; DOI: {r.get('doi', 'n/a')}) as an H5AD file. It contains \
{r.get('n_cells', 0):,} cells from {r.get('n_subjects', '?')} subjects across {r.get('n_celltypes', '?')} cell types. \
Cell-type labels (`{r.get('ct_col', 'cell_type')}`), subject identifiers (`{r.get('donor_col', 'donor_id')}`), disease, \
tissue and sex annotations, and the low-dimensional embedding (`{r.get('embedding', 'X_umap')}`) were used exactly as \
curated on CELLxGENE; no re-clustering, re-annotation or re-embedding was performed.{grp_txt}

**Gene expression.** {r.get('gene', '')} was quantified per cell. {norm} {det} The gene's expression was read directly \
from the H5AD file on disk.

**Cellular level.** For every cell type, the analysis recorded the number of cells and the number and proportion of \
cells expressing the gene. It also recorded the mean expression across all cells and within expressing cells only, \
and the median within expressing cells. The same summaries were produced for disease, tissue and sex where these \
varied. Expression was displayed on the authors' embedding: non-expressing cells are shown in grey, and expressing \
cells are coloured by expression level, with the colour scale capped at the 99th percentile. Distributions among \
expressing cells were shown as violin plots, and the proportion and mean level per cell type and condition as dot plots.

**Cell-specific level.** {pop_txt} Within each population, the same cell counts, proportions and expression \
summaries were computed, broken down by disease, tissue and sex, and repeated per subject.

**Sample (subject) level.** Each subject was summarised by its total cells, expressing cells and proportion \
expressing. Where a donor contributed samples to more than one `{GROUP_COL}` group, each donor-group pair was \
treated as a separate unit. A subject was classified as expressing the gene when at least {p['min_cells']} of its \
cells, and at least {p['min_pct']}% of its cells, expressed it. The number and proportion of expressing subjects \
were reported overall and per `{GROUP_COL}` group. To show how sensitive this classification is to the choice of \
threshold, it was repeated under alternative cut-offs: at least 1, 3 or 10 expressing cells, and at least 1%, 5% or \
10% of cells. Per-cell-type subject counts only included subjects with at least {p['min_group']} cells of that type. \
`{GROUP_COL}` groups were compared on per-subject proportions, using a two-sided Mann-Whitney U test for two groups or a \
Kruskal-Wallis test for more. This uses each subject as one replicate rather than treating individual cells as \
independent. These tests are exploratory and not corrected for multiple testing.

{de_methods()}**Descriptive report.** On request, a narrative report was drafted by the open-weight language model {LLM_MODEL}, \
run locally, using only the summary statistics above. It should be checked against the tables and figures."""


def _records(path, cols=None, n=15):
    path = Path(path)
    if not path.exists():
        return None
    d = pd.read_csv(path)
    if cols:
        d = d[[c for c in cols if c in d.columns]]
    return d.head(n).round(4).to_dict("records")


def build_digest(out):
    """Compact numeric summary of each level, used as the only input to the interpretation model."""
    out = Path(out)
    keep = ["n_cells", "n_expressing", "pct_expressing", "mean_expr_expressing_cells"]
    c1, c2, c3 = out / "1_cellular_level", out / "2_cell_specific_level", out / "3_sample_level"
    ctx = {k: RUN_INFO.get(k) for k in ["gene", "dataset", "collection", "n_cells", "n_subjects", "n_celltypes",
                                       "tissues", "diseases", "assays"]}
    ctx["grouping_column"] = GROUP_COL
    ctx["groups_analysed"] = RUN_INFO.get("groups")
    ctx["subject_rule"] = f">= {PARAMS['min_cells']} expressing cells and >= {PARAMS['min_pct']}% of cells"
    lvl1 = {"overall": _records(c1 / "overall_summary.csv", n=1),
            "cell_types_ranked_by_pct_expressing": _records(c1 / "celltype_summary.csv",
                                                            [RUN_INFO.get("ct_col", "cell_type")] + keep, 20)}
    for cov in dict.fromkeys([GROUP_COL, "disease", "tissue", "sex"]):
        r = _records(c1 / f"summary_by_{cov}.csv", [cov] + keep, 12)
        if r:
            lvl1[f"by_{cov}"] = r
    lvl2 = {}
    if c2.exists():
        lvl2["populations"] = _records(c2 / "selected_populations_summary.csv",
                                       ["population", "n_cells_total", "n_cells_expressing", "pct_expressing",
                                        "mean_expr_expressing_cells"], 10)
        for d in sorted(p for p in c2.iterdir() if p.is_dir()):
            e = {}
            for cov in dict.fromkeys([GROUP_COL, "disease", "tissue", "sex"]):
                r = _records(d / f"summary_by_{cov}.csv", [cov] + keep, 10)
                if r:
                    e[f"by_{cov}"] = r
            h = _records(d / "sample_level/subjects_expressing_headline.csv",
                         ["group", "n_subjects_expressing", "n_subjects_total", "pct_expressing"], 10)
            if h:
                e["subjects_expressing"] = h
            lvl2[d.name] = e
    lvl3 = {}
    if c3.exists():
        lvl3["subjects_expressing_chosen_rule"] = _records(
            c3 / "subjects_expressing_headline.csv",
            ["group", "n_subjects_expressing", "n_subjects_total", "pct_expressing"], 10)
        lvl3["sensitivity_to_threshold"] = _records(c3 / "subjects_expressing_overview.csv", n=10)
        lvl3["per_subject"] = _records(c3 / "subject_summary.csv",
                                       [RUN_INFO.get("donor_col", "donor_id"), GROUP_COL, "n_cells", "n_expressing",
                                        "pct_expressing", "expressing_subject"], 40)
        lvl3["subjects_expressing_per_cell_type"] = _records(
            c3 / "celltype_n_subjects_expressing.csv",
            [RUN_INFO.get("ct_col", "cell_type"), "n_subjects_with_celltype", "n_subjects_expressing",
             "pct_subjects_expressing"], 20)
    if STATS:
        lvl3["group_tests_on_subjects"] = [s for s in STATS if s["level"].startswith("3.")]
        if lvl2:
            lvl2["group_tests_on_subjects"] = [s for s in STATS if s["level"].startswith("2.")]
    lvl4 = {}
    c4 = out / "4_DE"
    if c4.exists():
        lvl4["nebula_runs"] = _records(c4 / "nebula_summary.csv", n=40)
        lvl4["fisher_celltype_vs_rest"] = _records(c4 / "fisher_celltype_vs_rest.csv",
                                                   ["scope", RUN_INFO.get("ct_col", "cell_type"), "pct_pos_celltype",
                                                    "pct_pos_other_cells", "odds_ratio", "fdr"], 40)
        lvl4["fisher_between_groups"] = _records(c4 / "fisher_between_groups.csv",
                                                 [RUN_INFO.get("ct_col", "cell_type"), "group_1", "group_2",
                                                  "pct_pos_group_1", "pct_pos_group_2", "odds_ratio_g1_vs_g2", "fdr"], 40)
    return {"context": ctx, "1_cellular": lvl1, "2_cell_specific": lvl2, "3_sample": lvl3, "4_de": lvl4}


def ask_llm(messages, model=None, url=None, timeout=900):
    model, url = model or LLM_MODEL, (url or LLM_URL).rstrip("/")
    if url.endswith("/v1"):   # OpenAI-compatible server (vLLM, LM Studio, llama.cpp server, ...)
        r = requests.post(f"{url}/chat/completions", timeout=timeout,
                          json={"model": model, "messages": messages, "temperature": 0.2})
        r.raise_for_status()
        txt = r.json()["choices"][0]["message"]["content"]
    else:                     # Ollama
        r = requests.post(f"{url}/api/chat", timeout=timeout,
                          json={"model": model, "messages": messages, "stream": False,
                                "options": {"temperature": 0.2, "num_ctx": 16384}})
        r.raise_for_status()
        txt = r.json()["message"]["content"]
    txt = re.sub(r"<think>.*?</think>", "", txt or "", flags=re.S)
    return txt.split("assistantfinal")[-1].strip()


def check_llm(model=None, url=None):
    """Return (ok, message) for the configured model server."""
    model, url = model or LLM_MODEL, (url or LLM_URL).rstrip("/")
    try:
        if url.endswith("/v1"):
            names = [m["id"] for m in requests.get(f"{url}/models", timeout=5).json().get("data", [])]
        else:
            names = [m["name"] for m in requests.get(f"{url}/api/tags", timeout=5).json().get("models", [])]
    except Exception:
        return False, (f"No model server at {url}. Install Ollama (https://ollama.com), then run "
                       f"<code>ollama pull {model}</code> (and <code>ollama serve</code> if it is not running).")
    if not any(n == model or n.startswith(model + ":") or n.split(":")[0] == model.split(":")[0] and
               model.split(":")[-1] in n for n in names):
        return False, f"Server found but <code>{model}</code> is not installed. Run <code>ollama pull {model}</code>."
    return True, "ok"


def md_to_html(md):
    md = re.sub(r"(?m)^(?![-*] )(.+)\n(?=[-*] )", r"\1\n\n", md or "")   # blank line before lists
    try:
        import markdown as _md
        return _md.markdown(md)
    except ImportError:
        h = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", md)
        h = re.sub(r"`(.+?)`", r"<code>\1</code>", h)
        h = re.sub(r"(?m)^[-*] (.*)$", r"<li>\1</li>", h)
        return "<p>" + h.replace("\n\n", "</p><p>").replace("\n", " ") + "</p>"


# ----------------------------------------------------------------------------- "How it's calculated" (QC page)
def _f(v, d=4):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    return f"{int(v):,}" if v == int(v) and abs(v) >= 1 else f"{v:.{d}g}"


def _rows(items):
    """items: list of (column, what it is, how it is calculated, worked example)"""
    h = ('<table class="calc"><tr><th>Column</th><th>What it means</th><th>How it is calculated</th>'
         '<th>Worked example from your data</th></tr>')
    for col, what, how, ex in items:
        h += f"<tr><td><code>{col}</code></td><td>{what}</td><td>{how}</td><td>{ex}</td></tr>"
    return h + "</table>"


def _group_items(r, label, total_pos, gene, unit="group"):
    n, k = r["n_cells"], r["n_expressing"]
    mp, ma = r["mean_expr_expressing_cells"], r["mean_expr_all_cells"]
    return [
        ("n_cells", f"Number of cells in the {unit}.", f"Count of cells whose label is the {unit} name.",
         f"cells labelled <b>{label}</b> = <b>{_f(n)}</b>"),
        ("n_expressing", f"Cells in which {gene} was detected.",
         f"Count of cells in the {unit} with at least 1 count of {gene}.", f"<b>{_f(k)}</b> of the {_f(n)} cells"),
        ("n_not_expressing", f"Cells with no {gene} detected.", "n_cells − n_expressing",
         f"{_f(n)} − {_f(k)} = <b>{_f(r['n_not_expressing'])}</b>"),
        ("proportion_expressing", "Fraction of the cells that express the gene (0–1).", "n_expressing ÷ n_cells",
         f"{_f(k)} ÷ {_f(n)} = <b>{_f(r['proportion_expressing'])}</b>"),
        ("pct_expressing", "Same as above, as a percentage.", "proportion_expressing × 100",
         f"{_f(r['proportion_expressing'])} × 100 = <b>{_f(r['pct_expressing'])}%</b>"),
        ("mean_expr_all_cells", f"Average {gene} level over <i>all</i> cells of the {unit}, counting non-expressing "
         "cells as 0.", "sum of expression values of all cells ÷ n_cells",
         f"sum over {_f(n)} cells ÷ {_f(n)} = <b>{_f(ma)}</b><br><small>check: ≈ proportion × mean in expressing = "
         f"{_f(r['proportion_expressing'])} × {_f(mp)} = {_f(r['proportion_expressing'] * mp)}</small>"),
        ("mean_expr_expressing_cells", f"Average {gene} level among the expressing cells only (how high it is "
         "where it is present).", "sum of expression values of expressing cells ÷ n_expressing",
         f"average over the {_f(k)} expressing cell(s) = <b>{_f(mp)}</b>"),
        ("median_expr_expressing_cells", "Middle value of expression among the expressing cells (less affected by "
         "a few very high cells).", "sort the expressing cells' values and take the middle one",
         f"median of the {_f(k)} expressing cell(s) = <b>{_f(r['median_expr_expressing_cells'])}</b>"),
        ("pct_of_all_expressing_cells", f"Share of <i>all</i> {gene}-expressing cells in the dataset that belong to "
         f"this {unit}.", "n_expressing ÷ (expressing cells in whole dataset) × 100",
         f"{_f(k)} ÷ {_f(total_pos)} × 100 = <b>{_f(r['pct_of_all_expressing_cells'])}%</b>"),
    ]


def qc_checks(out):
    """Automatic consistency checks on the saved tables."""
    out, res = Path(out), []
    c1 = out / "1_cellular_level"
    ov = pd.read_csv(c1 / "overall_summary.csv").iloc[0]
    ct = pd.read_csv(c1 / "celltype_summary.csv")
    add = lambda name, ok, detail: res.append((name, bool(ok), detail))
    add("Expressing + not expressing = total (every cell type)",
        (ct.n_expressing + ct.n_not_expressing == ct.n_cells).all(), f"{len(ct)} cell types checked")
    add("Cells summed over cell types = total cells", ct.n_cells.sum() == ov.n_cells_total,
        f"{ct.n_cells.sum():,} vs {int(ov.n_cells_total):,}")
    add("Expressing cells summed over cell types = total expressing", ct.n_expressing.sum() == ov.n_cells_expressing,
        f"{ct.n_expressing.sum():,} vs {int(ov.n_cells_expressing):,}")
    add("% expressing = 100 × expressing ÷ cells (every cell type)",
        np.allclose(ct.pct_expressing, 100 * ct.n_expressing / ct.n_cells, atol=1e-2), "tolerance 0.01 percentage points")
    add("pct_of_all_expressing_cells sums to 100%", abs(ct.pct_of_all_expressing_cells.sum() - 100) < 0.1
        or ov.n_cells_expressing == 0, f"sum = {ct.pct_of_all_expressing_cells.sum():.2f}%")
    approx = ct.proportion_expressing * ct.mean_expr_expressing_cells
    add("mean (all cells) ≈ proportion × mean (expressing cells)",
        np.allclose(ct.mean_expr_all_cells, approx, atol=2e-3, rtol=0.02),
        "holds when non-expressing cells have value 0 (small rounding differences are expected)")
    ss = out / "3_sample_level/subject_summary.csv"
    if ss.exists():
        s = pd.read_csv(ss)
        add("Cells summed over subjects = total cells", s.n_cells.sum() == ov.n_cells_total,
            f"{s.n_cells.sum():,} vs {int(ov.n_cells_total):,}")
        rule = (s.n_expressing >= PARAMS["min_cells"]) & (s.pct_expressing >= PARAMS["min_pct"])
        add("expressing_subject matches the rule for every subject", (rule == s.expressing_subject).all(),
            f">= {PARAMS['min_cells']} cells and >= {PARAMS['min_pct']}%")
    return res


def de_calc_html(d4, gene, d):
    h = [f"<h3>DE: Fisher's exact tests <small>(fisher_celltype_vs_rest.csv, fisher_between_groups.csv)</small></h3>"]
    vr = d4 / "fisher_celltype_vs_rest.csv"
    if vr.exists():
        v = pd.read_csv(vr).iloc[0]
        ctc = [c for c in v.index if c not in ("scope",)][0]
        h.append(_rows([
            ("2×2 table (cell type vs rest)", "Rows: this cell type / all other cells. Columns: gene+ / gene- cells.",
             "a = gene+ in cell type, b = gene- in cell type, c = gene+ in other cells, d = gene- in other cells",
             f"{v[ctc]} ({v.scope}): a={_f(v.pos_in_celltype)}, b={_f(v.neg_in_celltype)}, "
             f"c={_f(v.pos_in_other_cells)}, d={_f(v.neg_in_other_cells)}"),
            ("odds_ratio", "How much more likely a cell is gene+ in this cell type than in the other cells.",
             "(a × d) ÷ (b × c) (conditional estimate from Fisher's test)", f"<b>{_f(v.odds_ratio)}</b>"),
            ("or_ci95_low / high", "95% confidence interval of the odds ratio.",
             "Woolf: exp(ln OR ± 1.96 × √(1/a + 1/b + 1/c + 1/d)), adding 0.5 to each cell",
             f"{_f(v.or_ci95_low)} – {_f(v.or_ci95_high)}"),
            ("p_value / fdr", "Two-sided Fisher's exact test; FDR = Benjamini-Hochberg within each table.",
             "fisher_exact([[a, b], [c, d]])", f"p = {_f(v.p_value)}, FDR = {_f(v.fdr)}"),
        ]))
        h.append("<p>Between groups: same test, with rows = group 1 / group 2 within one cell type "
                 "(odds ratio > 1 = more gene+ cells in group 1). These tests count cells, so they ignore "
                 "donor-to-donor variation and can be very significant with many cells.</p>")
    h.append(f"<h3>DE: NEBULA <small>(4_DE/&lt;cell type&gt;/&lt;scope&gt;/nebula_de_results.csv)</small></h3>")
    h.append(_rows([
        ("cells used", f"gene+ = ≥1 count of {gene}; gene- = 0 counts. At most {d['max_cells']:,} per class "
                       "(random, fixed seed).", "", ""),
        ("model", "Negative binomial mixed model per gene (NEBULA-LN).",
         "counts ~ status + z_log_ngenes [+ group], random intercept per donor, offset = log(total counts of the "
         "cell)", "status = gene+ vs gene- (reference: gene-)"),
        ("log2FC", "Fold change in expression, gene+ vs gene- cells, adjusted for the covariates.",
         "NEBULA coefficient for status (natural log) ÷ ln 2", ""),
        ("p_value / fdr", "Wald test of the status coefficient; FDR = Benjamini-Hochberg over reliable fits.", "", ""),
        ("reliable_fit / convergence", "NEBULA's convergence code; codes of -20 or lower indicate a possible "
                                       "convergence failure, and those genes are excluded from the FDR.", "", ""),
        ("pct_detected_pos / neg", "Share of gene+ / gene- cells with ≥1 count of this gene.", "", ""),
        ("mean_lognorm_pos / neg", "Average log1p(counts ÷ total counts × 10,000) in gene+ / gene- cells "
                                   "(descriptive only; the test uses NEBULA).", "", ""),
        ("genes tested", f"Detected in ≥ {d['min_det'] * 100:g}% of gene+ or gene- cells, excluding {gene}.", "", ""),
    ]))
    return "\n".join(h)


def calc_guide_html(out, gene):
    out = Path(out)
    c1, c3 = out / "1_cellular_level", out / "3_sample_level"
    ov = pd.read_csv(c1 / "overall_summary.csv").iloc[0]
    total_pos = int(ov.n_cells_expressing)
    ct = pd.read_csv(c1 / "celltype_summary.csv")
    ctc = ct.columns[0]
    r = ct.iloc[0]
    p = PARAMS
    h = [f'<div class="calcpage"><h2>How each number is calculated</h2>'
         f"<p>Each table is explained column by column, with a worked example that uses the first row of your own "
         f"results, so you can recompute any number by hand. The analysis gene is <b>{gene}</b>.</p>",
         "<h3>Basic definitions</h3><ul>"
         f"<li><b>Expressing cell</b>: a cell with at least 1 count of {gene}"
         + (" in the raw counts (raw.X)" if RUN_INFO.get("has_raw") else "") + ".</li>"
         f"<li><b>Expression value</b>: {p['x_def']}. All means and medians use this value.</li>"
         f"<li><b>Groups</b>: CELLxGENE's own annotations (<code>{RUN_INFO.get('ct_col', 'cell_type')}</code>, "
         f"<code>{RUN_INFO.get('donor_col', 'donor_id')}</code>, disease, tissue, sex), used as deposited.</li>"
         f"<li><b>Expressing subject</b>: at least <b>{p['min_cells']}</b> expressing cells <b>and</b> at least "
         f"<b>{p['min_pct']}%</b> of that subject's cells expressing. Both conditions must hold.</li></ul>",
         "<h3>Automatic consistency checks</h3>"
         '<table class="calc"><tr><th>Check</th><th>Result</th><th>Detail</th></tr>']
    for name, ok, det in qc_checks(out):
        h.append(f"<tr><td>{name}</td><td>{'pass' if ok else 'FAIL'}</td><td>{det}</td></tr>")
    h.append("</table>")

    h.append(f"<h3>1. Overall summary <small>(overall_summary.csv, population_summary.csv)</small></h3>"
             f"<p>Same idea as the cell-type table below, with all cells (or all cells of the selected population) "
             f"treated as one group. Here: {int(ov.n_cells_expressing):,} expressing ÷ {int(ov.n_cells_total):,} "
             f"cells × 100 = <b>{ov.pct_expressing}%</b>.</p>")
    h.append(f"<h3>2. Per cell type <small>(celltype_summary.csv; the same columns appear in summary_by_disease / "
             f"tissue / sex and in the cell-specific folders)</small></h3>"
             f"<p>Example row: <b>{r[ctc]}</b>. The dataset has {total_pos:,} expressing cells in total.</p>")
    h.append(_rows(_group_items(r, r[ctc], total_pos, gene, "cell type")))

    ss = c3 / "subject_summary.csv"
    if ss.exists():
        s = pd.read_csv(ss)
        sc = RUN_INFO.get("donor_col", s.columns[0])
        sr = s.iloc[0]
        items = _group_items(sr, sr[sc], total_pos, gene, "subject")
        c_ok, p_ok = sr.n_expressing >= p["min_cells"], sr.pct_expressing >= p["min_pct"]
        items.append(("expressing_subject", "Whether this subject counts as expressing the gene.",
                      f"True if n_expressing ≥ {p['min_cells']} <b>and</b> pct_expressing ≥ {p['min_pct']}%",
                      f"{_f(sr.n_expressing)} ≥ {p['min_cells']} {'(met)' if c_ok else '(not met)'} and "
                      f"{_f(sr.pct_expressing)}% ≥ {p['min_pct']}% {'(met)' if p_ok else '(not met)'}: "
                      f"<b>{bool(sr.expressing_subject)}</b>"))
        h.append(f"<h3>3. Per subject <small>(subject_summary.csv)</small></h3><p>Example row: subject <b>{sr[sc]}</b>"
                 + (f" ({sr[GROUP_COL]})" if GROUP_COL in s.columns else "") + ". Columns are calculated as for cell "
                 "types, using all cells of that subject.</p>" + _rows(items))

        hl = pd.read_csv(c3 / "subjects_expressing_by_threshold.csv")
        e = hl[(hl.criterion.str.startswith("chosen")) & (hl.group == "all subjects")].iloc[0]
        h.append("<h3>4. Subjects expressing <small>(subjects_expressing_headline / overview / by_threshold.csv)"
                 "</small></h3>" + _rows([
                     ("criterion", "The rule used to call a subject 'expressing'.",
                      "The chosen rule, plus alternatives (≥1, ≥3 or ≥10 cells; ≥1, ≥5 or ≥10%) to show sensitivity.",
                      f"{e.criterion}"),
                     ("group", "Which subjects are counted.", f"All subjects, then each {GROUP_COL} group separately.",
                      "all subjects"),
                     ("n_subjects_expressing", "Subjects meeting the criterion.",
                      "Count of subjects whose row in subject_summary meets the criterion.",
                      f"<b>{_f(e.n_subjects_expressing)}</b>"),
                     ("n_subjects_total", "Subjects in the group.", "Count of subjects in the group.",
                      f"<b>{_f(e.n_subjects_total)}</b>"),
                     ("proportion_expressing / pct_expressing", "Fraction of subjects expressing.",
                      "n_subjects_expressing ÷ n_subjects_total (× 100 for %)",
                      f"{_f(e.n_subjects_expressing)} ÷ {_f(e.n_subjects_total)} = {_f(e.proportion_expressing)} "
                      f": <b>{_f(e.pct_expressing)}%</b>"),
                     ("summary", "Ready-to-quote text.", "'expressing/total (pct%)'", f"<b>{e.summary}</b>")]))
        cp = c3 / "celltype_n_subjects_expressing.csv"
        if cp.exists():
            c = pd.read_csv(cp).iloc[0]
            h.append(f"<h3>5. Subjects expressing per cell type <small>(celltype_n_subjects_expressing.csv)</small></h3>"
                     f"<p>Example row: <b>{c.iloc[0]}</b>.</p>" + _rows([
                         ("n_subjects_with_celltype", "Subjects that have enough cells of this type to judge.",
                          f"Count of subjects with ≥ {p['min_group']} cells of this type.",
                          f"<b>{_f(c.n_subjects_with_celltype)}</b>"),
                         ("n_subjects_expressing", "Of those, subjects expressing the gene in this cell type.",
                          f"Within this cell type: ≥ {p['min_cells']} expressing cells and ≥ {p['min_pct']}% of the "
                          "subject's cells of this type.", f"<b>{_f(c.n_subjects_expressing)}</b>"),
                         ("n_subjects_not_expressing", "The rest.", "n_subjects_with_celltype − n_subjects_expressing",
                          f"{_f(c.n_subjects_with_celltype)} − {_f(c.n_subjects_expressing)} = "
                          f"<b>{_f(c.n_subjects_not_expressing)}</b>"),
                         ("proportion / pct_subjects_expressing", "Fraction of eligible subjects expressing.",
                          "n_subjects_expressing ÷ n_subjects_with_celltype (× 100 for %)",
                          f"{_f(c.n_subjects_expressing)} ÷ {_f(c.n_subjects_with_celltype)}: "
                          f"<b>{_f(c.pct_subjects_expressing)}%</b>")]))
        h.append("<h3>6. Subject × cell type tables <small>(subject_x_celltype_*.csv, heatmap)</small></h3><ul>"
                 "<li><code>n_cells</code>: cells of that type in that subject.</li>"
                 f"<li><code>n_expressing</code>: of those, cells with ≥1 count of {gene}.</li>"
                 "<li><code>pct_expressing</code>: n_expressing ÷ n_cells × 100. Left blank (grey in the heatmap) "
                 f"when the subject has fewer than {p['min_group']} cells of that type, because a percentage of very "
                 "few cells is unreliable.</li>"
                 "<li><code>mean_expr_expressing</code>: average expression among that subject's expressing cells "
                 "of that type.</li></ul>")
    if (out / "4_DE").exists() and RUN_INFO.get("de"):
        d = RUN_INFO["de"]
        h.append(de_calc_html(out / "4_DE", gene, d))
    h.append("<h3>7. Per-cell table <small>(per_cell_expression.csv)</small></h3><ul>"
             "<li>One row per cell, with its CELLxGENE annotations.</li>"
             f"<li><code>expression</code>: the cell's expression value for {gene}.</li>"
             f"<li><code>detected</code>: True if the cell has ≥1 count of {gene}.</li>"
             "<li><code>emb_1</code>, <code>emb_2</code>: the authors' UMAP coordinates.</li></ul>"
             "<p>Every summary table can be rebuilt from this file: group by a column, then count and average.</p>")
    h.append("<h3>8. How the plots are made</h3><ul>"
             "<li><b>UMAPs</b>: the authors' coordinates from CELLxGENE (not recomputed). Expression plots show "
             "non-expressing cells in grey and colour expressing cells by expression level, capped at the 99th "
             "percentile so a few extreme cells don't wash out the colour scale.</li>"
             "<li><b>Bar plots</b>: pct_expressing (or counts) per group, labelled 'expressing / total'.</li>"
             "<li><b>Violin plots</b>: expression values of expressing cells only. Groups with fewer than 3 "
             "expressing cells are not drawn.</li>"
             "<li><b>Dot plots</b>: dot size = pct_expressing; colour = mean_expr_expressing_cells.</li>"
             "<li><b>Subject box plot</b>: one dot per subject (that subject's pct_expressing). The p-value compares "
             "disease groups using subjects as replicates: two-sided Mann-Whitney U for 2 groups, Kruskal-Wallis for "
             "more than 2.</li>"
             "<li><b>Subject heatmap</b>: pct_expressing per subject × cell type; grey = too few cells.</li></ul>"
             "<p><small>All numbers in the CSV files are rounded (counts exact; proportions 5 decimals; "
             "percentages 3 decimals; means 4 decimals), so recomputed values can differ in the last digit.</small>"
             "</p></div>")
    return "\n".join(h)


CALC_CSS = """
.calcpage{background:#fff;color:#222;max-width:1100px;padding:20px 28px;border:1px solid #e5e5e5;
 font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;line-height:1.45}
.calcpage h2{margin-top:0} .calcpage h3{margin-top:26px;border-bottom:1px solid #eee;padding-bottom:3px}
table.calc{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}
table.calc th{background:#f4f6f8;text-align:left} table.calc td,table.calc th{border:1px solid #ddd;padding:5px 8px;
 vertical-align:top} table.calc td:first-child{white-space:nowrap}
"""


def write_calc_page(out, gene):
    p = Path(out) / "HOW_IT_IS_CALCULATED.html"
    p.write_text(f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{gene}: how it is calculated</title>"
                 f"<style>body{{background:#fff;margin:20px}}{CALC_CSS}</style></head><body>"
                 f"{calc_guide_html(out, gene)}</body></html>")
    return p


# ----------------------------------------------------------------------------- descriptive report (on demand)
REPORT_SYSTEM = """You are an experienced single-cell genomics analyst writing a descriptive results report for a \
scientist. You receive summary statistics from a single-gene analysis of a public CELLxGENE dataset.
Rules:
- Use ONLY the numbers provided. Never invent values, cell types, groups or p-values. Quote numbers as given, with \
denominators (e.g. "48 of 113,304 cells, 0.04%"; "5/6 patients, 83.3%").
- You may add brief, well-established background on the gene's known biology, explicitly marked as background.
- Be careful with uncertainty: detection means at least one count; single-cell data has dropout; very low \
percentages (well under 1%) may reflect ambient RNA or doublets; few subjects means low power; p-values are \
exploratory. Do not over-claim and do not give medical advice.
- Write in Markdown, formal scientific prose, past tense for results."""

SECTION_ASK = {
    "1_cellular": "Write the 'Cellular level' results subsection. Describe overall detection, which cell types "
                  "express the gene most and least (proportion and level among expressing cells), and differences by "
                  "group/tissue/sex if given, referring to the figures (UMAP, bar plot, dot plot, violin plot). "
                  "Then add a paragraph starting with '**Interpretation.**'. 150-300 words, no headings.",
    "2_cell_specific": "Write the 'Cell-specific level' results subsection for the selected population(s): cells, "
                       "detection, differences by condition, and how many subjects express it in that population. "
                       "Then add a paragraph starting with '**Interpretation.**'. 120-250 words, no headings.",
    "3_sample": "Write the 'Sample (subject) level' results subsection: number and proportion of subjects expressing "
                "the gene overall and per group (grouping column in context), how much this depends on the threshold (sensitivity table), "
                "which cell types are expressed across most subjects, and any subject-level group comparison. Then add a "
                "paragraph starting with '**Interpretation.**'. 150-300 words, no headings.",
}

SECTION_ASK["4_de"] = ("Write the 'Differential expression' results subsection: for each cell type, the NEBULA "
                       "comparison of gene-positive vs gene-negative cells (number of genes up/down, top genes; note "
                       "skipped cell types), and the Fisher's exact results for cell proportions (cell type vs rest; "
                       "between groups). Mention that gene-level associations are co-expression patterns, not "
                       "causal effects, and that Fisher tests treat cells as independent. Then add a paragraph "
                       "starting with '**Interpretation.**'. 150-300 words, no headings.")

SUMMARY_ASK = """Below are the results subsections of the report. Write two parts in Markdown:
1. '**Summary.**' followed by one paragraph (4-6 sentences) giving the main findings with key numbers.
2. '**Conclusions and limitations.**' followed by 3-5 bullets: what can and cannot be concluded from this \
analysis.
Use only numbers that appear in the subsections. No other headings."""


def generate_descriptive_report(res, model=None, url=None, progress=print):
    out, gene, digest = Path(res["folder"]), res["symbol"], res["digest"]
    texts = {}
    for key in ["1_cellular", "2_cell_specific", "3_sample", "4_de"]:
        if not digest.get(key):
            continue
        progress(f"writing the {key.split('_', 1)[1].replace('_', '-')} section")
        payload = {"context": digest["context"], key: digest[key]}
        texts[key] = ask_llm([{"role": "system", "content": REPORT_SYSTEM},
                              {"role": "user", "content": SECTION_ASK[key] + "\n\nData (JSON):\n" +
                               json.dumps(payload, default=str)}], model, url)
    progress("writing the summary and conclusions")
    joined = "\n\n".join(f"### {k}\n{t}" for k, t in texts.items())
    summary = ask_llm([{"role": "system", "content": REPORT_SYSTEM},
                       {"role": "user", "content": SUMMARY_ASK + "\n\nContext: " +
                        json.dumps(digest["context"], default=str) + "\n\n" + joined}], model, url)

    titles = {"1_cellular": "Cellular level (all cell types)", "2_cell_specific": "Cell-specific level",
              "3_sample": "Sample level (subjects / patients)", "4_de": "Differential expression and cell proportions"}
    figs = {"1_cellular": ["1_cellular_level/umap_cell_types.png", "1_cellular_level/umap_expression.png",
                           "1_cellular_level/bar_pct_expressing_by_celltype.png",
                           f"1_cellular_level/dotplot_celltype_by_{safe_name(GROUP_COL)}.png", "1_cellular_level/dotplot_celltype.png"],
            "2_cell_specific": [str(p.relative_to(out)) for d in sorted((out / "2_cell_specific_level").glob("*/"))
                                for p in [d / "umap_population_expression.png", d / f"bar_pct_expressing_by_{GROUP_COL}.png"]
                                if p.exists()] if (out / "2_cell_specific_level").exists() else [],
            "3_sample": ["3_sample_level/bar_n_and_proportion_subjects_expressing.png",
                         "3_sample_level/heatmap_subject_x_celltype.png",
                         f"3_sample_level/box_subject_pct_by_{safe_name(GROUP_COL)}.png"],
            "4_de": ["4_DE/heatmap_top_de_genes.png", "4_DE/fisher_between_groups.png"]}
    tables = {"1_cellular": ("1_cellular_level/celltype_summary.csv", "Expression per cell type (top 15)"),
              "3_sample": ("3_sample_level/subjects_expressing_overview.csv",
                           "Subjects expressing the gene under different thresholds")}
    r = RUN_INFO
    md = [f"# {gene} expression in {r.get('dataset')}", "",
          f"_Descriptive report generated with {model or LLM_MODEL} from the analysis summary statistics on "
          f"{time.strftime('%Y-%m-%d')}. AI-written text: verify against the tables and figures._", "",
          summary, "", "## Methods", "", methods_markdown(), "", "## Results", ""]
    html = [f"<h1>{gene} expression in {r.get('dataset')}</h1>",
            f"<p class='sub'>Collection: {r.get('collection')} · DOI: {r.get('doi')} · {r.get('n_cells', 0):,} cells · "
            f"{r.get('n_subjects')} subjects</p>",
            f"<p class='ai'>Text written by {model or LLM_MODEL} from the computed summary statistics. "
            "Verify against the tables, figures and the 'How it's calculated' page.</p>",
            md_to_html(summary), "<h2>Methods</h2>", md_to_html(methods_markdown()), "<h2>Results</h2>"]
    for key, t in texts.items():
        md += [f"### {titles[key]}", "", t, ""]
        html += [f"<h3>{titles[key]}</h3>", md_to_html(t)]
        if key in tables and (out / tables[key][0]).exists():
            tdf = pd.read_csv(out / tables[key][0]).head(15)
            html.append(f"<p class='cap'><b>Table.</b> {tables[key][1]}.</p>"
                        f"<div class='tbl'>{tdf.to_html(index=False, classes='df', border=0)}</div>")
        fl = [f for f in figs.get(key, []) if (out / f).exists()][:4]
        if fl:
            html.append('<div class="grid">' + "".join(
                f'<figure><img src="data:image/png;base64,{img_b64(out / f)}"><figcaption>{f}</figcaption></figure>'
                for f in fl) + "</div>")
    css = ("body{font-family:Georgia,'Times New Roman',serif;max-width:980px;margin:auto;padding:28px;color:#222;"
           "line-height:1.55;background:#fff} h1{font-size:26px} h2{border-bottom:2px solid #ddd;margin-top:34px}"
           ".sub{color:#555}.ai{background:#fbf8ef;border-left:4px solid #c9a227;padding:6px 10px;font-size:14px}"
           ".grid{display:flex;flex-wrap:wrap;gap:14px}figure{margin:0;flex:1 1 430px}figure img{width:100%;"
           "border:1px solid #eee}figcaption,.cap{font-size:12px;color:#555;font-family:Helvetica,Arial,sans-serif}"
           "table.df{border-collapse:collapse;font-size:12px;font-family:Helvetica,Arial,sans-serif}"
           "table.df td,table.df th{border:1px solid #ddd;padding:3px 6px}.tbl{overflow-x:auto}")
    body = "\n".join(html)
    (out / "DESCRIPTIVE_REPORT.md").write_text("\n".join(md))
    (out / "DESCRIPTIVE_REPORT.html").write_text(f"<!DOCTYPE html><html><head><meta charset='utf-8'>"
                                                 f"<title>{gene} descriptive report</title><style>{css}</style></head>"
                                                 f"<body>{body}</body></html>")
    return body, css


# ----------------------------------------------------------------------------- plotting helpers
def palette(n):
    if n <= 10:
        base = list(plt.get_cmap("tab10").colors)
    elif n <= 20:
        base = list(plt.get_cmap("tab20").colors)
    elif n <= 60:
        base = (list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
                + list(plt.get_cmap("tab20c").colors))
    else:
        base = [plt.get_cmap("gist_rainbow")(i / n)[:3] for i in range(n)]
    return [tuple(c[:3]) for c in base[:n]]


def point_size(n):
    return float(np.clip(60000 / max(n, 1), 0.1, 12))


def style_embedding(ax, title, key):
    lab = (key or "").replace("X_", "").upper() or "EMB"
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(f"{lab} 1")
    ax.set_ylabel(f"{lab} 2")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)


def savefig(fig, path, section, caption):
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    FIGS.append((section, caption, Path(path)))


def save_table(df, path, section, caption, index=True, preview_rows=25):
    df.to_csv(path, index=index)
    TABLES.append((section, caption, Path(path), df.head(preview_rows)))


def plot_embedding_categorical(xy, labels, path, title, key, section, caption, max_legend=40, order=None):
    labels = pd.Series(labels).astype(str).values
    cats = order or pd.Series(labels).value_counts().index.tolist()
    cols = np.array(palette(len(cats)))
    codes = pd.Categorical(labels, categories=cats).codes
    perm = np.random.default_rng(0).permutation(len(labels))
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(xy[perm, 0], xy[perm, 1], c=cols[codes[perm]], s=point_size(len(labels)), lw=0, rasterized=True)
    if len(cats) <= max_legend:
        handles = [Line2D([], [], marker="o", ls="", color=cols[i], label=c) for i, c in enumerate(cats)]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
                  fontsize=8 if len(cats) < 20 else 6, markerscale=1.2, ncol=1 if len(cats) <= 25 else 2)
    else:
        for i, c in enumerate(cats):
            m = codes == i
            if m.sum():
                ax.text(*np.median(xy[m], axis=0), c, fontsize=5, ha="center")
    style_embedding(ax, title, key)
    savefig(fig, path, section, caption)


def plot_embedding_expression(xy, expr, detected, path, title, key, section, caption, background_xy=None):
    fig, ax = plt.subplots(figsize=(7.5, 7))
    if background_xy is not None:
        ax.scatter(background_xy[:, 0], background_xy[:, 1], c="#eeeeee", s=point_size(len(background_xy)),
                   lw=0, rasterized=True)
    s = point_size(len(expr))
    neg = ~detected
    ax.scatter(xy[neg, 0], xy[neg, 1], c="#d3d3d3", s=s, lw=0, rasterized=True)
    order = np.argsort(expr[detected])
    pos_xy, pos_e = xy[detected][order], expr[detected][order]
    vmax = float(np.percentile(pos_e, 99)) if pos_e.size else 1.0
    sc = ax.scatter(pos_xy[:, 0], pos_xy[:, 1], c=pos_e, cmap="viridis", s=s * 1.5, lw=0,
                    vmin=0, vmax=max(vmax, 1e-6), rasterized=True)
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("expression (expressing cells)")
    style_embedding(ax, title, key)
    savefig(fig, path, section, caption)


def plot_embedding_binary(xy, detected, path, title, key, section, caption):
    fig, ax = plt.subplots(figsize=(7.5, 7))
    s = point_size(len(detected))
    ax.scatter(xy[~detected, 0], xy[~detected, 1], c="#d3d3d3", s=s, lw=0, rasterized=True)
    ax.scatter(xy[detected, 0], xy[detected, 1], c="#d62728", s=s * 1.5, lw=0, rasterized=True)
    n, k = len(detected), int(detected.sum())
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#d62728", label=f"expressing ({k:,})"),
                       Line2D([], [], marker="o", ls="", color="#d3d3d3", label=f"not detected ({n - k:,})")],
              loc="upper right", frameon=False, fontsize=8)
    style_embedding(ax, title, key)
    savefig(fig, path, section, caption)


def plot_embedding_highlight(xy, mask, path, title, key, section, caption):
    fig, ax = plt.subplots(figsize=(7.5, 7))
    s = point_size(len(mask))
    ax.scatter(xy[~mask, 0], xy[~mask, 1], c="#e0e0e0", s=s, lw=0, rasterized=True)
    ax.scatter(xy[mask, 0], xy[mask, 1], c="#1f77b4", s=s * 1.5, lw=0, rasterized=True)
    style_embedding(ax, title, key)
    savefig(fig, path, section, caption)


def plot_barh(values, path, xlabel, title, section, caption, colors=None, annotate=None, max_items=80):
    values = values.iloc[:max_items]
    h = max(2.5, 0.28 * len(values) + 1)
    fig, ax = plt.subplots(figsize=(8, h))
    y = np.arange(len(values))
    ax.barh(y, values.values, color=colors if colors is not None else "#4c72b0")
    ax.set_yticks(y)
    ax.set_yticklabels([str(i)[:60] for i in values.index], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=11)
    if annotate is not None:
        annotate = annotate.iloc[:max_items]
        xmax = values.max() if len(values) else 1
        for yi, (v, t) in enumerate(zip(values.values, annotate.values)):
            ax.text(v + 0.01 * (xmax or 1), yi, t, va="center", fontsize=6.5)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    savefig(fig, path, section, caption)


def plot_violin(df, group, order, path, title, section, caption, min_cells=3, max_items=60):
    order = [g for g in order if g in set(df[group])][:max_items]
    data, labels = [], []
    for g in order:
        v = df.loc[(df[group] == g) & df["detected"], "expression"].values
        if len(v) >= min_cells:
            if np.ptp(v) == 0:
                v = v + np.random.default_rng(0).normal(0, 1e-3, len(v))
            data.append(v)
            labels.append(g)
    if not data:
        return
    w = max(6, 0.35 * len(data) + 2)
    fig, ax = plt.subplots(figsize=(w, 4.8))
    parts = ax.violinplot(data, showmedians=True, widths=0.8)
    for pc in parts["bodies"]:
        pc.set_facecolor("#8172b2")
        pc.set_alpha(0.7)
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels([l[:45] for l in labels], rotation=75, ha="right", fontsize=7)
    ax.set_ylabel("expression (expressing cells only)")
    ax.set_title(title, fontsize=11)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    savefig(fig, path, section, caption)


def plot_dotplot(df, row, col, row_order, path, title, section, caption, max_rows=80):
    row_order = row_order[:max_rows]
    col_order = df[col].value_counts().index.tolist()
    g = df.groupby([row, col], observed=True)
    pct = (100 * g["detected"].mean()).unstack().reindex(index=row_order, columns=col_order)
    mean = df[df.detected].groupby([row, col], observed=True)["expression"].mean().unstack()
    mean = mean.reindex(index=row_order, columns=col_order)
    yy, xx = np.meshgrid(np.arange(len(row_order)), np.arange(len(col_order)), indexing="ij")
    fig, ax = plt.subplots(figsize=(max(4, 1.1 * len(col_order) + 3.5), max(3, 0.28 * len(row_order) + 1.5)))
    sizes = np.nan_to_num(pct.values) / 100 * 250
    sc = ax.scatter(xx.ravel(), yy.ravel(), s=sizes.ravel(), c=np.nan_to_num(mean.values).ravel(),
                    cmap="Reds", edgecolors="#555555", linewidths=0.3)
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels([c[:30] for c in col_order], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels([r[:55] for r in row_order], fontsize=7)
    ax.set_xlim(-0.6, len(col_order) - 0.4)
    ax.invert_yaxis()
    cb = fig.colorbar(sc, ax=ax, shrink=0.4, pad=0.02)
    cb.set_label("mean expr. in expressing cells")
    for p in [10, 50, 100]:
        ax.scatter([], [], s=p / 100 * 250, c="white", edgecolors="#555555", label=f"{p}%")
    ax.legend(title="% expressing", loc="upper left", bbox_to_anchor=(1.25, 1), frameon=False, fontsize=7)
    ax.set_title(title, fontsize=11)
    savefig(fig, path, section, caption)


def plot_heatmap(mat, path, title, section, caption, cbar_label="% expressing"):
    rows, cols = mat.shape
    fig, ax = plt.subplots(figsize=(min(30, max(6, 0.3 * cols + 4)), min(40, max(4, 0.25 * rows + 2))))
    cmap = plt.get_cmap("magma_r").copy()
    cmap.set_bad("#f0f0f0")
    im = ax.imshow(np.ma.masked_invalid(mat.values.astype(float)), aspect="auto", cmap=cmap)
    ax.set_xticks(range(cols))
    ax.set_xticklabels([str(c)[:40] for c in mat.columns], rotation=75, ha="right", fontsize=6)
    ax.set_yticks(range(rows))
    ax.set_yticklabels([str(r)[:45] for r in mat.index], fontsize=6)
    cb = fig.colorbar(im, ax=ax, shrink=0.5)
    cb.set_label(cbar_label)
    ax.set_title(title + "\n(grey = too few cells)", fontsize=10)
    savefig(fig, path, section, caption)


# ----------------------------------------------------------------------------- summaries
def summarize(df, by):
    g = df.groupby(by, observed=True)
    out = pd.DataFrame({"n_cells": g.size(), "n_expressing": g["detected"].sum().astype(int)})
    out["n_not_expressing"] = out["n_cells"] - out["n_expressing"]
    out["proportion_expressing"] = (out["n_expressing"] / out["n_cells"]).round(5)
    out["pct_expressing"] = (100 * out["proportion_expressing"]).round(3)
    out["mean_expr_all_cells"] = g["expression"].mean().round(4)
    pos = df[df["detected"]].groupby(by, observed=True)["expression"]
    out["mean_expr_expressing_cells"] = pos.mean().reindex(out.index).fillna(0).round(4)
    out["median_expr_expressing_cells"] = pos.median().reindex(out.index).fillna(0).round(4)
    tot = max(int(df["detected"].sum()), 1)
    out["pct_of_all_expressing_cells"] = (100 * out["n_expressing"] / tot).round(3)
    return out.sort_values(["pct_expressing", "n_expressing"], ascending=False)


def overall_row(df, label):
    n, k = len(df), int(df["detected"].sum())
    return pd.DataFrame([{
        "population": label, "n_cells_total": n, "n_cells_expressing": k, "n_cells_not_expressing": n - k,
        "proportion_expressing": round(k / max(n, 1), 5), "pct_expressing": round(100 * k / max(n, 1), 3),
        "mean_expr_all_cells": round(float(df["expression"].mean()), 4) if n else 0,
        "mean_expr_expressing_cells": round(float(df.loc[df.detected, "expression"].mean()), 4) if k else 0,
    }])


def covariates_present(df, candidates=("disease", "tissue", "sex", "development_stage", "assay",
                                       "self_reported_ethnicity", "tissue_type")):
    cands = list(dict.fromkeys([GROUP_COL] + list(candidates)))
    return [c for c in cands if c in df.columns and df[c].nunique() > 1]


# ----------------------------------------------------------------------------- level 1: cellular
def level_cellular(df, xy, key, gene, ct_col, out):
    sec = "1. Cellular level (all cell types)"
    out.mkdir(parents=True, exist_ok=True)
    header(sec)

    ov = overall_row(df, "all cells")
    ov["n_cell_types"] = df[ct_col].nunique()
    save_table(ov, out / "overall_summary.csv", sec, "Overall summary", index=False)
    r = ov.iloc[0]
    say(f"   Total cells: {r.n_cells_total:,} | expressing {gene}: {r.n_cells_expressing:,} ({r.pct_expressing:.2f}%)")

    ct = summarize(df, ct_col)
    ct.index.name = ct_col
    save_table(ct, out / "celltype_summary.csv", sec, f"{gene} per cell type")
    order = ct.index.tolist()

    if xy is not None:
        plot_embedding_categorical(xy, df[ct_col], out / "umap_cell_types.png", f"Cell types ({ct_col})", key,
                                   sec, "Embedding coloured by cell type")
        plot_embedding_expression(xy, df["expression"].values, df["detected"].values, out / "umap_expression.png",
                                  f"{gene} expression", key, sec, f"{gene} expression on the embedding")
        plot_embedding_binary(xy, df["detected"].values, out / "umap_expressing_vs_not.png",
                              f"{gene}+ vs {gene}- cells", key, sec, "Expressing vs non-expressing cells")

    ann = ct.apply(lambda x: f"{int(x.n_expressing):,}/{int(x.n_cells):,}", axis=1)
    plot_barh(ct["pct_expressing"], out / "bar_pct_expressing_by_celltype.png", f"% cells expressing {gene}",
              f"% {gene}+ cells per cell type", sec, "% expressing per cell type (label: expressing/total)",
              annotate=ann)
    nsort = ct.sort_values("n_cells", ascending=False)
    fig_df = nsort[["n_expressing", "n_not_expressing"]]
    h = max(2.5, 0.28 * len(fig_df) + 1)
    fig, ax = plt.subplots(figsize=(8, h))
    y = np.arange(len(fig_df))
    ax.barh(y, fig_df["n_expressing"], color="#d62728", label=f"{gene}+")
    ax.barh(y, fig_df["n_not_expressing"], left=fig_df["n_expressing"], color="#c7c7c7", label=f"{gene}-")
    ax.set_yticks(y)
    ax.set_yticklabels([str(i)[:60] for i in fig_df.index], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("number of cells")
    ax.legend(frameon=False)
    ax.set_title(f"Cell counts per cell type ({gene}+ / {gene}-)", fontsize=11)
    savefig(fig, out / "bar_cell_counts_by_celltype.png", sec, "Number of expressing / non-expressing cells")
    plot_barh(ct["mean_expr_expressing_cells"].sort_values(ascending=False), out / "bar_mean_expression_by_celltype.png",
              "mean expression in expressing cells", f"{gene} level per cell type", sec,
              "Mean expression among expressing cells", colors="#55a868")
    plot_violin(df, ct_col, order, out / "violin_by_celltype.png", f"{gene} in expressing cells", sec,
                "Distribution of expression among expressing cells")
    dot_col = GROUP_COL if GROUP_COL in df.columns and df[GROUP_COL].nunique() > 1 else None
    if dot_col:
        plot_dotplot(df, ct_col, dot_col, order, out / f"dotplot_celltype_by_{safe_name(dot_col)}.png",
                     f"{gene}: cell type x {dot_col}", sec, f"Dot plot (size = % expressing) by {dot_col}")
    else:
        tmp = df.assign(_gene=gene)
        plot_dotplot(tmp, ct_col, "_gene", order, out / "dotplot_celltype.png", f"{gene} per cell type", sec,
                     "Dot plot (size = % expressing)")

    for cov in covariates_present(df):
        s = summarize(df, cov)
        s.index.name = cov
        save_table(s, out / f"summary_by_{cov}.csv", sec, f"{gene} by {cov}")
        if s.shape[0] <= 60:
            plot_barh(s["pct_expressing"], out / f"bar_pct_expressing_by_{cov}.png", f"% cells expressing {gene}",
                      f"% {gene}+ cells by {cov}", sec, f"% expressing by {cov}", colors="#dd8452",
                      annotate=s.apply(lambda x: f"{int(x.n_expressing):,}/{int(x.n_cells):,}", axis=1))

    cols = ["cell_id", ct_col] + [c for c in dict.fromkeys(["donor_id", GROUP_COL, "disease", "tissue", "sex"])
                                  if c in df.columns and c != ct_col]
    cols += ["expression", "detected"] + (["emb_1", "emb_2"] if "emb_1" in df.columns else [])
    df[cols].to_csv(out / "per_cell_expression.csv", index=False)
    NOTES.append((sec, f"Per-cell table: <code>{(out / 'per_cell_expression.csv').name}</code> ({len(df):,} rows). "
                       f"Columns: expression level, detected (True if the gene was detected in that cell), embedding coordinates."))
    return ct


# ----------------------------------------------------------------------------- level 2: cell-specific
def pick_populations(ct, preset=None):
    say("\nCell types (sorted by % expressing):")
    for i, (name, r) in enumerate(ct.iterrows(), 1):
        say(f"  [{i:>3}] {name[:60]:<60} {int(r.n_cells):>9,} cells  {r.pct_expressing:6.2f}% +")
    raw = preset or ask("Choose population(s): numbers like 1,4,7  or names separated by ';'  ('none' to skip)",
                        "1")
    if raw is None or str(raw).lower() in ("none", "skip", ""):
        return []
    names = ct.index.tolist()
    tokens = [t.strip() for t in re.split(r"[,\s]+", raw) if t.strip()]
    if tokens and all(t.isdigit() for t in tokens):
        return [names[int(t) - 1] for t in tokens if 0 < int(t) <= len(names)]
    chosen = []
    for t in [t.strip() for t in re.split(r"[;|]", raw) if t.strip()]:
        exact = [n for n in names if n.lower() == t.lower()]
        sub = exact or [n for n in names if t.lower() in n.lower()]
        if not sub:
            say(f"    '{t}' not found; skipping.")
        chosen += [n for n in sub if n not in chosen]
    return chosen


def level_cell_specific(df, xy, key, gene, ct_col, populations, donor_col, out, min_cells, min_pct, min_group):
    PARAMS.update(min_cells=min_cells, min_pct=min_pct, min_group=min_group)
    sec = "2. Cell-specific level"
    header(sec)
    groups = [(p, [p]) for p in populations]
    if len(populations) > 1:
        groups.append(("combined_selection", populations))
    results = []
    for label, members in groups:
        sub_sec = f"{sec}: {label}"
        d = out / safe_name(label)
        d.mkdir(parents=True, exist_ok=True)
        mask = df[ct_col].isin(members).values
        sub = df[mask]
        ov = overall_row(sub, label)
        ov["members"] = "; ".join(members)
        results.append(ov)
        save_table(ov, d / "population_summary.csv", sub_sec, "Population summary", index=False)
        r = ov.iloc[0]
        say(f"   {label}: {r.n_cells_total:,} cells, {r.n_cells_expressing:,} express {gene} ({r.pct_expressing:.2f}%)")
        if xy is not None:
            plot_embedding_highlight(xy, mask, d / "umap_population_highlight.png", f"{label} on full embedding",
                                     key, sub_sec, "Population highlighted on the full embedding")
            plot_embedding_expression(xy[mask], sub["expression"].values, sub["detected"].values,
                                      d / "umap_population_expression.png", f"{gene} in {label}", key, sub_sec,
                                      f"{gene} expression within the population (zoomed)")
            plot_embedding_binary(xy[mask], sub["detected"].values, d / "umap_population_expressing_vs_not.png",
                                  f"{gene}+ vs {gene}- in {label}", key, sub_sec, "Expressing vs non-expressing")
        if len(members) > 1:
            s = summarize(sub, ct_col)
            save_table(s, d / "summary_by_member_celltype.csv", sub_sec, "Breakdown by member cell type")
        for cov in covariates_present(sub):
            s = summarize(sub, cov)
            s.index.name = cov
            save_table(s, d / f"summary_by_{cov}.csv", sub_sec, f"{gene} in {label} by {cov}")
            if s.shape[0] <= 60:
                plot_barh(s["pct_expressing"], d / f"bar_pct_expressing_by_{cov}.png", f"% expressing {gene}",
                          f"{label}: % {gene}+ by {cov}", sub_sec, f"% expressing by {cov}", colors="#dd8452",
                          annotate=s.apply(lambda x: f"{int(x.n_expressing):,}/{int(x.n_cells):,}", axis=1))
                plot_violin(sub, cov, s.index.tolist(), d / f"violin_by_{cov}.png",
                            f"{label}: {gene} by {cov} (expressing cells)", sub_sec, f"Expression by {cov}")
        sub.to_csv(d / "per_cell_expression.csv", index=False)
        if donor_col:
            level_sample(sub, gene, ct_col, donor_col, d / "sample_level", min_cells, min_pct, min_group,
                         section=f"{sub_sec} - per subject", per_celltype=len(members) > 1)
    if results:
        allp = pd.concat(results, ignore_index=True)
        save_table(allp, out / "selected_populations_summary.csv", sec, "All selected populations", index=False)


# ----------------------------------------------------------------------------- level 3: sample / subject
def level_sample(df, gene, ct_col, donor_col, out, min_cells, min_pct, min_group,
                 section="3. Sample level (subjects / patients)", per_celltype=True):
    PARAMS.update(min_cells=min_cells, min_pct=min_pct, min_group=min_group)
    out.mkdir(parents=True, exist_ok=True)
    if section.startswith("3."):
        header(section)
    has_dis = GROUP_COL in df.columns and df[GROUP_COL].nunique() > 1
    unit_cols = [donor_col] + ([GROUP_COL] if has_dis and GROUP_COL != donor_col else [])

    ds = summarize(df, unit_cols).reset_index()
    ds["expressing_subject"] = (ds["n_expressing"] >= min_cells) & (ds["pct_expressing"] >= min_pct)
    for c in ["sex", "tissue", "development_stage"]:
        if c in df.columns and c not in unit_cols:
            lut = df.groupby(unit_cols, observed=True)[c].agg(lambda s: "; ".join(sorted(set(s))))
            keys = ds[unit_cols[0]] if len(unit_cols) == 1 else ds[unit_cols].apply(tuple, axis=1)
            ds[c] = keys.map(lut).values
    ds = ds.sort_values("pct_expressing", ascending=False)
    save_table(ds, out / "subject_summary.csv", section, f"{gene} per subject", index=False)

    n_sub = ds[donor_col].nunique()
    n_pos = ds.loc[ds["expressing_subject"], donor_col].nunique()
    say(f"   {n_pos} of {n_sub} subjects express {gene} ({100 * n_pos / max(n_sub, 1):.1f}%) "
        f"[rule: >= {min_cells} expressing cells and >= {min_pct}% of cells]")

    criteria = [(">=1 expressing cell", ds.n_expressing >= 1),
                (">=3 expressing cells", ds.n_expressing >= 3),
                (">=10 expressing cells", ds.n_expressing >= 10),
                (">=1% cells expressing", ds.pct_expressing >= 1),
                (">=5% cells expressing", ds.pct_expressing >= 5),
                (">=10% cells expressing", ds.pct_expressing >= 10),
                (f"chosen rule (>={min_cells} cells & >={min_pct}%)", ds.expressing_subject)]
    groups = [("all subjects", ds)] + ([(d, g) for d, g in ds.groupby(GROUP_COL)] if has_dis else [])
    long_rows = []
    for label, m in criteria:
        for gname, g in groups:
            tot = g[donor_col].nunique()
            k = g.loc[m[g.index], donor_col].nunique()
            long_rows.append({"criterion": label, "group": gname, "n_subjects_expressing": k,
                              "n_subjects_not_expressing": tot - k, "n_subjects_total": tot,
                              "proportion_expressing": round(k / max(tot, 1), 4),
                              "pct_expressing": round(100 * k / max(tot, 1), 2),
                              "summary": f"{k}/{tot} ({100 * k / max(tot, 1):.1f}%)"})
    thr_long = pd.DataFrame(long_rows)
    save_table(thr_long, out / "subjects_expressing_by_threshold.csv", section,
               "Number and proportion of subjects expressing the gene (several thresholds, per group)", index=False)
    wide = thr_long.pivot(index="criterion", columns="group", values="summary")
    wide = wide.reindex(index=[c for c, _ in criteria], columns=[g for g, _ in groups])
    save_table(wide, out / "subjects_expressing_overview.csv", section,
               "Overview: expressing / total subjects (proportion)")

    headline = thr_long[thr_long.criterion.str.startswith("chosen rule")].copy()
    save_table(headline.drop(columns="criterion"), out / "subjects_expressing_headline.csv", section,
               "Headline: subjects expressing under the chosen rule", index=False)
    note = (f"<b>{n_pos} of {n_sub} subjects ({100 * n_pos / max(n_sub, 1):.1f}%)</b> express {gene} "
            f"(&ge;{min_cells} expressing cells and &ge;{min_pct}% of that subject's cells)")
    if has_dis:
        note += ": " + "; ".join(f"{r.group}: {r.summary}" for r in headline.itertuples() if r.group != "all subjects")
    NOTES.append((section, note + f". Subject unit = <code>{' + '.join(unit_cols)}</code>."))

    # stacked bar: number of expressing vs non-expressing subjects per group, labelled with proportion
    hl = headline.set_index("group")
    fig, ax = plt.subplots(figsize=(max(4, 1.3 * len(hl) + 2), 4.5))
    x = np.arange(len(hl))
    ax.bar(x, hl.n_subjects_expressing, color="#d62728", label=f"{gene}+ subjects")
    ax.bar(x, hl.n_subjects_not_expressing, bottom=hl.n_subjects_expressing, color="#c7c7c7",
           label=f"{gene}- subjects")
    for xi, r in zip(x, hl.itertuples()):
        ax.text(xi, r.n_subjects_total, f"{r.n_subjects_expressing}/{r.n_subjects_total}\n({r.pct_expressing:.1f}%)",
                ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([str(g)[:30] for g in hl.index], rotation=30, ha="right")
    ax.set_ylabel("number of subjects")
    ax.set_ylim(0, hl.n_subjects_total.max() * 1.25)
    ax.set_title(f"Subjects expressing {gene} (chosen rule)", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    savefig(fig, out / "bar_n_and_proportion_subjects_expressing.png", section,
            "Number and proportion of subjects expressing the gene")

    # proportion of subjects across thresholds (robustness)
    piv = thr_long.pivot(index="criterion", columns="group", values="pct_expressing")
    piv = piv.reindex(index=[c for c, _ in criteria], columns=[g for g, _ in groups])
    fig, ax = plt.subplots(figsize=(9, 4.5))
    w = 0.8 / piv.shape[1]
    cols_ = palette(piv.shape[1])
    for j, gname in enumerate(piv.columns):
        ax.bar(np.arange(len(piv)) + j * w - 0.4 + w / 2, piv[gname], w, label=gname, color=cols_[j])
    ax.set_xticks(np.arange(len(piv)))
    ax.set_xticklabels(piv.index, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("% of subjects expressing")
    ax.set_ylim(0, 105)
    ax.set_title(f"Proportion of subjects expressing {gene} under different thresholds", fontsize=10)
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    savefig(fig, out / "bar_proportion_subjects_by_threshold.png", section,
            "Proportion of subjects expressing the gene across thresholds")

    lab = ds[unit_cols].astype(str).agg(" | ".join, axis=1)
    vals = pd.Series(ds["pct_expressing"].values, index=lab.values)
    if has_dis:
        dis_levels = ds[GROUP_COL].unique().tolist()
        cmap = dict(zip(dis_levels, palette(len(dis_levels))))
        colors = [cmap[d] for d in ds[GROUP_COL]]
    else:
        colors = ["#d62728" if e else "#9e9e9e" for e in ds["expressing_subject"]]
    plot_barh(vals, out / "bar_pct_expressing_by_subject.png", f"% cells expressing {gene}",
              f"% {gene}+ cells per subject", section,
              "% expressing per subject" + (f" (coloured by {GROUP_COL})" if has_dis else " (red = expressing subject)"),
              colors=colors[:80], annotate=pd.Series(
                  [f"{a:,}/{b:,}" for a, b in zip(ds.n_expressing, ds.n_cells)], index=lab.values), max_items=80)

    if has_dis and ds[GROUP_COL].nunique() >= 2:
        groups = {k: g["pct_expressing"].values for k, g in ds.groupby(GROUP_COL) if len(g) >= 2}
        fig, ax = plt.subplots(figsize=(max(4, 1.3 * len(groups) + 2), 4.5))
        keys = list(groups)
        if keys:
            ax.boxplot([groups[k] for k in keys], widths=0.5, showfliers=False)
            for i, k in enumerate(keys, 1):
                jitter = np.random.default_rng(i).uniform(-0.12, 0.12, len(groups[k]))
                ax.scatter(np.full(len(groups[k]), i) + jitter, groups[k], s=18, alpha=0.8, zorder=3)
            ax.set_xticks(range(1, len(keys) + 1))
            ax.set_xticklabels([k[:30] for k in keys], rotation=30, ha="right")
            stat_txt = ""
            if len(keys) == 2:
                p = stats.mannwhitneyu(groups[keys[0]], groups[keys[1]]).pvalue
                stat_txt = f"Mann-Whitney U p = {p:.3g}"
            elif len(keys) > 2:
                p = stats.kruskal(*[groups[k] for k in keys]).pvalue
                stat_txt = f"Kruskal-Wallis p = {p:.3g}"
            ax.set_title(f"Subject-level % {gene}+ by {GROUP_COL}\n{stat_txt}", fontsize=10)
            ax.set_ylabel(f"% cells expressing {gene} (per subject)")
            for sp in ["top", "right"]:
                ax.spines[sp].set_visible(False)
            savefig(fig, out / f"box_subject_pct_by_{safe_name(GROUP_COL)}.png", section,
                    "Subject-level comparison (each dot = one subject; test on subjects, not cells)")
            if stat_txt:
                STATS.append({"level": section, "test": stat_txt.split(" p")[0], "p_value": float(f"{p:.3g}"),
                              "groups": {k: {"n_subjects": int(len(groups[k])),
                                             "median_pct_expressing": round(float(np.median(groups[k])), 3)}
                                         for k in keys}})
                NOTES.append((section, f"Subject-level test of % expressing across {GROUP_COL} groups: {stat_txt} "
                                       f"(exploratory; one value per subject avoids treating cells as replicates)."))
        else:
            plt.close(fig)

    if per_celltype and df[ct_col].nunique() > 1:
        g = df.groupby([donor_col, ct_col], observed=True)
        n = g.size().unstack()
        k = g["detected"].sum().unstack()
        pct = (100 * k / n).where(n >= min_group)
        mean = df[df.detected].groupby([donor_col, ct_col], observed=True)["expression"].mean().unstack()
        n.to_csv(out / "subject_x_celltype_n_cells.csv")
        k.to_csv(out / "subject_x_celltype_n_expressing.csv")
        save_table(pct.round(2), out / "subject_x_celltype_pct_expressing.csv", section,
                   f"% expressing per subject x cell type (blank = < {min_group} cells)")
        mean.round(4).to_csv(out / "subject_x_celltype_mean_expr_expressing.csv")
        col_order = pct.mean().sort_values(ascending=False).index
        row_order = pct.mean(axis=1).sort_values(ascending=False).index
        plot_heatmap(pct.loc[row_order, col_order], out / "heatmap_subject_x_celltype.png",
                     f"% {gene}+ cells: subject x cell type", section, "Subject x cell type heatmap")
        pos = ((k >= min_cells) & (100 * k / n >= min_pct) & (n >= min_group))
        cts = pd.DataFrame({
            "n_subjects_with_celltype": (n >= min_group).sum(),
            "n_subjects_expressing": pos.sum(),
        })
        cts["n_subjects_not_expressing"] = cts.n_subjects_with_celltype - cts.n_subjects_expressing
        cts["proportion_subjects_expressing"] = (cts.n_subjects_expressing
                                                 / cts.n_subjects_with_celltype.replace(0, np.nan)).round(4)
        cts["pct_subjects_expressing"] = (100 * cts["proportion_subjects_expressing"]).round(1)
        cts = cts.sort_values(["n_subjects_expressing", "pct_subjects_expressing"], ascending=False)
        cts.index.name = ct_col
        save_table(cts, out / "celltype_n_subjects_expressing.csv", section,
                   f"Per cell type: subjects expressing {gene} (subject needs >= {min_group} cells of that type)")
        plot_barh(cts["n_subjects_expressing"], out / "bar_n_subjects_expressing_by_celltype.png",
                  f"number of subjects expressing {gene}", f"Subjects with {gene}+ cells, per cell type", section,
                  "Subjects expressing the gene in each cell type (label: expressing/with cell type)",
                  colors="#8172b2", annotate=cts.apply(
                      lambda x: f"{int(x.n_subjects_expressing)}/{int(x.n_subjects_with_celltype)} "
                                f"({x.pct_subjects_expressing:.0f}%)", axis=1))
    return ds


# ----------------------------------------------------------------------------- DE: NEBULA (expression) + Fisher (proportions)
import shutil
import subprocess
import tempfile
import scipy.sparse as sp_sparse
from scipy.io import mmwrite

DE_SECTION = "4. Differential expression (gene+ vs gene- cells)"

NEBULA_R = r"""
args <- commandArgs(trailingOnly = TRUE); d <- args[1]; ncore <- as.integer(args[2])
suppressPackageStartupMessages({library(Matrix); library(nebula)})
cnt <- as(readMM(file.path(d, "counts.mtx")), "CsparseMatrix")          # genes x cells
rownames(cnt) <- readLines(file.path(d, "genes.txt"))
md <- read.csv(file.path(d, "meta.csv"), stringsAsFactors = FALSE)
colnames(cnt) <- md$cell
md$status <- factor(md$status, levels = c("neg", "pos"))
covs <- readLines(file.path(d, "covariates.txt")); covs <- covs[nchar(covs) > 0]
if ("grp" %in% covs) md$grp <- factor(md$grp)
f <- as.formula(paste(c("~ status", covs), collapse = " + "))
pred <- model.matrix(f, data = md)
o <- order(md$donor); cnt <- cnt[, o, drop = FALSE]; md <- md[o, ]; pred <- pred[o, , drop = FALSE]
fit <- nebula(cnt, md$donor, pred = pred, offset = md$libsize, method = "LN", ncore = ncore, verbose = FALSE)
s <- fit$summary
s$convergence <- fit$convergence
s$overdispersion_subject <- fit$overdispersion$Subject
s$overdispersion_cell <- fit$overdispersion$Cell
write.csv(s, file.path(d, "nebula_summary.csv"), row.names = FALSE)
writeLines(as.character(packageVersion("nebula")), file.path(d, "nebula_version.txt"))
"""

R_INSTALL = r"""
repos <- "https://cloud.r-project.org"
# nebula's declared minimum versions are too low for the API it uses; upgrade these if older
need <- function(p, v) !requireNamespace(p, quietly = TRUE) || packageVersion(p) < v
up <- c(parallelly = "1.38.0", future = "1.40.0", future.apply = "1.11.0", globals = "0.16.3", doFuture = "1.0.0")
todo <- names(up)[mapply(need, names(up), up)]
if (length(todo)) install.packages(todo, repos = repos)
if (!requireNamespace("nebula", quietly = TRUE)) {
  ok <- tryCatch({install.packages("nebula", repos = repos); TRUE}, error = function(e) FALSE)
  if (!requireNamespace("nebula", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)
    remotes::install_github("lhe17/nebula", upgrade = "never")
  }
}
ok <- tryCatch({
  suppressPackageStartupMessages(library(nebula)); data(sample_data, package = "nebula")
  re <- nebula(sample_data$count, sample_data$sid, pred = model.matrix(~X1, data = sample_data$pred),
               ncore = 1, verbose = FALSE)
  nrow(re$summary) > 0
}, error = function(e) {message("smoke test failed: ", conditionMessage(e)); FALSE})
cat("NEBULA_OK", ok, "\n")
"""

R_SMOKE = r"""
ok <- tryCatch({
  suppressPackageStartupMessages(library(nebula)); data(sample_data, package = "nebula")
  re <- nebula(sample_data$count, sample_data$sid, pred = model.matrix(~X1, data = sample_data$pred),
               ncore = 1, verbose = FALSE)
  "OK"
}, error = function(e) conditionMessage(e))
cat("SMOKE:", ok, "\n")
"""


def find_rscript(user_path=None):
    if user_path and Path(user_path).exists():
        return str(user_path)
    p = shutil.which("Rscript")
    if p:
        return p
    cands = list(Path("C:/Program Files/R").glob("R-*/bin/Rscript.exe")) + \
        [Path("/Library/Frameworks/R.framework/Resources/bin/Rscript"), Path("/opt/homebrew/bin/Rscript"),
         Path("/usr/local/bin/Rscript")]
    cands = [c for c in cands if c.exists()]
    return str(sorted(cands)[-1]) if cands else None


def r_status(rscript):
    """Return (ok, message)."""
    if not rscript:
        return False, "Rscript not found. Install R (https://cran.r-project.org) or enter the path to Rscript."
    try:
        r = subprocess.run([rscript, "-e", 'cat(as.character(getRversion()), requireNamespace("nebula", quietly=TRUE),'
                                           ' if (requireNamespace("nebula", quietly=TRUE)) as.character('
                                           'packageVersion("nebula")) else "")'],
                           capture_output=True, text=True, timeout=120)
        parts = r.stdout.split()
        if len(parts) >= 2 and parts[1] == "TRUE":
            sm = subprocess.run([rscript, "-e", R_SMOKE], capture_output=True, text=True, timeout=600)
            line = next((l for l in sm.stdout.splitlines() if l.startswith("SMOKE:")), "SMOKE: no output")
            if "SMOKE: OK" in line:
                return True, f"R {parts[0]} with nebula {parts[2] if len(parts) > 2 else ''} (test fit OK)"
            return False, (f"nebula {parts[2] if len(parts) > 2 else ''} is installed but a test fit failed "
                           f"({line[7:].strip()[:200]}). Press <b>Install nebula</b> to update its dependencies.")
        return False, f"R {parts[0] if parts else '?'} found, but the nebula package is not installed."
    except Exception as e:
        return False, f"Could not run Rscript: {e}"


def r_install_nebula(rscript):
    r = subprocess.run([rscript, "-e", R_INSTALL], capture_output=True, text=True, timeout=3600)
    ok = "NEBULA_OK TRUE" in r.stdout
    return ok, (r.stdout + r.stderr)[-3000:]


def bh(p):
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    m = ~np.isnan(p)
    if m.sum() == 0:
        return out
    pv = p[m]
    o = np.argsort(pv)
    ranked = pv[o] * m.sum() / np.arange(1, m.sum() + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty_like(pv)
    res[o] = np.minimum(q, 1)
    out[m] = res
    return out


def _fisher_row(a, b, c, d):
    o, p = stats.fisher_exact([[a, b], [c, d]])
    la = np.log((a + .5) * (d + .5) / ((b + .5) * (c + .5)))
    se = np.sqrt(1 / (a + .5) + 1 / (b + .5) + 1 / (c + .5) + 1 / (d + .5))
    return o, p, np.exp(la - 1.96 * se), np.exp(la + 1.96 * se)


def fisher_tables(df, ct_col, group_col):
    """(1) each cell type vs all other cells; (2) each pair of groups within each cell type."""
    scopes = [("all groups", df)]
    if group_col in df.columns and df[group_col].nunique() > 1:
        scopes += [(g, d) for g, d in df.groupby(group_col)]
    rows = []
    for scope, d in scopes:
        tot_pos, tot_n = int(d.detected.sum()), len(d)
        part = []
        for ct, x in d.groupby(ct_col):
            a = int(x.detected.sum()); b = len(x) - a
            c = tot_pos - a; e = (tot_n - len(x)) - c
            o, p, lo, hi = _fisher_row(a, b, c, e)
            part.append({"scope": scope, ct_col: ct, "pos_in_celltype": a, "neg_in_celltype": b,
                         "pos_in_other_cells": c, "neg_in_other_cells": e,
                         "pct_pos_celltype": round(100 * a / max(a + b, 1), 3),
                         "pct_pos_other_cells": round(100 * c / max(c + e, 1), 3),
                         "odds_ratio": o, "or_ci95_low": lo, "or_ci95_high": hi, "p_value": p})
        part = pd.DataFrame(part)
        if len(part):
            part["fdr"] = bh(part.p_value)
        rows.append(part)
    vs_rest = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    between = []
    if group_col in df.columns and df[group_col].nunique() > 1:
        gs = sorted(df[group_col].unique())
        for ct, x in df.groupby(ct_col):
            for i in range(len(gs)):
                for j in range(i + 1, len(gs)):
                    g1, g2 = gs[i], gs[j]
                    x1, x2 = x[x[group_col] == g1], x[x[group_col] == g2]
                    if len(x1) == 0 or len(x2) == 0:
                        continue
                    a, c = int(x1.detected.sum()), int(x2.detected.sum())
                    b, e = len(x1) - a, len(x2) - c
                    o, p, lo, hi = _fisher_row(a, b, c, e)
                    between.append({ct_col: ct, "group_1": g1, "group_2": g2, "pos_group_1": a, "neg_group_1": b,
                                    "pos_group_2": c, "neg_group_2": e,
                                    "pct_pos_group_1": round(100 * a / max(a + b, 1), 3),
                                    "pct_pos_group_2": round(100 * c / max(c + e, 1), 3),
                                    "odds_ratio_g1_vs_g2": o, "or_ci95_low": lo, "or_ci95_high": hi, "p_value": p})
    between = pd.DataFrame(between)
    if len(between):
        between["fdr"] = bh(between.p_value)
    return vs_rest, between


def read_count_rows(h5, path, rows, max_nnz=20_000_000, gap=200_000):
    """Rows (sorted, unique) of a CSR/CSC/dense matrix on disk -> scipy CSR (len(rows) x n_genes)."""
    node = h5[path]
    rows = np.asarray(rows, dtype=np.int64)
    if isinstance(node, h5py.Dataset):
        blocks = [node[rows[i:i + 2000]] for i in range(0, len(rows), 2000)]
        return sp_sparse.csr_matrix(np.vstack(blocks))
    enc = node.attrs.get("encoding-type", "csr_matrix")
    enc = enc.decode() if isinstance(enc, bytes) else enc
    shape = tuple(node.attrs["shape"])
    data, indices, indptr = node["data"], node["indices"], node["indptr"][:]
    if enc == "csc_matrix":
        pos = np.full(shape[0], -1, dtype=np.int64); pos[rows] = np.arange(len(rows))
        r_, c_, v_ = [], [], []
        g0 = 0
        while g0 < shape[1]:
            g1 = int(np.searchsorted(indptr, indptr[g0] + max_nnz, side="right")) - 1
            g1 = min(max(g1, g0 + 1), shape[1])
            a, b = indptr[g0], indptr[g1]
            ind, dat = indices[a:b], data[a:b]
            cols = np.repeat(np.arange(g0, g1), np.diff(indptr[g0:g1 + 1]))
            m = pos[ind] >= 0
            r_.append(pos[ind[m]]); c_.append(cols[m]); v_.append(dat[m])
            g0 = g1
        return sp_sparse.csr_matrix((np.concatenate(v_), (np.concatenate(r_), np.concatenate(c_))),
                                    shape=(len(rows), shape[1]))
    starts, ends = indptr[rows], indptr[rows + 1]
    lens = ends - starts
    out_data, out_ind = [], []
    i = 0
    while i < len(rows):
        j = i
        while j + 1 < len(rows) and starts[j + 1] - ends[j] <= gap and ends[j + 1] - starts[i] <= max_nnz:
            j += 1
        a, b = starts[i], ends[j]
        dat, ind = data[a:b], indices[a:b]
        sel = np.concatenate([np.arange(s - a, e - a) for s, e in zip(starts[i:j + 1], ends[i:j + 1])]) \
            if j > i else np.arange(0, b - a)
        out_data.append(dat[sel]); out_ind.append(ind[sel])
        i = j + 1
    ptr = np.concatenate([[0], np.cumsum(lens)])
    return sp_sparse.csr_matrix((np.concatenate(out_data) if out_data else np.array([]),
                                 np.concatenate(out_ind) if out_ind else np.array([], dtype=np.int64), ptr),
                                shape=(len(rows), shape[1]))


def counts_source(D):
    """Where integer counts live: (h5 path, gene ids, gene symbols) or None."""
    if D.has_raw:
        try:
            from anndata.io import read_elem
        except ImportError:
            from anndata.experimental import read_elem
        rv = read_elem(D.h5["raw"]["var"])
        ids = rv.index.astype(str).values
        sym = (rv["feature_name"].astype(str).values if "feature_name" in rv.columns else
               pd.Series(ids).map(D.symbols().to_dict()).fillna(pd.Series(ids)).values)
        return "raw/X", ids, sym
    if looks_like_counts(D.h5, "X"):
        return "X", D.var.index.astype(str).values, D.symbols().values
    return None


def plot_volcano(t, path, title, fdr_thr, section, label_n=12):
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    x = t["log2FC"].values
    y = -np.log10(np.clip(t["p_value"].values, 1e-300, 1))
    sig = (t["fdr"] < fdr_thr).values
    up, dn = sig & (x > 0), sig & (x < 0)
    ax.scatter(x[~sig], y[~sig], s=6, c="#bbbbbb", lw=0, rasterized=True, label="not significant")
    ax.scatter(x[up], y[up], s=9, c="#c0392b", lw=0, label=f"higher in gene+ ({up.sum()})")
    ax.scatter(x[dn], y[dn], s=9, c="#2471a3", lw=0, label=f"higher in gene- ({dn.sum()})")
    top = t[sig].nsmallest(label_n, "p_value")
    for _, r in top.iterrows():
        ax.text(r.log2FC, -np.log10(max(r.p_value, 1e-300)), str(r.gene)[:15], fontsize=6.5)
    ax.axvline(0, color="#888", lw=0.6)
    ax.set_xlabel("log2 fold change (gene+ vs gene- cells, NEBULA)")
    ax.set_ylabel("-log10 p")
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    for s_ in ["top", "right"]:
        ax.spines[s_].set_visible(False)
    savefig(fig, path, section, title)


def plot_between_heatmap(bt, ct_col, path, title, section):
    if bt.empty:
        return
    bt = bt.copy()
    bt["comparison"] = bt.group_1.astype(str) + " vs " + bt.group_2.astype(str)
    m = bt.pivot_table(index=ct_col, columns="comparison", values="odds_ratio_g1_vs_g2")
    q = bt.pivot_table(index=ct_col, columns="comparison", values="fdr")
    lm = np.log2(m.clip(lower=1e-3, upper=1e3))
    v = np.nanmax(np.abs(lm.values)) if np.isfinite(lm.values).any() else 1
    fig, ax = plt.subplots(figsize=(max(5, 1.3 * lm.shape[1] + 4), max(3, 0.3 * lm.shape[0] + 1.5)))
    cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#eeeeee")
    im = ax.imshow(np.ma.masked_invalid(lm.values), cmap=cmap, vmin=-v, vmax=v, aspect="auto")
    for i in range(lm.shape[0]):
        for j in range(lm.shape[1]):
            if q.values[i, j] == q.values[i, j] and q.values[i, j] < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=10)
    ax.set_xticks(range(lm.shape[1])); ax.set_xticklabels(lm.columns, rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(lm.shape[0])); ax.set_yticklabels([s[:45] for s in lm.index], fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.5).set_label("log2 odds ratio (group 1 vs group 2)")
    ax.set_title(title + "\n(* FDR < 0.05; red = more gene+ cells in group 1)", fontsize=10)
    savefig(fig, path, section, title)


def run_de(D, df, gene_symbol, gene_id, ct_col, donor_col, group_col, out, celltypes=None, scopes=("pooled", "within"),
           min_cells=20, min_donors=3, min_cells_per_donor=3, max_cells=5000, min_det=0.10, fdr_thr=0.05,
           rscript=None, ncore=1, progress=print, seed=0):
    """Fisher tests + NEBULA DE of gene+ vs gene- cells per cell type. df must contain '_row' (row in the H5AD)."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    sec = DE_SECTION
    for lst in (FIGS, TABLES, NOTES):
        lst[:] = [x for x in lst if x[0] != sec]
    result = {"folder": out, "jobs": [], "skipped": []}

    # ---- Fisher
    progress("Fisher's exact tests")
    vr, bt = fisher_tables(df, ct_col, group_col)
    save_table(vr, out / "fisher_celltype_vs_rest.csv", sec,
               "Fisher: gene+ fraction in each cell type vs all other cells", index=False)
    if len(bt):
        save_table(bt, out / "fisher_between_groups.csv", sec,
                   "Fisher: gene+ fraction between groups, within each cell type", index=False)
        plot_between_heatmap(bt, ct_col, out / "fisher_between_groups.png",
                             f"{gene_symbol}+ cells between groups, per cell type", sec)
    result["fisher_vs_rest"], result["fisher_between"] = vr, bt

    # ---- NEBULA
    src = counts_source(D)
    ok, msg = r_status(rscript)
    if src is None:
        NOTES.append((sec, "NEBULA skipped: this dataset has no integer counts (raw.X absent and X is normalized)."))
        result["nebula_msg"] = "no counts"
        return result
    if not ok:
        NOTES.append((sec, f"NEBULA skipped: {msg}"))
        result["nebula_msg"] = msg
        return result
    path, gids, gsyms = src
    rng = np.random.default_rng(seed)
    cts = celltypes or sorted(df[ct_col].unique())
    has_groups = group_col in df.columns and df[group_col].nunique() > 1
    summary = []
    for ct in cts:
        x = df[df[ct_col] == ct]
        jobs = []
        if "pooled" in scopes:
            jobs.append(("pooled", x))
        if "within" in scopes and has_groups:
            jobs += [(f"within: {g}", xg) for g, xg in x.groupby(group_col)]
        for scope, xs in jobs:
            npos, nneg = int(xs.detected.sum()), int((~xs.detected).sum())
            dd = xs.groupby(donor_col).detected.agg(["sum", "size"])
            nboth = int(((dd["sum"] >= min_cells_per_donor) & (dd["size"] - dd["sum"] >= min_cells_per_donor)).sum())
            base = {ct_col: ct, "scope": scope, "n_pos": npos, "n_neg": nneg, "donors_with_both": nboth}
            if npos < min_cells or nneg < min_cells or nboth < min_donors:
                summary.append({**base, "status": f"skipped (needs ≥{min_cells} gene+ and gene- cells and "
                                                  f"≥{min_donors} donors with ≥{min_cells_per_donor} of each)"})
                continue
            pos_rows, neg_rows = xs.loc[xs.detected, "_row"].values, xs.loc[~xs.detected, "_row"].values
            if len(pos_rows) > max_cells:
                pos_rows = rng.choice(pos_rows, max_cells, replace=False)
            if len(neg_rows) > max_cells:
                neg_rows = rng.choice(neg_rows, max_cells, replace=False)
            rows = np.sort(np.concatenate([pos_rows, neg_rows]))
            progress(f"NEBULA: {ct} [{scope}] - reading counts for {len(rows):,} cells")
            M = read_count_rows(D.h5, path, rows).tocsr()
            M.data = np.rint(M.data)
            lib = np.asarray(M.sum(1)).ravel()
            ngenes = np.diff(M.indptr)
            meta = df.set_index("_row").loc[rows]
            status = np.where(np.isin(rows, pos_rows), "pos", "neg")
            det_pos = np.asarray((M[status == "pos"] > 0).mean(0)).ravel()
            det_neg = np.asarray((M[status == "neg"] > 0).mean(0)).ravel()
            keep = ((det_pos >= min_det) | (det_neg >= min_det)) & (gids != gene_id) & (gsyms != gene_symbol)
            if keep.sum() == 0:
                summary.append({**base, "status": "skipped (no genes pass the detection filter)"})
                continue
            covs = ["z_log_ngenes"]
            if scope == "pooled" and has_groups and meta[group_col].nunique() > 1:
                covs.append("grp")
            lng = np.log(np.maximum(ngenes, 1))
            md = pd.DataFrame({"cell": [f"c{r}" for r in rows], "donor": meta[donor_col].astype(str).values,
                               "status": status, "libsize": np.maximum(lib, 1),
                               "z_log_ngenes": (lng - lng.mean()) / (lng.std() or 1),
                               "grp": meta[group_col].astype(str).values if group_col in meta.columns else "all"})
            work = Path(tempfile.mkdtemp(prefix="nebula_"))
            mmwrite(str(work / "counts.mtx"), M[:, keep].T.astype(np.int64).tocoo(), field="integer")
            (work / "genes.txt").write_text("\n".join(gids[keep]) + "\n")
            md.to_csv(work / "meta.csv", index=False)
            (work / "covariates.txt").write_text("\n".join(covs) + "\n")
            (work / "run_nebula.R").write_text(NEBULA_R)
            progress(f"NEBULA: {ct} [{scope}] - fitting {keep.sum():,} genes on {len(rows):,} cells "
                     f"({md.donor.nunique()} donors)")
            r = subprocess.run([rscript, str(work / "run_nebula.R"), str(work), str(int(ncore))],
                               capture_output=True, text=True)
            if r.returncode != 0 or not (work / "nebula_summary.csv").exists():
                err = (r.stderr or r.stdout)
                err = err[err.find("Error"):] if "Error" in err else err[-400:]
                summary.append({**base, "status": "NEBULA error: " + err[:400].replace("\n", " ")})
                continue
            s = pd.read_csv(work / "nebula_summary.csv")
            nv = (work / "nebula_version.txt").read_text().strip() if (work / "nebula_version.txt").exists() else ""
            shutil.rmtree(work, ignore_errors=True)
            idmap = dict(zip(gids, gsyms))
            norm = M.multiply(1e4 / np.maximum(lib, 1)[:, None]).tocsr()
            norm.data = np.log1p(norm.data)
            mpos = np.asarray(norm[status == "pos"].mean(0)).ravel()
            mneg = np.asarray(norm[status == "neg"].mean(0)).ravel()
            gidx = {g: i for i, g in enumerate(gids)}
            if "gene" not in s.columns:      # nebula omits names in rare edge cases: map by position
                s["gene"] = gids[keep][s["gene_id"].astype(int).values - 1]
            t = pd.DataFrame({"gene_id": s["gene"].astype(str)})
            t["gene"] = t.gene_id.map(idmap)
            t["log2FC"] = s["logFC_statuspos"] / np.log(2)
            t["se_log2FC"] = s["se_statuspos"] / np.log(2)
            t["p_value"] = s["p_statuspos"]
            t["convergence"] = s["convergence"]
            t["reliable_fit"] = t.convergence > -20
            t["fdr"] = np.nan
            t.loc[t.reliable_fit, "fdr"] = bh(t.loc[t.reliable_fit, "p_value"])
            ii = t.gene_id.map(gidx).values
            t["pct_detected_pos"] = (100 * det_pos[ii]).round(2)
            t["pct_detected_neg"] = (100 * det_neg[ii]).round(2)
            t["mean_lognorm_pos"] = mpos[ii].round(4)
            t["mean_lognorm_neg"] = mneg[ii].round(4)
            t["overdispersion_subject"] = s.get("overdispersion_subject")
            t = t.sort_values(["fdr", "p_value"])
            jd = out / safe_name(ct) / safe_name(scope.replace("within: ", "within_"))
            jd.mkdir(parents=True, exist_ok=True)
            t.to_csv(jd / "nebula_de_results.csv", index=False)
            plot_volcano(t.dropna(subset=["fdr"]), jd / "volcano.png", f"{ct} [{scope}]: {gene_symbol}+ vs "
                         f"{gene_symbol}- cells", fdr_thr, sec)
            nup = int(((t.fdr < fdr_thr) & (t.log2FC > 0)).sum()); ndn = int(((t.fdr < fdr_thr) & (t.log2FC < 0)).sum())
            summary.append({**base, "status": "done", "n_cells_used": len(rows), "n_donors": md.donor.nunique(),
                            "n_genes_tested": int(keep.sum()), "n_unreliable_fits": int((~t.reliable_fit).sum()),
                            f"n_up_fdr<{fdr_thr:g}": nup, f"n_down_fdr<{fdr_thr:g}": ndn,
                            "covariates": " + ".join(["status"] + covs), "nebula_version": nv,
                            "top_up": ", ".join(t[(t.fdr < fdr_thr) & (t.log2FC > 0)].gene.astype(str).head(8)),
                            "top_down": ", ".join(t[(t.fdr < fdr_thr) & (t.log2FC < 0)].gene.astype(str).head(8))})
            result["jobs"].append({"celltype": ct, "scope": scope, "table": t, "folder": jd})
    sm = pd.DataFrame(summary)
    save_table(sm, out / "nebula_summary.csv", sec, "NEBULA runs per cell type and scope", index=False)
    result["summary"] = sm
    if not result["jobs"]:
        errs = sm[sm.status.str.startswith("NEBULA error")] if len(sm) else sm
        result["nebula_msg"] = (f"NEBULA failed for every cell type: {errs.status.iloc[0][14:]}" if len(errs) else
                                "every selected cell type was below the thresholds (see the NEBULA runs table)")

    pooled = [j for j in result["jobs"] if j["scope"] == "pooled"] or result["jobs"]
    if pooled:
        top = set()
        for j in pooled:
            t = j["table"]
            sig = t[t.fdr < fdr_thr]
            top |= set(sig[sig.log2FC > 0].gene.astype(str).head(5)) | set(sig[sig.log2FC < 0].gene.astype(str).head(5))
        top = sorted(top)[:60]
        if top:
            lm = pd.DataFrame({f"{j['celltype']}" + ("" if j["scope"] == "pooled" else f" [{j['scope']}]"):
                               j["table"].drop_duplicates("gene").set_index("gene").reindex(top)["log2FC"]
                               for j in pooled}).T
            qm = pd.DataFrame({f"{j['celltype']}" + ("" if j["scope"] == "pooled" else f" [{j['scope']}]"):
                               j["table"].drop_duplicates("gene").set_index("gene").reindex(top)["fdr"]
                               for j in pooled}).T
            lm.to_csv(out / "top_genes_log2FC_matrix.csv")
            v = np.nanpercentile(np.abs(lm.values), 98) if np.isfinite(lm.values).any() else 1
            fig, ax = plt.subplots(figsize=(min(30, max(6, 0.32 * lm.shape[1] + 3)), max(2.5, 0.35 * lm.shape[0] + 2)))
            cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#eeeeee")
            im = ax.imshow(np.ma.masked_invalid(lm.values.astype(float)), cmap=cmap, vmin=-v, vmax=v, aspect="auto")
            for i in range(lm.shape[0]):
                for k in range(lm.shape[1]):
                    if qm.values[i, k] == qm.values[i, k] and qm.values[i, k] < fdr_thr:
                        ax.text(k, i, "*", ha="center", va="center", fontsize=8)
            ax.set_xticks(range(lm.shape[1])); ax.set_xticklabels(lm.columns, rotation=75, ha="right", fontsize=7)
            ax.set_yticks(range(lm.shape[0])); ax.set_yticklabels([s[:40] for s in lm.index], fontsize=7)
            fig.colorbar(im, ax=ax, shrink=0.5).set_label("log2FC gene+ vs gene- (NEBULA)")
            ax.set_title(f"Top DE genes per cell type ({gene_symbol}+ vs {gene_symbol}- cells)\n"
                         f"(* FDR < {fdr_thr:g}; grey = not tested)", fontsize=10)
            savefig(fig, out / "heatmap_top_de_genes.png", sec, "Top DE genes across cell types")
    done = sm[sm.status == "done"] if len(sm) else sm
    NOTES.append((sec, f"NEBULA ran for {len(done)} cell type × scope combination(s); "
                       f"{len(sm) - len(done)} were skipped or failed (see table). Fisher tests ran for all cell types."))
    return result


# ----------------------------------------------------------------------------- column choosers
def choose_column(obs, preferred, hints, what, preset=None):
    if preset:
        if preset not in obs.columns:
            sys.exit(f"Column '{preset}' not in obs. Available: {list(obs.columns)}")
        return preset
    if preferred in obs.columns:          # standard CELLxGENE annotation: use it, no question
        return preferred
    say(f"   Standard CELLxGENE column '{preferred}' not found.")
    cands = []
    cands += [c for c in obs.columns if any(h in c.lower() for h in hints) and c not in cands
              and not c.endswith("ontology_term_id") and 1 < obs[c].nunique() < max(2, len(obs) // 2)]
    if not cands:
        return None
    say(f"\nColumns usable as {what}:")
    for i, c in enumerate(cands, 1):
        say(f"  [{i}] {c}  ({obs[c].nunique()} levels)")
    ans = ask(f"Which column for {what}? (number or name)", "1")
    if ans.isdigit():
        return cands[int(ans) - 1]
    return ans if ans in obs.columns else cands[0]


# ----------------------------------------------------------------------------- report
def img_b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


FIRST_SECTION = {}


def write_report(out, meta, interpretations=None):
    interpretations = interpretations or {}
    sections = list(dict.fromkeys([f[0] for f in FIGS] + [t[0] for t in TABLES] + [n[0] for n in NOTES]))
    FIRST_SECTION.clear()
    for s in sections:
        k = {"1.": "1_cellular", "2.": "2_cell_specific", "3.": "3_sample"}.get(s[:2])
        FIRST_SECTION.setdefault(k, s)
    html = [f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{meta.get('Gene','Gene')} - CELLxGENE report</title><style>
body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1150px;margin:auto;padding:24px;color:#222}}
h1{{font-size:24px}} h2{{border-bottom:2px solid #ddd;padding-bottom:4px;margin-top:40px}}
.meta td{{padding:2px 10px;vertical-align:top}} .grid{{display:flex;flex-wrap:wrap;gap:18px}}
figure{{margin:0;flex:1 1 480px}} figure img{{max-width:100%;border:1px solid #eee}}
figcaption{{font-size:13px;color:#555}} table.df{{border-collapse:collapse;font-size:12px}}
table.df td,table.df th{{border:1px solid #ddd;padding:3px 6px}} .tbl{{overflow-x:auto;margin:10px 0}}
.note{{background:#f6f8fa;padding:8px 12px;border-left:4px solid #4c72b0;margin:8px 0}}
code{{background:#f2f2f2;padding:1px 4px}}
.interp{{background:#fbf8ef;border:1px solid #e6dcc0;border-left:5px solid #c9a227;padding:10px 14px;margin:10px 0}}
.methods p{{line-height:1.5}}</style></head><body>
<h1>{meta.get('Gene','Gene')} expression report</h1><table class="meta">"""]
    for k, v in meta.items():
        html.append(f"<tr><td><b>{k}</b></td><td>{v}</td></tr>")
    html.append("</table>")
    html.append('<h2>Methods</h2><div class="methods">' + md_to_html(methods_markdown()) + "</div>")
    for s in sections:
        html.append(f"<h2>{s}</h2>")
        key = {"1.": "1_cellular", "2.": "2_cell_specific", "3.": "3_sample"}.get(s[:2])
        if key and interpretations.get(key) and s in (FIRST_SECTION.get(key),):
            html.append(f'<div class="interp"><b>Interpretation ({LLM_MODEL}, AI-generated; verify against the data)</b>'
                        + md_to_html(interpretations[key]) + "</div>")
        for sec, txt in NOTES:
            if sec == s:
                html.append(f'<div class="note">{txt}</div>')
        for sec, cap, path, prev in TABLES:
            if sec == s:
                html.append(f"<h4>{cap} <small>(<code>{path.relative_to(out)}</code>)</small></h4>"
                            f'<div class="tbl">{prev.to_html(classes="df", border=0)}</div>')
        html.append('<div class="grid">')
        for sec, cap, path in FIGS:
            if sec == s:
                html.append(f'<figure><img src="data:image/png;base64,{img_b64(path)}">'
                            f"<figcaption>{cap} - <code>{path.relative_to(out)}</code></figcaption></figure>")
        html.append("</div>")
    html.append("</body></html>")
    p = out / "report.html"
    p.write_text("\n".join(html))
    return p


