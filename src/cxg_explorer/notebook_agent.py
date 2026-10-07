# cellxgeneExplorer
# Copyright (c) 2026 TeraWeave Technologies. Licensed under the PolyForm Noncommercial License 1.0.0 (see LICENSE.md).
# Required Notice: Copyright (c) 2026 TeraWeave Technologies
"""Interactive Jupyter panel for the single-study gene explorer. Use: `from cxg_explorer import explorer; explorer()`."""
from .agent import *  # noqa: F401,F403  (engine functions and shared state)
from .agent import _f  # noqa: F401
from . import agent as A
from . import interactive as IX

import ipywidgets as W
from IPython.display import display, HTML, clear_output

DATA_DIR = Path("cellxgene_data"); DATA_DIR.mkdir(exist_ok=True)
OUT_DIR = Path("cellxgene_results")
STATE = {}
L = lambda w="95%": W.Layout(width=w)

title_in = W.Textarea(placeholder="Paste the paper title as listed on CELLxGENE (or a DOI / CELLxGENE collection link)",
                      layout=W.Layout(width="95%", height="64px"))
search_btn = W.Button(description="Search CELLxGENE", button_style="primary", layout=L("200px"))
coll_dd = W.Dropdown(description="Collection", layout=L(), disabled=True, style={"description_width": "90px"})
ds_btn = W.Button(description="Show datasets", disabled=True, layout=L("200px"))
ds_dd = W.Dropdown(description="Dataset", layout=L(), disabled=True, style={"description_width": "90px"})
dl_btn = W.Button(description="Download & open", button_style="info", disabled=True, layout=L("200px"))
local_in = W.Text(placeholder="...or path to an H5AD already downloaded from CELLxGENE", layout=L("70%"))
local_btn = W.Button(description="Open local file", layout=L("200px"))
study_in = W.Text(description="Study folder", placeholder="name of the results folder for this study",
                  layout=L("70%"), style={"description_width": "90px"})
study_in.observe(lambda ch: STATE.update(study_custom=True) if not STATE.get("_setting_study") else None, "value")


def default_study_name(coll=None, fallback=""):
    """First author + year (e.g. 'Smith2021'), else a shortened collection name / file name."""
    if coll:
        pm = coll.get("publisher_metadata") or {}
        authors = pm.get("authors") or []
        fam = ""
        if authors:
            a0 = authors[0]
            fam = a0.get("family") or a0.get("name", "") if isinstance(a0, dict) else str(a0).split()[-1]
        if fam and pm.get("published_year"):
            return safe_name(f"{fam}{pm['published_year']}", 40)
        return safe_name(coll.get("name", ""), 40)
    return safe_name(fallback, 40)


def set_study(name):
    if STATE.get("study_custom") and study_in.value.strip():
        return                       # keep what the user typed
    STATE["_setting_study"] = True
    study_in.value = name
    STATE["_setting_study"] = False


group_col_dd = W.Dropdown(description="Group by", disabled=True, layout=L("420px"), style={"description_width": "90px"})
group_sel = W.SelectMultiple(description="Groups", rows=6, layout=L(), disabled=True, style={"description_width": "90px"})
grp_all_btn = W.Button(description="Select all groups", disabled=True, layout=L("180px"))
group_info = W.HTML()


def group_columns(obs):
    """Columns usable as condition/group: 'disease' first, then other low-cardinality annotations."""
    pref = ["disease", "condition", "disease_state", "group", "treatment", "stage", "severity", "tissue",
            "development_stage", "sex", "tissue_type"]
    ok = lambda c: (not c.endswith("ontology_term_id") and 1 < obs[c].nunique() <= 30
                    and c not in ("cell_type", "donor_id", "assay", "suspension_type", "is_primary_data"))
    cols = [c for c in pref if c in obs.columns and ok(c)]
    cols += [c for c in obs.columns if c not in cols and ok(c)
             and any(h in c.lower() for h in ["disease", "condition", "status", "group", "severity", "treat", "outcome"])]
    return cols


def refresh_groups(_=None):
    D = STATE.get("D")
    col = group_col_dd.value
    if D is None or not col:
        group_sel.options, group_info.value = [], ""
        return
    vc = D.obs[col].astype(str).value_counts()
    nd = D.obs.groupby(col, observed=True)["donor_id"].nunique() if "donor_id" in D.obs.columns else None
    group_sel.options = [(f"{k}  ({v:,} cells" + (f", {int(nd.get(k, 0))} donors)" if nd is not None else ")"), k)
                         for k, v in vc.items()]
    group_sel.value = tuple(vc.index)
    update_group_info()


def update_group_info(_=None):
    n = len(group_sel.value)
    tot = len(group_sel.options)
    msg = ("all groups" if n == tot else f"{n} of {tot} groups") + " selected"
    if n == 2:
        msg += ": direct two-group comparison (Mann-Whitney U on subjects)"
    elif n > 2:
        msg += ": multi-group comparison (Kruskal-Wallis on subjects)"
    elif n == 1:
        msg += ": single group, no between-group comparison"
    else:
        msg = "no group selected"
    group_info.value = f"<small>{msg}</small>"


group_col_dd.observe(refresh_groups, "value")
group_sel.observe(update_group_info, "value")
grp_all_btn.on_click(lambda _: setattr(group_sel, "value", tuple(v for _, v in group_sel.options)))


gene_in = W.Text(description="Gene", placeholder="e.g. CD274", disabled=True, style={"description_width": "90px"})
ct_sel = W.SelectMultiple(description="Populations", rows=10, layout=L(), disabled=True,
                          style={"description_width": "90px"})
min_cells_in = W.BoundedIntText(value=3, min=1, max=10000, description="min cells", layout=L("190px"),
                                tooltip="a patient expresses the gene if it has >= this many expressing cells ...")
min_pct_in = W.BoundedFloatText(value=1.0, min=0, max=100, step=0.5, description="min %", layout=L("190px"),
                                tooltip="... and >= this % of its cells express it")
llm_model_in = W.Text(value=LLM_MODEL, description="Model", layout=L("260px"))
llm_url_in = W.Text(value=LLM_URL, description="Server", layout=L("340px"))
report_btn = W.Button(description="Generate descriptive report", button_style="warning",
                      disabled=True, layout=L("260px"), tooltip="Writes a narrative report with gpt-oss:20b")
rscript_in = W.Text(value=find_rscript() or "", description="Rscript", layout=L("520px"),
                    placeholder="path to Rscript (auto-detected)", style={"description_width": "90px"})
r_check_btn = W.Button(description="Check R / nebula", layout=L("170px"))
r_install_btn = W.Button(description="Install nebula", layout=L("170px"))
r_status_html = W.HTML()
de_ct_sel = W.SelectMultiple(description="Cell types", rows=7, layout=L(), disabled=True,
                             style={"description_width": "90px"})
de_pooled = W.Checkbox(value=True, description="Pooled groups (group as covariate)", indent=False, layout=L("300px"))
de_within = W.Checkbox(value=True, description="Within each group separately", indent=False, layout=L("300px"))
de_min = W.BoundedIntText(value=20, min=3, max=10000, description="min cells/class", layout=L("200px"),
                          style={"description_width": "110px"})
de_max = W.BoundedIntText(value=5000, min=50, max=1000000, description="max cells/class", layout=L("220px"),
                          style={"description_width": "110px"})
de_det = W.BoundedFloatText(value=10, min=0, max=100, step=1, description="min detect %", layout=L("200px"),
                            style={"description_width": "110px"})
de_fdr = W.BoundedFloatText(value=0.05, min=0.0001, max=1, step=0.01, description="FDR", layout=L("160px"),
                            style={"description_width": "50px"})
de_cores = W.BoundedIntText(value=1, min=1, max=128, description="R cores", layout=L("160px"),
                            style={"description_width": "60px"})
de_btn = W.Button(description="Run DE", button_style="danger", disabled=True, layout=L("200px"))
run_btn = W.Button(description="Run analysis", button_style="success", disabled=True, layout=L("200px"))
status = W.HTML()
log = W.Output(layout=W.Layout(max_height="260px", overflow="auto", border="1px solid #ddd"))
results = W.Output()


def set_status(msg, kind="info"):
    color = {"info": "#1f5f9e", "ok": "#2e7d32", "warn": "#b26a00", "err": "#c62828"}[kind]
    status.value = f'<div style="padding:6px 10px;border-left:4px solid {color};background:#f7f9fb">{msg}</div>'


def guarded(fn):
    def wrapper(_=None):
        try:
            with log:
                fn()
        except Exception as e:
            set_status(f"{type(e).__name__}: {e}", "err")
            if "403" in str(e):
                set_status("CELLxGENE was blocked by the network (403). Check internet/proxy access to "
                           "api.cellxgene.cziscience.com and datasets.cellxgene.cziscience.com.", "err")
    return wrapper


@guarded
def on_search():
    q = title_in.value.strip()
    if not q:
        return set_status("Paste a paper title first.", "warn")
    set_status("Searching CELLxGENE Discover ...")
    hits = find_collections(q, DATA_DIR)
    STATE["hits"] = hits
    coll_dd.options = [(f"[{min(s, 1):.2f} · {why}] {c.get('name')}  ({len(c.get('datasets') or [])} datasets)", i)
                       for i, (s, c, why) in enumerate(hits)]
    coll_dd.disabled, ds_btn.disabled = False, False
    best = hits[0][0] if hits else 0
    set_status("Top match selected; check it is your paper." if best >= 0.45 else
               "No strong match. Try the DOI or the CELLxGENE collection link.", "ok" if best >= 0.45 else "warn")


@guarded
def on_datasets():
    c = STATE["hits"][coll_dd.value][1]
    if not c.get("datasets") or "assets" not in (c["datasets"][0] or {}):
        c = get_json(f"{CURATION_API}/collections/{c['collection_id']}")
    STATE["coll"] = c
    set_study(default_study_name(c))
    ds = sorted([d for d in c.get("datasets", []) if d], key=lambda d: -(d.get("cell_count") or 0))
    STATE["datasets"] = ds
    def lab(d):
        h5 = next((a for a in d.get("assets", []) if a.get("filetype") == "H5AD"), {})
        return (f"{d.get('title')} | {d.get('cell_count', 0):,} cells | {labels_of(d.get('assay'))} | "
                f"H5AD {human_size(h5.get('filesize'))}")
    ds_dd.options = [(lab(d), i) for i, d in enumerate(ds)]
    ds_dd.disabled, dl_btn.disabled = False, False
    set_status(f"{len(ds)} dataset(s) in <b>{c.get('name')}</b>. "
               f'<a href="{COLLECTION_URL}{c["collection_id"]}" target="_blank">open on CELLxGENE</a>')


def open_dataset(path, meta):
    set_status("Opening dataset ...")
    D = CxgDataset(path)
    side = Path(str(path) + ".cellxgene.json")
    cx = json.loads(side.read_text()) if side.exists() else {}
    STATE.update(D=D, meta=meta, cx=cx)
    ct_col = "cell_type" if "cell_type" in D.obs.columns else D.obs.columns[0]
    vc = D.obs[ct_col].astype(str).value_counts()
    ct_sel.options = [(f"{k}  ({v:,} cells)", k) for k, v in vc.items()]
    gcols = group_columns(D.obs)
    group_col_dd.options = gcols or []
    group_col_dd.disabled = group_sel.disabled = grp_all_btn.disabled = not gcols
    if gcols:
        group_col_dd.value = gcols[0]
        refresh_groups()
    else:
        group_info.value = "<small>No group/condition column with 2-30 levels found; all cells are analysed.</small>"
    gene_in.disabled = ct_sel.disabled = run_btn.disabled = False
    donors = D.obs["donor_id"].nunique() if "donor_id" in D.obs.columns else "?"
    set_status(f"Loaded <b>{D.n_obs:,} cells</b>, {vc.size} cell types, {donors} donors. "
               "Now type a gene, select population(s) and press <b>Run analysis</b>.", "ok")


@guarded
def on_download():
    d = STATE["datasets"][ds_dd.value]
    c = STATE["coll"]
    set_status("Downloading from CELLxGENE (progress in the log) ...")
    path = download_h5ad(d, DATA_DIR)
    try:   # keep CELLxGENE collection/dataset metadata next to the file (citation, description, links)
        Path(str(path) + ".cellxgene.json").write_text(json.dumps({"collection": c, "dataset": d}, default=str))
    except Exception:
        pass
    open_dataset(path, {"Collection": c.get("name"), "DOI": c.get("doi") or "-",
                        "Collection URL": f'<a href="{COLLECTION_URL}{c["collection_id"]}">{COLLECTION_URL}{c["collection_id"]}</a>',
                        "Dataset": d.get("title"), "Dataset ID": d.get("dataset_id")})


@guarded
def on_local():
    p = Path(local_in.value.strip())
    if not p.exists():
        return set_status(f"File not found: {p}", "warn")
    set_study(default_study_name(fallback=p.stem.split("__")[0]))
    open_dataset(p, {"Source": f"local file {p.name}"})


@guarded
def on_run():
    D = STATE["D"]
    if not group_col_dd.disabled and not group_sel.value:
        return set_status("Select at least one group.", "warn")
    var_id, symbol, sugg = D.resolve_gene(gene_in.value)
    if var_id is None:
        return set_status(f"'{gene_in.value}' not found." + (f" Did you mean: {', '.join(sugg)}?" if sugg else
                          " Use the official symbol (e.g. CD274 for PD-L1)."), "warn")
    set_status(f"Analysing {symbol} ...")
    res = run_analysis(D, var_id, symbol, list(ct_sel.value), min_cells_in.value, min_pct_in.value, STATE["meta"],
                       study=study_in.value.strip() or None,
                       group_col=group_col_dd.value if not group_col_dd.disabled else None,
                       groups=list(group_sel.value) if not group_col_dd.disabled else None)
    set_status(f"Analysis done. Files saved in <code>{res['folder']}</code>", "ok")
    show_results(res)


SAMPLE_HINTS = re.compile(r"sample|specimen|biopsy|procedure|collection|sampling|source|site|location|method|"
                          r"protocol|resection|surgery|fixation|preserv|dissociat|sort|enrich|time_?point|visit|"
                          r"treatment|region|anatomical|compartment|fraction|isolation|digest", re.I)
SCHEMA_COLS = {"cell_type", "donor_id", "disease", "tissue", "sex", "assay", "organism", "development_stage",
               "self_reported_ethnicity", "suspension_type", "tissue_type", "is_primary_data", "observation_joinid"}


def author_sample_columns(obs, n_donors):
    """Author-provided columns that describe samples or how they were obtained."""
    desc, sample_col = {}, None
    for c in obs.columns:
        if c in SCHEMA_COLS or c.endswith("ontology_term_id") or not SAMPLE_HINTS.search(c):
            continue
        nun = obs[c].nunique()
        if nun <= 1:
            continue
        if sample_col is None and re.search(r"sample|specimen|library", c, re.I) and n_donors <= nun <= 5000:
            sample_col = c
        if nun <= 50 and c != sample_col:
            desc[c] = {str(k): int(v) for k, v in obs[c].astype(str).value_counts().head(10).items()}
    if sample_col:
        desc.pop(sample_col, None)
    return desc, sample_col


def parse_uns_citation(text):
    """CELLxGENE H5ADs carry uns['citation'] with the publication DOI and collection link."""
    t = str(text or "")
    doi = re.search(r"doi\.org/(10\.[^\s]+)", t)
    cid = re.search(r"collections/([0-9a-f-]{36})", t)
    return (doi.group(1).rstrip(".") if doi else None), (cid.group(1) if cid else None)


def write_study_info(out, D, symbol, var_id, study, meta, keep, all_groups, pops):
    """Machine-readable description of this run (+ samples_used.csv), used by the cross-study summary notebook."""
    obs = D.obs if keep is None else D.obs[keep]
    vc = lambda c: obs[c].astype(str).value_counts().to_dict() if c in obs.columns else {}
    onto = {}
    if "cell_type_ontology_term_id" in obs.columns and RUN_INFO.get("ct_col") in obs.columns:
        pairs = obs[[RUN_INFO["ct_col"], "cell_type_ontology_term_id"]].astype(str).drop_duplicates()
        onto = dict(zip(pairs.iloc[:, 0], pairs.iloc[:, 1]))
    prim = (float(obs["is_primary_data"].astype(str).str.lower().isin(["true", "1"]).mean())
            if "is_primary_data" in obs.columns else None)
    cx = STATE.get("cx") or {}
    coll, dset = cx.get("collection") or {}, cx.get("dataset") or {}
    uns_cit = str(D.adata.uns.get("citation", "")) if "citation" in D.adata.uns else ""
    u_doi, u_cid = parse_uns_citation(uns_cit)
    cid = re.search(r"collections/([0-9a-f-]{36})", str(meta.get("Collection URL", "")))
    pm = coll.get("publisher_metadata") or {}
    donor_col = RUN_INFO.get("donor_col") if RUN_INFO.get("donor_col") in obs.columns else None
    n_donors = obs[donor_col].nunique() if donor_col else 0
    author_cols, sample_col = author_sample_columns(obs, n_donors)

    info = {"gene": symbol, "gene_id": var_id, "study": study, "dataset": RUN_INFO.get("dataset"),
            "dataset_id": meta.get("Dataset ID") or dset.get("dataset_id"),
            "collection": RUN_INFO.get("collection") if RUN_INFO.get("collection") != "n/a" else coll.get("name"),
            "collection_id": (cid.group(1) if cid else None) or coll.get("collection_id") or u_cid,
            "doi": (RUN_INFO.get("doi") if RUN_INFO.get("doi") not in (None, "n/a", "-") else None)
                   or coll.get("doi") or u_doi,
            "citation": {"authors": [{"family": a.get("family"), "given": a.get("given"), "name": a.get("name")}
                                     for a in (pm.get("authors") or []) if isinstance(a, dict)],
                         "journal": pm.get("journal"), "year": pm.get("published_year"),
                         "month": pm.get("published_month"), "cellxgene_citation": uns_cit},
            "collection_description": coll.get("description"), "consortia": coll.get("consortia"),
            "contact": coll.get("contact_name"),
            "links": [{"type": l.get("link_type"), "name": l.get("link_name"), "url": l.get("link_url")}
                      for l in (coll.get("links") or [])],
            "n_cells": int(len(obs)), "n_subjects": RUN_INFO.get("n_subjects"),
            "n_samples": int(obs[sample_col].nunique()) if sample_col else None, "sample_column": sample_col,
            "celltype_col": RUN_INFO.get("ct_col"), "donor_col": RUN_INFO.get("donor_col"),
            "group_col": A.GROUP_COL, "groups": RUN_INFO.get("groups"), "all_groups": all_groups,
            "groups_subset": bool(RUN_INFO.get("groups_subset")), "populations": pops,
            "tissue": vc("tissue"), "tissue_type": vc("tissue_type"), "assay": vc("assay"),
            "suspension_type": vc("suspension_type"), "organism": vc("organism"), "sex": vc("sex"),
            "development_stage": vc("development_stage"), "self_reported_ethnicity": vc("self_reported_ethnicity"),
            "author_sample_annotations": author_cols,
            "is_primary_fraction": prim, "expression_values": PARAMS["x_def"], "has_raw_counts": bool(D.has_raw),
            "min_cells": PARAMS["min_cells"], "min_pct": PARAMS["min_pct"], "celltype_ontology": onto,
            "created": time.strftime("%Y-%m-%d %H:%M")}
    (Path(out) / "study_info.json").write_text(json.dumps(info, indent=2, default=str))

    # one row per donor (x sample, if the authors provide a sample column) actually used in this run
    keys = [c for c in [donor_col, sample_col] if c]
    if keys:
        cols = [c for c in dict.fromkeys([A.GROUP_COL, "tissue", "tissue_type", "suspension_type", "assay", "sex",
                                          "development_stage", "self_reported_ethnicity"] + list(author_cols))
                if c in obs.columns and c not in keys]
        o = obs[keys + cols].astype(str)
        su = o.groupby(keys, observed=True).agg({c: (lambda s: "; ".join(sorted(set(s)))) for c in cols})
        su.insert(0, "n_cells", o.groupby(keys, observed=True).size())
        su.reset_index().to_csv(Path(out) / "samples_used.csv", index=False)


def run_analysis(D, var_id, symbol, pops, min_cells=3, min_pct=1.0, meta=None, min_group=10, study=None,
                 group_col=None, groups=None):
    A.GROUP_COL = group_col or "disease"
    FIGS.clear(); TABLES.clear(); NOTES.clear(); STATS.clear()
    PARAMS.update(min_cells=min_cells, min_pct=min_pct, min_group=min_group)
    expr, detected, value_desc = D.gene_expression(var_id)
    ct_col = "cell_type" if "cell_type" in D.obs.columns else choose_column(D.obs, "cell_type", ["cell"], "cell types")
    donor_col = "donor_id" if "donor_id" in D.obs.columns else None
    df = pd.DataFrame({"cell_id": D.obs.index.astype(str)})
    for c in set([ct_col, donor_col, A.GROUP_COL] + covariates_present(D.obs) + ["donor_id", "disease", "tissue", "sex", "development_stage"]):
        if c and c in D.obs.columns:
            df[c] = D.obs[c].astype(str).replace({"nan": "unknown", "": "unknown"}).values
    df["expression"], df["detected"] = expr, detected
    df["_row"] = np.arange(len(df))     # row in the H5AD, used by the DE step
    xy = D.embedding()
    if xy is not None:
        df["emb_1"], df["emb_2"] = xy[:, 0], xy[:, 1]
    all_groups = sorted(df[A.GROUP_COL].unique()) if A.GROUP_COL in df.columns else []
    subset = bool(groups) and A.GROUP_COL in df.columns and set(groups) != set(all_groups)
    keep = None
    if subset:
        keep = df[A.GROUP_COL].isin(groups).values
        df = df[keep].reset_index(drop=True)
        xy = xy[keep] if xy is not None else None
        pops = [p for p in pops if p in set(df[ct_col])]
        if df.empty:
            raise ValueError("No cells left after selecting groups.")
    study_dir = OUT_DIR / safe_name(study or default_study_name(fallback=D.title), 60)
    tag = ("__" + safe_name("_vs_".join(groups), 50)) if subset else ""
    out = study_dir / f"{safe_name(symbol)}__{safe_name(D.title, 40)}{tag}"
    out.mkdir(parents=True, exist_ok=True)
    m = dict(meta or {})
    m.update({"Gene": f"{symbol} ({var_id})", "Expression values": value_desc,
              "'Expressing' means": "count > 0 in that cell" + (" (from raw.X)" if D.has_raw else ""),
              "Annotations used": f"CELLxGENE {ct_col}, {donor_col}, {A.GROUP_COL}, {D.embedding_key}",
              "Groups analysed": (f"{A.GROUP_COL}: " + ", ".join(groups if subset else all_groups)) if all_groups else "all cells",
              "Cells": f"{len(df):,}"})
    meta_ = meta or {}
    lab = lambda c: ", ".join(sorted(df[c].unique())[:12]) if c in df.columns else "n/a"
    RUN_INFO.clear()
    RUN_INFO.update(gene=symbol, dataset=meta_.get("Dataset", D.title), collection=meta_.get("Collection", "n/a"),
                    doi=meta_.get("DOI", "n/a"), n_cells=len(df), n_celltypes=df[ct_col].nunique(),
                    n_subjects=df[donor_col].nunique() if donor_col else "n/a", ct_col=ct_col,
                    donor_col=donor_col or "n/a", embedding=D.embedding_key, has_raw=D.has_raw,
                    populations=pops, groups=list(groups) if subset else all_groups, groups_subset=subset,
                    tissues=lab("tissue"), diseases=lab("disease"),
                    assays=", ".join(sorted(D.obs["assay"].astype(str).unique())) if "assay" in D.obs.columns else "n/a")
    level_cellular(df, xy, D.embedding_key, symbol, ct_col, out / "1_cellular_level")
    if pops:
        level_cell_specific(df, xy, D.embedding_key, symbol, ct_col, pops, donor_col, out / "2_cell_specific_level",
                            min_cells, min_pct, min_group)
    if donor_col:
        level_sample(df, symbol, ct_col, donor_col, out / "3_sample_level", min_cells, min_pct, min_group)
    m["Generated"] = time.strftime("%Y-%m-%d %H:%M")
    (out / "METHODS.md").write_text("# Methods\n\n" + methods_markdown())
    write_study_info(out, D, symbol, var_id, study_dir.name, meta_, keep if subset else None, all_groups, pops)
    write_calc_page(out, symbol)
    write_report(out, m)
    return {"folder": out, "symbol": symbol, "df": df, "ct_col": ct_col, "donor_col": donor_col, "meta": m,
            "digest": build_digest(out)}


def _grid(paths, width=470):
    return HTML('<div style="display:flex;flex-wrap:wrap;gap:12px">' + "".join(
        f'<figure style="margin:0;width:{width}px"><img src="data:image/png;base64,{img_b64(p)}" style="width:100%">'
        f'<figcaption style="font-size:11px;color:#666">{cap}</figcaption></figure>' for cap, p in paths) + "</div>")


def interactive_cellular(res):
    """Interactive views at the top of the cellular-level tab (also saved as standalone HTML)."""
    df, sym, ct, out = res["df"], res["symbol"], res["ct_col"], res["folder"]
    display(HTML("<h3>Interactive views</h3><p><small>Hover for details, drag to zoom, double-click to reset. "
                 "Saved as standalone HTML in <code>interactive/</code>.</small></p>"))
    figs = {}
    if "emb_1" in df.columns:
        figs["umap_cell_types"] = IX.umap_celltype(df, ct, "Cell types (CELLxGENE annotation)")
        figs["umap_expression"] = IX.umap_expression(df, sym, ct)
    cs = pd.read_csv(out / "1_cellular_level" / "celltype_summary.csv")
    figs["bar_pct_expressing"] = IX.bar_celltype(cs, sym, name_col=cs.columns[0])
    for name, fig in figs.items():
        IX.save(fig, out / "interactive" / f"{name}.html")
        display(fig)


def show_results(res):
    out, sym = res["folder"], res["symbol"]
    tabs, names = [], ["1. Cellular level", "2. Cell-specific level", "3. Sample level"]
    for prefix in ["1.", "2.", "3."]:
        o = W.Output()
        with o:
            if prefix == "1.":
                interactive_cellular(res)
            secs = list(dict.fromkeys([s for s, *_ in FIGS + TABLES] + [s for s, _ in NOTES]))
            for s in [s for s in secs if s.startswith(prefix)]:
                display(HTML(f"<h3 style='margin-top:18px'>{s}</h3>"))
                for sec, txt in NOTES:
                    if sec == s:
                        display(HTML(f'<div style="background:#f3f6fa;padding:6px 10px;border-left:4px solid #4c72b0">{txt}</div>'))
                for sec, cap, path, prev in TABLES:
                    if sec == s:
                        display(HTML(f"<b>{cap}</b> <small>({path.name})</small>"))
                        display(prev)
                figs = [(cap, p) for sec, cap, p in FIGS if sec == s]
                if figs:
                    display(_grid(figs))
            if not any(s.startswith(prefix) for s in secs):
                display(HTML("<i>Not run (no population selected / no donor column).</i>"))
        tabs.append(o)
    co = W.Output()
    with co:
        display(HTML(f"<style>{CALC_CSS}</style>" + calc_guide_html(out, sym)
                     + f"<p><small>Saved as <code>{out / 'HOW_IT_IS_CALCULATED.html'}</code></small></p>"))
    tabs.append(co); names.append("How it's calculated")
    mo = W.Output()
    with mo:
        display(HTML('<div style="max-width:900px;line-height:1.5"><h3>Methods</h3>' + md_to_html(methods_markdown())
                     + f"<p><small>Saved as <code>{out / 'METHODS.md'}</code>.</small></p></div>"))
    tabs.append(mo); names.append("Methods")
    STATE["de_out"] = W.Output()
    with STATE["de_out"]:
        display(HTML("<p style='padding:12px'>Choose cell types and settings in section 5. of the panel, then press "
                     "<b>Run DE</b>. Fisher's exact tests run instantly; NEBULA needs R with the nebula package.</p>"))
    tabs.insert(3, STATE["de_out"]); names.insert(3, "4. DE (NEBULA + Fisher)")
    STATE["calc_out"] = co
    STATE["report_out"] = W.Output()
    with STATE["report_out"]:
        display(HTML("<p style='padding:12px'>Review the results first. When you're ready, press "
                     "<b>Generate descriptive report</b> in the panel above; "
                     f"<code>{llm_model_in.value}</code> will then write the report here.</p>"))
    tabs.append(STATE["report_out"]); names.append("Descriptive report")
    tab = W.Tab(children=tabs)
    for i, n in enumerate(names):
        tab.set_title(i, n)
    STATE["tab"] = tab
    ov = pd.read_csv(out / "1_cellular_level/overall_summary.csv").iloc[0]
    head = (f"<h2>{sym}</h2><p><b>{int(ov.n_cells_total):,}</b> cells in total · "
            f"<b>{int(ov.n_cells_expressing):,}</b> express {sym} "
            f"(<b>{ov.pct_expressing:.2f}%</b>, proportion {ov.proportion_expressing:.4f})")
    hp = out / "3_sample_level/subjects_expressing_headline.csv"
    if hp.exists():
        h = pd.read_csv(hp)
        head += "<br>Patients expressing: " + " · ".join(f"<b>{r.group}</b> {r.summary}" for r in h.itertuples())
    head += f'</p><p>Data report: <code>{out / "report.html"}</code></p>'
    with results:
        clear_output()
        display(HTML(head), tab)
    STATE["res"] = res
    report_btn.disabled = False
    d = res["df"]
    cnt = d.groupby(res["ct_col"]).detected.agg(["sum", "size"])
    de_ct_sel.options = [(f"{k}  ({int(r['sum']):,} + / {int(r['size'] - r['sum']):,} −)", k) for k, r in
                         cnt.sort_values("sum", ascending=False).iterrows()]
    de_ct_sel.value = tuple(k for k, r in cnt.iterrows()
                            if r["sum"] >= de_min.value and r["size"] - r["sum"] >= de_min.value)
    de_ct_sel.disabled = de_btn.disabled = False


def on_r_check():
    ok, msg = r_status(rscript_in.value.strip() or None)
    r_status_html.value = f"<small>{msg}</small>"
    return ok


def on_r_install():
    rs = rscript_in.value.strip()
    if not rs:
        return set_status("Enter the path to Rscript first.", "warn")
    set_status("Installing nebula in R (CRAN, or GitHub as fallback). This can take several minutes ...")
    ok, log_txt = r_install_nebula(rs)
    with log:
        print(log_txt)
    on_r_check()
    set_status("nebula installed." if ok else "nebula installation failed; see the Log.", "ok" if ok else "err")


def on_de():
    res = STATE.get("res")
    if not res:
        return set_status("Run the analysis first.", "warn")
    scopes = [s for s, w in [("pooled", de_pooled), ("within", de_within)] if w.value]
    if not scopes:
        return set_status("Tick at least one DE mode (pooled / within each group).", "warn")
    D = STATE["D"]
    on_r_check()
    gid = res["meta"]["Gene"].split("(")[-1].rstrip(")")
    STATE["tab"].selected_index = 3
    with STATE["de_out"]:
        clear_output()
        display(HTML("<p style='padding:12px'>Running ... progress is shown in the status line.</p>"))
    de_btn.disabled = True
    try:
        r = run_de(D, res["df"], res["symbol"], gid, res["ct_col"], res["donor_col"], A.GROUP_COL,
                   res["folder"] / "4_DE", celltypes=list(de_ct_sel.value) or None, scopes=scopes,
                   min_cells=de_min.value, max_cells=de_max.value, min_det=de_det.value / 100, fdr_thr=de_fdr.value,
                   rscript=rscript_in.value.strip() or None, ncore=de_cores.value,
                   progress=lambda m: set_status(f"{m} ..."))
    finally:
        de_btn.disabled = False
    nv = ""
    if "summary" in r and len(r["summary"]) and "nebula_version" in r["summary"].columns:
        nv = str(r["summary"]["nebula_version"].dropna().iloc[0]) if r["summary"]["nebula_version"].notna().any() else ""
    RUN_INFO["de"] = {"scopes": scopes, "min_cells": de_min.value, "min_donors": 3, "min_cells_per_donor": 3,
                      "max_cells": de_max.value, "min_det": de_det.value / 100, "fdr": de_fdr.value,
                      "nebula_version": nv}
    out = res["folder"]
    (out / "METHODS.md").write_text("# Methods\n\n" + methods_markdown())
    write_calc_page(out, res["symbol"])
    write_report(out, res["meta"])
    res["digest"] = build_digest(out)
    STATE["de_result"] = r
    show_de(r, res)
    with STATE["calc_out"]:
        clear_output()
        display(HTML(f"<style>{CALC_CSS}</style>" + calc_guide_html(out, res["symbol"])))
    msg = r.get("nebula_msg")
    if msg:
        set_status(f"Fisher done, but expression (NEBULA) did NOT run: {msg}. See the 4. DE tab for how to fix it.",
                   "warn")
    else:
        set_status(f"DE done: expression (NEBULA) and cell proportions (Fisher) saved in <code>{out / '4_DE'}</code>",
                   "ok")


def nebula_fix_hint(msg):
    m = (msg or "").lower()
    if "rscript not found" in m or "could not run rscript" in m:
        return ("Enter the path to Rscript in the <b>Rscript</b> box (Windows: <code>C:\\Program Files\\R\\R-4.4.x"
                "\\bin\\Rscript.exe</code>; Mac: <code>/usr/local/bin/Rscript</code> or "
                "<code>/Library/Frameworks/R.framework/Resources/bin/Rscript</code>), press <b>Check R / nebula</b>, "
                "then <b>Run DE</b> again.")
    if "not installed" in m or "test fit failed" in m:
        return ("Press <b>Install nebula</b> (it also updates parallelly/future), then <b>Check R / nebula</b> until it "
                "shows , then <b>Run DE</b> again. nebula 1.5.8 needs R ≥ 4.4.")
    if "no counts" in m or "integer counts" in m:
        return "This dataset has no raw counts on CELLxGENE, so NEBULA cannot be run on it. Fisher results are still valid."
    if "thresholds" in m:
        return "Lower <b>min cells/class</b>, or select cell types with more gene+ and gene- cells, then <b>Run DE</b> again."
    if "failed" in m:
        return "See the error in the NEBULA runs table and the Log; press <b>Check R / nebula</b> to test the R setup."
    return ""


def show_de(r, res):
    sec = DE_SECTION
    with STATE["de_out"]:
        clear_output()
        msg = r.get("nebula_msg")
        if msg:
            display(HTML(f'<div style="background:#fff3e0;border-left:6px solid #e65100;padding:10px 14px;margin:6px 0">'
                         f'<b>Expression (NEBULA) did not run: {msg}</b><br>{nebula_fix_hint(msg)}<br>'
                         '<small>The cell-proportion results (Fisher) below are complete.</small></div>'))
        else:
            n = len(r.get("jobs", []))
            display(HTML(f'<div style="background:#e8f5e9;border-left:6px solid #2e7d32;padding:8px 14px;margin:6px 0">'
                         f'Expression (NEBULA) ran for {n} cell type × mode combination(s); cell proportions '
                         '(Fisher) ran for all cell types.</div>'))
        display(HTML("<h2>A. Expression: NEBULA (gene+ vs gene− cells)</h2>"))
        for s_, txt in NOTES:
            if s_ == sec:
                display(HTML(f'<div style="background:#f3f6fa;padding:6px 10px;border-left:4px solid #4c72b0">{txt}</div>'))
        if "summary" in r and len(r["summary"]):
            display(HTML("<h3>NEBULA runs</h3>"))
            display(r["summary"])
        hm = [(cap, p) for s_, cap, p in FIGS if s_ == sec and p.name == "heatmap_top_de_genes.png"]
        if hm:
            display(_grid(hm, width=760))
        if r.get("jobs"):
            opts = [(f"{j['celltype']} [{j['scope']}]", i) for i, j in enumerate(r["jobs"])]
            pick = W.Dropdown(options=opts, description="Show", layout=L("520px"))
            box = W.Output()

            def render(_=None):
                j = r["jobs"][pick.value]
                with box:
                    clear_output()
                    fig = IX.volcano(j["table"], de_fdr.value, f"{j['celltype']} [{j['scope']}]: gene+ vs gene- cells")
                    IX.save(fig, j["folder"] / "volcano_interactive.html")
                    display(fig)
                    cols = ["gene", "log2FC", "p_value", "fdr", "pct_detected_pos", "pct_detected_neg",
                            "mean_lognorm_pos", "mean_lognorm_neg", "reliable_fit"]
                    display(HTML(f"<b>Top 30 genes</b> <small>(full table: {j['folder'] / 'nebula_de_results.csv'})</small>"))
                    display(j["table"][cols].head(30))
            pick.observe(render, "value")
            display(HTML("<h3>Per cell type</h3>"), pick, box)
            render()
        display(HTML("<h2>B. Cell proportions: Fisher's exact test</h2>"))
        if len(r["fisher_vs_rest"]):
            display(HTML("<b>Each cell type vs all other cells</b>"))
            display(r["fisher_vs_rest"].round(5))
        if len(r["fisher_between"]):
            display(HTML("<b>Between groups, within each cell type</b>"))
            display(r["fisher_between"].round(5))
        figs = [(cap, p) for s_, cap, p in FIGS if s_ == sec and p.name.startswith("fisher")]
        if figs:
            display(_grid(figs))


def on_report():
    res = STATE.get("res")
    if not res:
        return set_status("Run the analysis first.", "warn")
    model, url = llm_model_in.value.strip(), llm_url_in.value.strip()
    ok, msg = check_llm(model, url)
    if not ok:
        return set_status("" + msg, "warn")
    STATE["tab"].selected_index = len(STATE["tab"].children) - 1
    ro = STATE["report_out"]
    with ro:
        clear_output()
        display(HTML(f"<p style='padding:12px'>{model} is writing the report. This can take a few minutes "
                     "on a laptop.</p>"))
    report_btn.disabled = True
    try:
        body, css = generate_descriptive_report(
            res, model, url, progress=lambda m: set_status(f"{model}: {m} ..."))
    finally:
        report_btn.disabled = False
    scoped = css.replace("body{", ".rep{")
    with ro:
        clear_output()
        display(HTML(f"<style>{scoped}</style><div class='rep'>{body}</div>"))
    set_status(f"Report written: <code>{res['folder'] / 'DESCRIPTIVE_REPORT.html'}</code> "
               f"(and .md)", "ok")


search_btn.on_click(on_search); ds_btn.on_click(on_datasets); dl_btn.on_click(on_download)
local_btn.on_click(on_local); run_btn.on_click(on_run)
report_btn.on_click(guarded(on_report))
r_check_btn.on_click(guarded(on_r_check)); r_install_btn.on_click(guarded(on_r_install)); de_btn.on_click(guarded(on_de))

panel = W.VBox([
    W.HTML("<h3>1. Paper</h3>"), title_in, search_btn, coll_dd, ds_btn, ds_dd, dl_btn,
    W.HBox([local_in, local_btn]),
    study_in,
    W.HTML("<small>Results go to <code>OUT_DIR / study folder / GENE__dataset</code>. Filled in automatically "
           "(first author + year); type any name to override.</small>"),
    W.HTML("<h3>2. Groups to analyse</h3>"), group_col_dd, group_sel, W.HBox([grp_all_btn, group_info]),
    W.HTML("<small>All groups are selected by default. Ctrl/Cmd-click to keep only one or two (e.g. "
           "<i>normal</i> vs <i>COVID-19</i>). Only cells from the selected groups are analysed at every level.</small>"),
    W.HTML("<h3>3. Gene & population</h3>"), gene_in, ct_sel,
    W.HTML("<small>Hold Ctrl/Cmd to select several populations; none selected = skip the cell-specific level. "
           "A patient counts as <i>expressing</i> with ≥ <b>min cells</b> expressing cells and ≥ <b>min %</b> of its cells.</small>"),
    W.HBox([min_cells_in, min_pct_in]),
    run_btn,
    W.HTML("<h3>4. Differential expression: NEBULA (expression) + Fisher (proportions)</h3>"),
    W.HTML("<small>Within each cell type, compares cells with ≥1 count of your gene (+) against cells without it "
           "(−). Runs after the main analysis.</small>"),
    W.HBox([rscript_in, r_check_btn, r_install_btn]), r_status_html, de_ct_sel,
    W.HBox([de_pooled, de_within]), W.HBox([de_min, de_max, de_det]), W.HBox([de_fdr, de_cores]), de_btn,
    W.HTML("<h3>5. Descriptive report (optional, after you've reviewed the results)</h3>"),
    W.HBox([llm_model_in, llm_url_in]), report_btn, status,
    W.Accordion(children=[log]),
])
panel.children[-1].set_title(0, "Log (progress, details)")


def launch(out_dir="cellxgene_results", data_dir="cellxgene_data"):
    """Show the interactive panel. Results go to out_dir/<study>/<GENE>__<dataset>/."""
    global OUT_DIR, DATA_DIR
    OUT_DIR, DATA_DIR = Path(out_dir), Path(data_dir)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    display(panel, results)
    set_status("Paste a paper title above and press <b>Search CELLxGENE</b>.")
