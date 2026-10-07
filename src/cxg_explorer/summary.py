# cellxgeneExplorer
# Copyright (c) 2026 TeraWeave Technologies. Licensed under the PolyForm Noncommercial License 1.0.0 (see LICENSE.md).
# Required Notice: Copyright (c) 2026 TeraWeave Technologies
"""Cross-study summary engine: reads the outputs of the single-study explorer for many studies."""
import base64
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def safe_name(s, maxlen=60):
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_")
    return s[:maxlen] or "unnamed"


# ----------------------------------------------------------------------------- discovery
def discover_runs(root):
    """Every per-study run folder below root (identified by 1_cellular_level/overall_summary.csv)."""
    runs = []
    for oc in Path(root).rglob("1_cellular_level/overall_summary.csv"):
        run = oc.parent.parent
        if "_cross_study" in run.parts:
            continue
        info = {}
        sp = run / "study_info.json"
        if sp.exists():
            try:
                info = json.loads(sp.read_text())
            except Exception:
                info = {}
        parts = run.name.split("__")
        runs.append({
            "run": run, "gene": info.get("gene") or parts[0], "study": info.get("study") or run.parent.name,
            "dataset": info.get("dataset") or (parts[1] if len(parts) > 1 else run.name),
            "subset": bool(info.get("groups_subset", len(parts) > 2)),
            "groups": info.get("groups") or [], "info": info, "mtime": oc.stat().st_mtime,
            "label_run": run.name,
        })
    return runs


def study_keys(runs, gene):
    """Group runs of one gene by (study, dataset); default pick = all-groups run, newest first."""
    sel = [r for r in runs if r["gene"].upper() == gene.upper()]
    keys = {}
    for r in sel:
        keys.setdefault((r["study"], r["dataset"]), []).append(r)
    for k in keys:
        keys[k].sort(key=lambda r: (r["subset"], -r["mtime"]))
    n_ds = pd.Series([k[0] for k in keys]).value_counts().to_dict()
    labels = {k: (k[0] if n_ds[k[0]] == 1 else f"{k[0]} | {str(k[1])[:30]}") for k in keys}
    return keys, labels


# ----------------------------------------------------------------------------- cell type harmonization
BROAD_RULES = [
    ("Erythroid", ["erythro", "reticulocyte"]),
    ("Platelet / megakaryocyte", ["platelet", "megakaryocyte"]),
    ("Mast cell", ["mast cell"]),
    ("T cell", ["t cell", "t-cell", "thymocyte", "alpha-beta", "gamma-delta", "regulatory t", "mucosal invariant"]),
    ("NK / ILC", ["natural killer", "nk cell", "innate lymphoid", "nk t cell"]),
    ("B / plasma cell", ["b cell", "b-cell", "plasma cell", "plasmablast"]),
    ("Myeloid", ["monocyte", "macrophage", "dendritic", "neutrophil", "granulocyte", "myeloid", "basophil",
                 "eosinophil", "microglia", "kupffer", "langerhans", "phagocyte"]),
    ("Stem / progenitor", ["stem cell", "progenitor", "hematopoietic", "precursor"]),
    ("Endothelial", ["endothelial", "endothelium"]),
    ("Fibroblast / stromal", ["fibroblast", "stromal", "pericyte", "mesenchymal", "smooth muscle", "mesothelial",
                              "adipocyte", "chondrocyte", "osteoblast", "stellate"]),
    ("Epithelial", ["epithelial", "keratinocyte", "hepatocyte", "enterocyte", "goblet", "club", "ciliated", "basal",
                    "alveolar", "pneumocyte", "ionocyte", "secretory", "acinar", "duct", "cholangiocyte", "tuft",
                    "enteroendocrine", "paneth", "urothelial", "podocyte", "tubule", "luminal", "colonocyte",
                    "mucous", "serous", "squamous", "thyrocyte", "melanocyte", "kidney", "loop of henle",
                    "transit amplifying", "endocrine", "islet", "pancreatic a cell", "type a pancreatic",
                    "type b pancreatic", "pancreatic d cell", "pancreatic pp cell", "intestinal", "gastric",
                    "parietal", "chief cell", "mammary", "prostate", "bronchial", "tracheal", "nasal"]),
    ("Neural / glial", ["neuron", "astrocyte", "oligodendrocyte", "glial", "schwann", "neural", "photoreceptor",
                        "ependymal", "interneuron"]),
    ("Muscle", ["muscle", "cardiomyocyte", "myocyte", "myoblast", "satellite cell"]),
    ("Germ / reproductive", ["sperm", "oocyte", "germ", "granulosa", "sertoli", "leydig", "trophoblast",
                             "decidual"]),
]


def auto_broad_class(label):
    t = str(label).lower()
    for cls, kws in BROAD_RULES:
        if any(re.search(r"\b" + re.escape(k), t) for k in kws):
            return cls
    return "Other / unassigned"


def load_or_update_mapping(path, labels):
    """Editable CSV: cell_type -> broad_class. Existing rows (user edits) are kept; new labels are auto-classified."""
    path = Path(path)
    old = pd.read_csv(path) if path.exists() else pd.DataFrame(columns=["cell_type", "broad_class"])
    known = set(old["cell_type"].astype(str))
    new = [l for l in sorted(set(labels)) if l not in known]
    if new:
        old = pd.concat([old, pd.DataFrame({"cell_type": new, "broad_class": [auto_broad_class(l) for l in new]})],
                        ignore_index=True)
        path.parent.mkdir(parents=True, exist_ok=True)
        old.sort_values(["broad_class", "cell_type"]).to_csv(path, index=False)
    return dict(zip(old["cell_type"].astype(str), old["broad_class"].astype(str))), len(new)


# ----------------------------------------------------------------------------- loading one study
def wilson_ci(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h) * 100, min(1.0, c + h) * 100


def first_key(d):
    return ", ".join(f"{k}" for k, _ in sorted(d.items(), key=lambda x: -x[1])[:3]) if d else "n/a"


def load_study(r, label, min_cells, min_pct):
    run, info = r["run"], r["info"]
    c1, c3 = run / "1_cellular_level", run / "3_sample_level"
    ov = pd.read_csv(c1 / "overall_summary.csv").iloc[0]
    ct = pd.read_csv(c1 / "celltype_summary.csv")
    ct = ct.rename(columns={ct.columns[0]: "cell_type"})[["cell_type", "n_cells", "n_expressing"]]
    ct["cell_type"] = ct["cell_type"].astype(str)
    ct["study"] = label
    grp_col = info.get("group_col") or "disease"

    tissue = first_key(info.get("tissue", {}))
    if tissue == "n/a" and (c1 / "summary_by_tissue.csv").exists():
        tissue = ", ".join(pd.read_csv(c1 / "summary_by_tissue.csv").iloc[:3, 0].astype(str))

    donors, groups = pd.DataFrame(), pd.DataFrame()
    ss = c3 / "subject_summary.csv"
    if ss.exists():
        s = pd.read_csv(ss)
        dcol = info.get("donor_col") if info.get("donor_col") in s.columns else s.columns[0]
        s[dcol] = s[dcol].astype(str)
        d = s.groupby(dcol)[["n_cells", "n_expressing"]].sum().reset_index().rename(columns={dcol: "donor"})
        d["pct_expressing"] = 100 * d.n_expressing / d.n_cells
        d["expressing"] = (d.n_expressing >= min_cells) & (d.pct_expressing >= min_pct)
        d["study"] = label
        donors = d
        if grp_col in s.columns:
            u = s.copy()
            u["expressing"] = (u.n_expressing >= min_cells) & (u.pct_expressing >= min_pct)
            g = u.groupby(grp_col).agg(n_cells=("n_cells", "sum"), n_expressing=("n_expressing", "sum"),
                                       n_subjects=(dcol, "nunique"),
                                       n_subjects_expressing=("expressing", "sum")).reset_index()
            g = g.rename(columns={grp_col: "group"})
            g["study"] = label
            groups = g
    if groups.empty and (c1 / f"summary_by_{grp_col}.csv").exists():
        g = pd.read_csv(c1 / f"summary_by_{grp_col}.csv")
        groups = g.rename(columns={g.columns[0]: "group"})[["group", "n_cells", "n_expressing"]].assign(study=label)

    if groups.empty:   # single-group study: one row with the study's totals
        gname = (info.get("groups") or info.get("all_groups") or ["all cells"])[0]
        groups = pd.DataFrame([{"group": str(gname), "n_cells": int(ov.n_cells_total),
                                "n_expressing": int(ov.n_cells_expressing), "study": label}])
        if not donors.empty:
            groups["n_subjects"] = len(donors)
            groups["n_subjects_expressing"] = int(donors.expressing.sum())
    n_pat = len(donors) if not donors.empty else np.nan
    k_pat = int(donors.expressing.sum()) if not donors.empty else np.nan
    prim = info.get("is_primary_fraction")
    inv = {
        "study": label, "study_folder": r["study"], "dataset": r["dataset"], "tissue": tissue,
        "assay": first_key(info.get("assay", {})), "organism": first_key(info.get("organism", {})),
        "group_column": grp_col, "groups": ", ".join(map(str, info.get("groups") or (groups["group"].astype(str).tolist()
                                                                                   if not groups.empty else []))),
        "n_cells": int(ov.n_cells_total), "n_cells_expressing": int(ov.n_cells_expressing),
        "pct_cells_expressing": round(100 * ov.n_cells_expressing / max(ov.n_cells_total, 1), 3),
        "n_patients": n_pat, "n_patients_expressing": k_pat,
        "pct_patients_expressing": round(100 * k_pat / n_pat, 1) if n_pat and n_pat == n_pat else np.nan,
        "pct_primary_cells": round(100 * prim, 1) if prim is not None else np.nan,
        "doi": info.get("doi", ""), "run_folder": str(run),
    }
    return inv, ct, donors, groups


# ----------------------------------------------------------------------------- summary tables
def build_summary(selected, labels, gene, min_cells=3, min_pct=1.0, min_ct_cells=20, expr_thr=10.0,
                  mapping_path=None, online_citations=True):
    inv, cts, dons, grps = [], [], [], []
    for key, r in selected.items():
        a, b, c, d = load_study(r, labels[key], min_cells, min_pct)
        inv.append(a); cts.append(b)
        if not c.empty:
            dons.append(c)
        if not d.empty:
            grps.append(d)
    inv = pd.DataFrame(inv)
    ct = pd.concat(cts, ignore_index=True)
    mapping, n_new = load_or_update_mapping(mapping_path, ct["cell_type"].unique())
    ct["broad_class"] = ct["cell_type"].map(mapping).fillna("Other / unassigned")
    ct["pct_expressing"] = 100 * ct.n_expressing / ct.n_cells
    ct["enough_cells"] = ct.n_cells >= min_ct_cells

    # top cell type per study (among cell types with enough cells)
    ok = ct[ct.enough_cells]
    top = ok.sort_values(["study", "pct_expressing", "n_expressing"], ascending=[True, False, False]) \
        .groupby("study").head(1).set_index("study")
    inv["top_cell_type"] = inv.study.map(top["cell_type"])
    inv["top_cell_type_pct"] = inv.study.map(top["pct_expressing"].round(2))
    inv["patients_expressing"] = inv.apply(
        lambda x: f"{int(x.n_patients_expressing)}/{int(x.n_patients)} ({x.pct_patients_expressing:.1f}%)"
        if x.n_patients == x.n_patients else "n/a", axis=1)
    order = inv.sort_values(["tissue", "study"]).study.tolist()

    def matrix(level):
        if level == "fine":
            m = ct.pivot_table(index="study", columns="cell_type", values="pct_expressing")
            n = ct.pivot_table(index="study", columns="cell_type", values="n_cells")
        else:
            agg = ct.groupby(["study", "broad_class"])[["n_cells", "n_expressing"]].sum().reset_index()
            agg["pct"] = 100 * agg.n_expressing / agg.n_cells
            m = agg.pivot(index="study", columns="broad_class", values="pct")
            n = agg.pivot(index="study", columns="broad_class", values="n_cells")
        return m.where(n >= min_ct_cells).reindex(order), n.reindex(order)

    fine_pct, fine_n = matrix("fine")
    broad_pct, broad_n = matrix("broad")

    def consistency(pct, level_name):
        top_counts = pct.idxmax(axis=1).value_counts() if pct.notna().any().any() else pd.Series(dtype=int)
        rows = []
        for c in pct.columns:
            v = pct[c].dropna()
            if v.empty:
                continue
            rows.append({level_name: c, "n_studies_with_celltype": len(v),
                         "n_studies_>=1%": int((v >= 1).sum()), "n_studies_>=5%": int((v >= 5).sum()),
                         f"n_studies_>={expr_thr:g}%": int((v >= expr_thr).sum()), "n_studies_>=25%": int((v >= 25).sum()),
                         "n_studies_top_celltype": int(top_counts.get(c, 0)),
                         "median_pct": round(float(v.median()), 2), "min_pct": round(float(v.min()), 2),
                         "max_pct": round(float(v.max()), 2),
                         "summary": f"≥{expr_thr:g}% of cells in {int((v >= expr_thr).sum())} of {len(v)} studies"})
        out = pd.DataFrame(rows)
        if not out.empty:
            out = out.sort_values([f"n_studies_>={expr_thr:g}%", "median_pct"], ascending=False)
        return out

    cons_fine = consistency(fine_pct, "cell_type")
    cons_broad = consistency(broad_pct, "broad_class")
    donors = pd.concat(dons, ignore_index=True) if dons else pd.DataFrame()
    groups = pd.concat(grps, ignore_index=True) if grps else pd.DataFrame()
    if not groups.empty:
        groups["pct_cells_expressing"] = (100 * groups.n_expressing / groups.n_cells).round(3)
        if "n_subjects" in groups.columns:
            groups["pct_subjects_expressing"] = (100 * groups.n_subjects_expressing / groups.n_subjects).round(1)
            groups["subjects_expressing"] = groups.apply(
                lambda x: f"{int(x.n_subjects_expressing)}/{int(x.n_subjects)} ({x.pct_subjects_expressing:.1f}%)",
                axis=1)
    pat = inv[["study", "tissue", "n_patients", "n_patients_expressing", "pct_patients_expressing"]].dropna().copy()
    if not pat.empty:
        ci = pat.apply(lambda x: wilson_ci(x.n_patients_expressing, x.n_patients), axis=1)
        pat["ci95_low"], pat["ci95_high"] = [round(c[0], 1) for c in ci], [round(c[1], 1) for c in ci]

    warnings = []
    lowprim = inv[inv.pct_primary_cells < 100]
    for x in lowprim.itertuples():
        warnings.append(f"<b>{x.study}</b>: only {x.pct_primary_cells:.0f}% of cells are original to this dataset "
                        "(the rest are re-used from other CELLxGENE datasets), so the same cells or patients may "
                        "appear in another study.")
    if not donors.empty:
        dd = donors.groupby("donor")["study"].nunique()
        shared = dd[dd > 1].index
        if len(shared):
            pairs = donors[donors.donor.isin(shared)].groupby("donor")["study"].apply(lambda s: " & ".join(sorted(set(s))))
            warnings.append(f"{len(shared)} donor ID(s) appear in more than one study (possible overlap; or just "
                            f"generic IDs like 'donor1'): " + "; ".join(f"{d} ({s})" for d, s in pairs.head(10).items()))
    if inv.n_patients.isna().any():
        warnings.append("Some studies have no subject-level table, so they are missing from the patient summaries.")
    if inv.groups.str.len().gt(0).any() and selected:
        sub = [labels[k] for k, r in selected.items() if r["subset"]]
        if sub:
            warnings.append("These studies use a run restricted to selected groups: " + ", ".join(sub))
    if n_new:
        warnings.append(f"{n_new} new cell type label(s) were auto-assigned to broad classes. Check/edit "
                        f"<code>{mapping_path}</code> and press Build again to apply your edits.")

    det, samples, ed_path, n_new_ed = build_study_details(selected, labels, Path(mapping_path).parent,
                                                          online=online_citations)
    det = det.set_index("study").loc[order].reset_index()
    det["ref"] = range(1, len(det) + 1)
    refmap = dict(zip(det.study, det.ref))
    inv.insert(0, "ref", inv.study.map(refmap))
    if n_new_ed:
        warnings.append(f"Added {n_new_ed} study row(s) to <code>{ed_path}</code>: fill in how samples were obtained "
                        "(and any notes or citation corrections), then press Build again.")
    missing = det[det.citation.str.startswith("citation not available")].study.tolist()
    if missing:
        warnings.append("No citation information for: " + ", ".join(missing) + ". Re-run these studies with the "
                        "current agent notebook, or add a citation in the editable study sheet.")
    return {"selected": selected, "labels": labels, "mapping": mapping,
            "details": det, "samples_used": samples, "details_editable": ed_path,
            "gene": gene, "inventory": inv.set_index("study").loc[order].reset_index(), "celltype_long": ct,
            "fine_pct": fine_pct, "fine_n": fine_n, "broad_pct": broad_pct, "broad_n": broad_n,
            "consistency_fine": cons_fine, "consistency_broad": cons_broad, "donors": donors, "groups": groups,
            "patients": pat.set_index("study").reindex([o for o in order if o in set(pat.study)]).reset_index()
            if not pat.empty else pat, "order": order, "warnings": warnings,
            "params": {"min_cells": min_cells, "min_pct": min_pct, "min_ct_cells": min_ct_cells,
                       "expr_thr": expr_thr, "mapping": str(mapping_path)}}


# ----------------------------------------------------------------------------- study details (citation, samples, technology)
import requests

CROSSREF = "https://api.crossref.org/works/"


def _authors_short(auth, n=6):
    names = []
    for a in auth:
        fam = a.get("family") or (a.get("name") or "").split(" ")[-1]
        giv = "".join(p[0] for p in re.split(r"[\s\-]+", a.get("given") or "") if p)
        if fam:
            names.append(f"{fam} {giv}".strip())
    if not names:
        return ""
    return ", ".join(names[:n]) + (", et al." if len(names) > n else "")


def fetch_citation(doi, cache, online=True):
    """Full reference for a DOI via Crossref (cached). Returns dict or None."""
    if not doi:
        return None
    doi = re.sub(r"^(https?://)?(dx\.)?doi\.org/", "", str(doi)).strip().rstrip(".")
    if doi in cache:
        return cache[doi]
    if not online:
        return None
    try:
        r = requests.get(CROSSREF + doi, timeout=20, headers={"User-Agent": "cellxgene-cross-study/1.0"})
        r.raise_for_status()
        m = r.json()["message"]
        yr = (m.get("published-print") or m.get("published-online") or m.get("issued") or {}).get("date-parts", [[None]])[0][0]
        cache[doi] = {"authors": _authors_short(m.get("author") or []), "title": (m.get("title") or [""])[0],
                      "journal": (m.get("container-title") or [""])[0], "year": yr, "volume": m.get("volume"),
                      "issue": m.get("issue"), "pages": m.get("page") or m.get("article-number"), "doi": doi}
    except Exception:
        cache[doi] = None
    return cache[doi]


def format_citation(info, cr):
    doi = info.get("doi")
    if cr:
        vol = f";{cr['volume']}" + (f"({cr['issue']})" if cr.get("issue") else "") if cr.get("volume") else ""
        pg = f":{cr['pages']}" if cr.get("pages") else ""
        return f"{cr['authors']}. {cr['title']}. {cr['journal']}. {cr['year']}{vol}{pg}. doi:{cr['doi']}"
    c = info.get("citation") or {}
    au = _authors_short(c.get("authors") or [], n=3)
    parts = [p for p in [au, f"({c.get('year')})" if c.get("year") else "", info.get("collection"),
                         c.get("journal"), f"doi:{doi}" if doi else ""] if p]
    return ". ".join(parts) if parts else "citation not available (re-run this study with the current agent notebook)"


def _counts_txt(d, n=4):
    if not d:
        return "n/a"
    tot = sum(d.values())
    items = sorted(d.items(), key=lambda x: -x[1])
    s = "; ".join(f"{k} ({100 * v / tot:.0f}%)" for k, v in items[:n])
    return s + (f"; +{len(items) - n} more" if len(items) > n else "")


def _links_txt(links):
    keep = [l for l in links or [] if l.get("url") and str(l.get("type", "")).upper() != "DOI"]
    return " | ".join(f"{(l.get('type') or 'link').replace('_', ' ').title()}: {l.get('name') or l['url']} <{l['url']}>"
                      for l in keep)


def _acquisition_auto(info):
    parts = []
    tt = info.get("tissue_type") or {}
    if tt:
        parts.append("CELLxGENE tissue type: " + _counts_txt(tt, 3))
    for col, vals in (info.get("author_sample_annotations") or {}).items():
        parts.append(f"{col}: " + ", ".join(f"{k} ({v:,} cells)" for k, v in list(vals.items())[:6]))
    return " · ".join(parts)


def build_study_details(selected, labels, out_root, online=True, desc_chars=700):
    cache_p = Path(out_root) / "citations_cache.json"
    cache = json.loads(cache_p.read_text()) if cache_p.exists() else {}
    rows, samples = [], []
    for key, r in selected.items():
        info, run, lab = r["info"], r["run"], labels[key]
        cr = fetch_citation(info.get("doi"), cache, online)
        desc = (info.get("collection_description") or "").strip()
        groups = info.get("groups") or []
        sp = run / "samples_used.csv"
        su = pd.read_csv(sp) if sp.exists() else pd.DataFrame()
        if not su.empty:
            su.insert(0, "study", lab)
            samples.append(su)
        gcol = info.get("group_col") or "disease"
        dcol = info.get("donor_col")
        if not su.empty and dcol in su.columns:
            per_g = (su.groupby(gcol)[dcol].nunique().to_dict() if gcol in su.columns else {})
            donors_txt = f"{su[dcol].nunique()} donors" + (" (" + ", ".join(f"{k}: {v}" for k, v in per_g.items()) + ")"
                                                         if len(per_g) > 1 else "")
        else:
            donors_txt = f"{info.get('n_subjects', 'n/a')} donors"
        if info.get("n_samples"):
            donors_txt += f"; {info['n_samples']} samples (column '{info.get('sample_column')}')"
        susp = info.get("suspension_type") or {}
        tech = _counts_txt(info.get("assay"), 3)
        if susp:
            tech += " · " + ", ".join(f"single-{'nucleus' if k == 'nucleus' else 'cell' if k == 'cell' else k}"
                                      for k in susp)
        rows.append({
            "ref": None, "study": lab, "citation": format_citation(info, cr),
            "doi_url": f"https://doi.org/{info['doi']}" if info.get("doi") else "",
            "collection_url": f"https://cellxgene.cziscience.com/collections/{info['collection_id']}"
            if info.get("collection_id") else "",
            "dataset": info.get("dataset") or r["dataset"],
            "technology": tech,
            "raw_counts_available": info.get("has_raw_counts"),
            "organism": _counts_txt(info.get("organism"), 2),
            "tissue": _counts_txt(info.get("tissue"), 4),
            "groups_used": f"{gcol}: " + ", ".join(map(str, groups)) + (" (subset)" if info.get("groups_subset") else "")
            if groups else "all cells",
            "donors_samples_used": donors_txt,
            "cells_used": info.get("n_cells"),
            "sex": _counts_txt(info.get("sex"), 3),
            "development_stage": _counts_txt(info.get("development_stage"), 3),
            "sample_acquisition_auto": _acquisition_auto(info),
            "author_description": desc[:desc_chars] + ("…" if len(desc) > desc_chars else ""),
            "author_description_full": desc,
            "data_links": _links_txt(info.get("links")),
            "consortia": ", ".join(info.get("consortia") or []),
        })
    cache_p.parent.mkdir(parents=True, exist_ok=True)
    cache_p.write_text(json.dumps(cache, indent=1))
    det = pd.DataFrame(rows)

    # merge with the user's editable sheet (their entries win)
    ed_p = Path(out_root) / "study_details_editable.csv"
    ed_cols = ["study", "citation_override", "sample_acquisition", "notes"]
    ed = pd.read_csv(ed_p, dtype=str).fillna("") if ed_p.exists() else pd.DataFrame(columns=ed_cols)
    for c in ed_cols:
        if c not in ed.columns:
            ed[c] = ""
    new = [s for s in det.study if s not in set(ed.study)]
    if new:
        ed = pd.concat([ed, pd.DataFrame({"study": new, "citation_override": "", "sample_acquisition": "",
                                          "notes": ""})], ignore_index=True)
        ed.sort_values("study").to_csv(ed_p, index=False)
    e = ed.set_index("study")
    det["citation"] = [e.at[s, "citation_override"] if s in e.index and e.at[s, "citation_override"] else c
                       for s, c in zip(det.study, det.citation)]
    det["sample_acquisition"] = [e.at[s, "sample_acquisition"] if s in e.index and e.at[s, "sample_acquisition"]
                                 else a for s, a in zip(det.study, det.sample_acquisition_auto)]
    det["sample_acquisition_source"] = ["your edit" if s in e.index and e.at[s, "sample_acquisition"] else
                                        ("CELLxGENE annotations" if a else "not available") for s, a in
                                        zip(det.study, det.sample_acquisition_auto)]
    det["notes"] = [e.at[s, "notes"] if s in e.index else "" for s in det.study]
    return det, (pd.concat(samples, ignore_index=True) if samples else pd.DataFrame()), str(ed_p), len(new)


def details_html(det):
    """One card per study: citation, samples used, acquisition, technology, links."""
    esc = lambda s: str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    cards = []
    for x in det.itertuples():
        link = lambda u, t: f"<a href='{u}' target='_blank'>{t}</a>" if u else ""
        links = esc(x.data_links)
        links = re.sub(r"&lt;(https?://[^&]+)&gt;", r"<a href='\1' target='_blank'>link</a>", links)
        acq = esc(x.sample_acquisition) or ("<i>not described in CELLxGENE annotations; see the authors' "
                                            "description and the paper</i>")
        short = esc(x.author_description[:220]) or "<i>none</i>"
        more = "…" if len(x.author_description) > 220 else ""
        raw = " · raw counts available" if x.raw_counts_available else ""
        notes = f"<tr><td>Notes</td><td>{esc(x.notes)}</td></tr>" if x.notes else ""
        cards.append(f"""<div class='card'><div class='ct'>[{x.ref}] {esc(x.study)}</div>
<p class='cit'>{esc(x.citation)} {link(x.doi_url, 'DOI')} {link(x.collection_url, 'CELLxGENE')}</p>
<table class='kv'>
<tr><td>Dataset</td><td>{esc(x.dataset)}</td></tr>
<tr><td>Technology</td><td>{esc(x.technology)}{raw}</td></tr>
<tr><td>Samples used in this analysis</td><td>{esc(x.groups_used)}<br>{esc(x.donors_samples_used)} · {x.cells_used:,} cells</td></tr>
<tr><td>Tissue / organism</td><td>{esc(x.tissue)} · {esc(x.organism)}</td></tr>
<tr><td>Donors</td><td>sex: {esc(x.sex)}<br>stage: {esc(x.development_stage)}</td></tr>
<tr><td>How samples were obtained</td><td>{acq}
<br><small>source: {x.sample_acquisition_source}</small></td></tr>
<tr><td>Authors' description (CELLxGENE)</td><td><details><summary>{short}{more}</summary>{esc(x.author_description_full)}</details></td></tr>
<tr><td>Data / protocol links</td><td>{links or '—'}</td></tr>
{notes}
</table></div>""")
    return "\n".join(cards)


DETAILS_CSS = """
.card{background:#fff;border:1px solid #e3e3e3;border-radius:6px;padding:12px 16px;margin:10px 0;max-width:1100px}
.card .ct{font-weight:600;font-size:15px}.card .cit{margin:4px 0 8px 0;font-size:13px}
table.kv{border-collapse:collapse;font-size:13px;width:100%}table.kv td{border-top:1px solid #f0f0f0;padding:4px 8px;
vertical-align:top}table.kv td:first-child{width:220px;color:#555;white-space:nowrap}
"""


# ----------------------------------------------------------------------------- plots
FIGS = []


def _save(fig, path, caption, tab):
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    FIGS.append((tab, caption, Path(path)))


def _palette(n):
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    return [base[i % len(base)] for i in range(n)]


def plot_heatmap(pct, path, title, caption, tab, max_cols=60, annotate_max=25,
                 note="grey = cell type absent or too few cells in that study"):
    if pct.empty or pct.shape[1] == 0:
        return
    m = pct.copy()
    if m.shape[1] > max_cols:
        keep = m.notna().sum().sort_values(ascending=False).index[:max_cols]
        m = m[keep]
    m = m[m.median().sort_values(ascending=False).index]
    r, c = m.shape
    fig, ax = plt.subplots(figsize=(min(32, max(6, 0.36 * c + 4)), min(30, max(3, 0.34 * r + 2))))
    cmap = plt.get_cmap("magma_r").copy()
    cmap.set_bad("#eeeeee")
    im = ax.imshow(np.ma.masked_invalid(m.values.astype(float)), aspect="auto", cmap=cmap, vmin=0,
                   vmax=max(1, np.nanpercentile(m.values.astype(float), 98)) if np.isfinite(m.values.astype(float)).any() else 1)
    ax.set_xticks(range(c)); ax.set_xticklabels([str(x)[:40] for x in m.columns], rotation=70, ha="right", fontsize=7)
    ax.set_yticks(range(r)); ax.set_yticklabels([str(x)[:40] for x in m.index], fontsize=8)
    if r <= annotate_max and c <= annotate_max:
        vmax = im.get_clim()[1]
        for i in range(r):
            for j in range(c):
                v = m.values[i, j]
                if v == v:
                    ax.text(j, i, f"{v:.0f}" if v >= 1 else f"{v:.1f}", ha="center", va="center", fontsize=6,
                            color="white" if v > 0.6 * vmax else "black")
    cb = fig.colorbar(im, ax=ax, shrink=0.5)
    cb.set_label("% cells expressing")
    ax.set_title(title + f"\n({note})", fontsize=10)
    _save(fig, path, caption, tab)


def plot_consistency(cons, name_col, thr, path, title, caption, tab, max_items=40):
    if cons.empty:
        return
    d = cons.head(max_items)
    h = max(2.5, 0.3 * len(d) + 1)
    fig, ax = plt.subplots(figsize=(8.5, h))
    y = np.arange(len(d))
    ax.barh(y, d["n_studies_with_celltype"], color="#dddddd", label="studies with this cell type")
    ax.barh(y, d[f"n_studies_>={thr:g}%"], color="#c0392b", label=f"studies with ≥{thr:g}% cells expressing")
    for yi, (a, b) in enumerate(zip(d[f"n_studies_>={thr:g}%"], d["n_studies_with_celltype"])):
        ax.text(b + 0.1, yi, f"{a}/{b}", va="center", fontsize=7)
    ax.set_yticks(y); ax.set_yticklabels([str(x)[:50] for x in d[name_col]], fontsize=8)
    ax.invert_yaxis(); ax.set_xlabel("number of studies")
    ax.legend(frameon=False, fontsize=7, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    ax.set_title(title, fontsize=10, pad=22)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    _save(fig, path, caption, tab)


def plot_study_overall(inv, path, gene, tab):
    d = inv.sort_values("pct_cells_expressing", ascending=True)
    tissues = d.tissue.unique().tolist()
    cmap = dict(zip(tissues, _palette(len(tissues))))
    fig, ax = plt.subplots(figsize=(8.5, max(2.5, 0.32 * len(d) + 1)))
    ax.barh(range(len(d)), d.pct_cells_expressing, color=[cmap[t] for t in d.tissue])
    for i, x in enumerate(d.itertuples()):
        ax.text(x.pct_cells_expressing, i, f"  {x.n_cells_expressing:,}/{x.n_cells:,}", va="center", fontsize=7)
    ax.set_yticks(range(len(d))); ax.set_yticklabels(d.study, fontsize=8)
    ax.set_xlabel(f"% of all cells expressing {gene}")
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=cmap[t], label=t[:35]) for t in tissues],
              title="tissue", fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.set_title(f"Overall detection of {gene} per study", fontsize=10)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    _save(fig, path, "Overall % of cells expressing per study (label: expressing/total cells), coloured by tissue", tab)


def plot_donor_strip(donors, order, path, gene, tab):
    if donors.empty:
        return
    order = [o for o in order if o in set(donors.study)]
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(order) + 3), 4.8))
    rng = np.random.default_rng(0)
    for i, s in enumerate(order):
        v = donors.loc[donors.study == s]
        x = i + rng.uniform(-0.18, 0.18, len(v))
        ax.scatter(x, v.pct_expressing, s=16, alpha=0.75,
                   c=np.where(v.expressing, "#c0392b", "#999999"), zorder=3)
        ax.hlines(v.pct_expressing.median(), i - 0.28, i + 0.28, color="black", lw=1.5)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel(f"% of a patient's cells expressing {gene}")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#c0392b", label="expressing patient"),
                       Line2D([], [], marker="o", ls="", color="#999999", label="not expressing"),
                       Line2D([], [], color="black", label="median")], frameon=False, fontsize=7)
    ax.set_title(f"Patient-level {gene} detection in each study (one dot = one patient)", fontsize=10)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    _save(fig, path, "Each patient's % expressing cells, by study; black line = study median", tab)


def plot_groups(groups, path, gene, tab):
    if groups.empty or groups.group.nunique() > 40:
        return
    m = groups.pivot_table(index="study", columns="group", values="pct_cells_expressing")
    plot_heatmap(m, path, f"% cells expressing {gene}: study × group", "Detection per group within each study "
                 "(compare groups within a row; rows are different studies)", tab, annotate_max=40,
                 note="grey = group not present in that study")


# ----------------------------------------------------------------------------- cell-type-specific cross-study report
CT_FIGS = []


def _save_ct(fig, path, caption, key):
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    CT_FIGS.append((key, caption, Path(path)))


def celltype_choices(S):
    """[(label, value)] for the selector: broad classes first, then CELLxGENE cell types, with study counts."""
    ct = S["celltype_long"]
    fine = ct[ct.enough_cells].groupby("cell_type").study.nunique().sort_values(ascending=False)
    broad = S["broad_pct"].notna().sum().sort_values(ascending=False)
    pl = lambda n: f"{int(n)} study" if int(n) == 1 else f"{int(n)} studies"
    opts = [(f"[class] {c}  ({pl(n)})", ("broad", c)) for c, n in broad.items() if n > 0]
    opts += [(f"{c}  ({pl(n)})", ("fine", c)) for c, n in fine.items()]
    return opts


_PC_CACHE = {}


def load_per_cell(run, info):
    run = Path(run)
    if run in _PC_CACHE:
        return _PC_CACHE[run]
    p = run / "1_cellular_level" / "per_cell_expression.csv"
    if not p.exists():
        _PC_CACHE[run] = None
        return None
    head = pd.read_csv(p, nrows=0).columns.tolist()
    ctc = info.get("celltype_col") if info.get("celltype_col") in head else head[1]
    dc = info.get("donor_col") if info.get("donor_col") in head else ("donor_id" if "donor_id" in head else None)
    gc = info.get("group_col") if info.get("group_col") in head else ("disease" if "disease" in head else None)
    cols = [c for c in dict.fromkeys([ctc, dc, gc, "detected"]) if c]
    d = pd.read_csv(p, usecols=cols, dtype={c: str for c in cols if c != "detected"})
    d = d.rename(columns={ctc: "cell_type", **({dc: "donor"} if dc else {}), **({gc: "group"} if gc else {})})
    d["detected"] = d["detected"].astype(str).str.lower().isin(["true", "1"])
    _PC_CACHE[run] = d
    return d


def build_celltype_summary(S, target, min_cells=3, min_pct=1.0, min_ct_cells=20, min_pat_cells=10, fdr_thr=0.05):
    kind, name = target
    mapping, labels, selected = S["mapping"], S["labels"], S["selected"]
    members_all = [c for c, b in mapping.items() if b == name] if kind == "broad" else [name]
    ov, donors, groups, f_rest, f_btw, de = [], [], [], [], [], []
    ctl = S["celltype_long"]
    for key, r in selected.items():
        lab, info = labels[key], r["info"]
        pc = load_per_cell(r["run"], info)
        st = ctl[ctl.study == lab]
        # rank of this cell type / class within the study (by % expressing, among those with enough cells)
        if kind == "fine":
            rk = st[st.enough_cells].sort_values("pct_expressing", ascending=False).reset_index(drop=True)
            pos = rk.index[rk.cell_type == name]
            n_ranked = len(rk)
        else:
            b = st.groupby("broad_class")[["n_cells", "n_expressing"]].sum()
            b = b[b.n_cells >= min_ct_cells]
            b["pct"] = b.n_expressing / b.n_cells
            rk = b.sort_values("pct", ascending=False).reset_index()
            pos = rk.index[rk.broad_class == name]
            n_ranked = len(rk)
        rank = int(pos[0]) + 1 if len(pos) else None
        if pc is not None:
            sub = pc[pc.cell_type.isin(members_all)]
        else:
            sub = None
        n = int(len(sub)) if sub is not None else int(st[st.cell_type.isin(members_all)].n_cells.sum())
        k = int(sub.detected.sum()) if sub is not None else int(st[st.cell_type.isin(members_all)].n_expressing.sum())
        row = {"study": lab, "tissue": S["inventory"].set_index("study").tissue.get(lab, ""),
               "members_in_study": ", ".join(sorted(set(sub.cell_type))) if sub is not None and len(sub) else "",
               "n_cells": n, "n_expressing": k, "pct_expressing": round(100 * k / n, 3) if n else np.nan,
               "present": n >= min_ct_cells,
               "rank_in_study": f"{rank} of {n_ranked}" if rank else "",
               "n_patients_eligible": np.nan, "n_patients_expressing": np.nan}
        if sub is not None and "donor" in sub.columns and n:
            d = sub.groupby("donor").detected.agg(n_cells="size", n_expressing="sum").reset_index()
            d["pct_expressing"] = 100 * d.n_expressing / d.n_cells
            d["eligible"] = d.n_cells >= min_pat_cells
            d["expressing"] = d.eligible & (d.n_expressing >= min_cells) & (d.pct_expressing >= min_pct)
            d["study"] = lab
            donors.append(d)
            e = d[d.eligible]
            row["n_patients_eligible"], row["n_patients_expressing"] = len(e), int(e.expressing.sum())
            if "group" in sub.columns:
                gd = sub.groupby(["group", "donor"]).detected.agg(n_cells="size", n_expressing="sum").reset_index()
                gd["pct"] = 100 * gd.n_expressing / gd.n_cells
                gd["eligible"] = gd.n_cells >= min_pat_cells
                gd["expressing"] = gd.eligible & (gd.n_expressing >= min_cells) & (gd.pct >= min_pct)
                g = gd.groupby("group").agg(n_cells=("n_cells", "sum"), n_expressing=("n_expressing", "sum"),
                                            n_patients_eligible=("eligible", "sum"),
                                            n_patients_expressing=("expressing", "sum")).reset_index()
                g["pct_cells_expressing"] = (100 * g.n_expressing / g.n_cells).round(3)
                g["study"] = lab
                groups.append(g)
        ov.append(row)
        d4 = Path(r["run"]) / "4_DE"
        if d4.exists():
            fr = d4 / "fisher_celltype_vs_rest.csv"
            if fr.exists():
                x = pd.read_csv(fr)
                x = x[(x.scope == "all groups") & x.iloc[:, 1].astype(str).isin(members_all)].copy()
                if len(x):
                    x.insert(0, "study", lab)
                    x = x.rename(columns={x.columns[2]: "cell_type"})
                    f_rest.append(x)
            fb = d4 / "fisher_between_groups.csv"
            if fb.exists():
                x = pd.read_csv(fb)
                x = x[x.iloc[:, 0].astype(str).isin(members_all)].copy()
                if len(x):
                    x.insert(0, "study", lab)
                    x = x.rename(columns={x.columns[1]: "cell_type"})
                    f_btw.append(x)
            for m in members_all:
                p = d4 / safe_name(m) / "pooled" / "nebula_de_results.csv"
                if p.exists():
                    t = pd.read_csv(p, usecols=["gene", "log2FC", "p_value", "fdr", "reliable_fit"])
                    t["study"], t["cell_type"] = lab, m
                    de.append(t)
    ov = pd.DataFrame(ov)
    ov["patients_expressing"] = ov.apply(
        lambda x: f"{int(x.n_patients_expressing)}/{int(x.n_patients_eligible)} "
                  f"({100 * x.n_patients_expressing / x.n_patients_eligible:.1f}%)"
        if x.n_patients_eligible == x.n_patients_eligible and x.n_patients_eligible > 0 else "n/a", axis=1)
    ovp = ov[ov.present].copy()
    pat = ovp[ovp.n_patients_eligible > 0][["study", "tissue", "n_patients_eligible", "n_patients_expressing"]].copy()
    if len(pat):
        pat["pct_patients_expressing"] = (100 * pat.n_patients_expressing / pat.n_patients_eligible).round(1)
        ci = [wilson_ci(a, b) for a, b in zip(pat.n_patients_expressing, pat.n_patients_eligible)]
        pat["ci95_low"], pat["ci95_high"] = [round(c[0], 1) for c in ci], [round(c[1], 1) for c in ci]
    dn = pd.concat(donors, ignore_index=True) if donors else pd.DataFrame()
    if len(dn):
        dn = dn[dn.study.isin(ovp.study) & dn.eligible]
    gr = pd.concat(groups, ignore_index=True) if groups else pd.DataFrame()
    if len(gr):
        gr = gr[gr.study.isin(ovp.study)]
        gr["patients_expressing"] = gr.apply(lambda x: f"{int(x.n_patients_expressing)}/{int(x.n_patients_eligible)}",
                                             axis=1)
    fr = pd.concat(f_rest, ignore_index=True) if f_rest else pd.DataFrame()
    fb = pd.concat(f_btw, ignore_index=True) if f_btw else pd.DataFrame()
    de = pd.concat(de, ignore_index=True) if de else pd.DataFrame()
    cons = pd.DataFrame()
    if len(de):
        de = de[de.reliable_fit.astype(str).str.lower().isin(["true", "1"])]
        de["up"] = (de.fdr < fdr_thr) & (de.log2FC > 0)
        de["down"] = (de.fdr < fdr_thr) & (de.log2FC < 0)
        best = de.sort_values("fdr").drop_duplicates(["gene", "study"])
        cons = best.groupby("gene").agg(n_studies_tested=("study", "nunique"), n_studies_up=("up", "sum"),
                                        n_studies_down=("down", "sum"), median_log2FC=("log2FC", "median"),
                                        min_fdr=("fdr", "min")).reset_index()
        cons["studies_up"] = cons.gene.map(best[best.up].groupby("gene").study.apply(lambda s: ", ".join(sorted(s))))
        cons["studies_down"] = cons.gene.map(best[best.down].groupby("gene").study.apply(lambda s: ", ".join(sorted(s))))
        cons["consistency"] = cons[["n_studies_up", "n_studies_down"]].max(axis=1)
        cons = cons.sort_values(["consistency", "min_fdr"], ascending=[False, True])
        cons["median_log2FC"] = cons.median_log2FC.round(3)
    return {"target": target, "name": name, "kind": kind, "members": members_all, "overview": ov,
            "patients": pat, "donors": dn, "groups": gr, "fisher_rest": fr, "fisher_between": fb,
            "de_long": de, "de_consistency": cons,
            "params": {"min_cells": min_cells, "min_pct": min_pct, "min_ct_cells": min_ct_cells,
                       "min_pat_cells": min_pat_cells, "fdr": fdr_thr}}


def make_celltype_plots(C, gene, out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    key, nm = C["name"], C["name"]
    ov = C["overview"][C["overview"].present].sort_values("pct_expressing")
    if len(ov):
        tissues = ov.tissue.unique().tolist()
        cm = dict(zip(tissues, _palette(len(tissues))))
        fig, ax = plt.subplots(figsize=(8.5, max(2.5, 0.32 * len(ov) + 1)))
        ax.barh(range(len(ov)), ov.pct_expressing, color=[cm[t] for t in ov.tissue])
        for i, x in enumerate(ov.itertuples()):
            ax.text(x.pct_expressing, i, f"  {x.n_expressing:,}/{x.n_cells:,}", va="center", fontsize=7)
        ax.set_yticks(range(len(ov))); ax.set_yticklabels(ov.study, fontsize=8)
        ax.set_xlabel(f"% of {nm} cells expressing {gene}")
        ax.legend(handles=[Line2D([], [], marker="s", ls="", color=cm[t], label=t[:35]) for t in tissues],
                  title="tissue", fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
        ax.set_title(f"{gene} in {nm}: detection per study", fontsize=10)
        for s_ in ["top", "right"]:
            ax.spines[s_].set_visible(False)
        _save_ct(fig, out / "bar_pct_by_study.png", "% of cells of this type expressing the gene, per study", key)
    dn = C["donors"]
    if len(dn):
        order = [s for s in C["overview"].sort_values("pct_expressing").study if s in set(dn.study)]
        fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(order) + 3), 4.6))
        rng = np.random.default_rng(0)
        for i, s in enumerate(order):
            v = dn[dn.study == s]
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v.pct_expressing, s=16, alpha=0.8,
                       c=np.where(v.expressing, "#c0392b", "#999999"), zorder=3)
            ax.hlines(v.pct_expressing.median(), i - 0.28, i + 0.28, color="black", lw=1.5)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=60, ha="right", fontsize=8)
        ax.set_ylabel(f"% of a patient's {nm} cells expressing {gene}")
        ax.set_title(f"Patient-level {gene} in {nm} (one dot = one patient; red = expressing)", fontsize=10)
        for s_ in ["top", "right"]:
            ax.spines[s_].set_visible(False)
        _save_ct(fig, out / "strip_patients.png", "Each patient's % expressing cells of this type, by study", key)
    gr = C["groups"]
    if len(gr) and gr.group.nunique() <= 40:
        m = gr.pivot_table(index="study", columns="group", values="pct_cells_expressing")
        r_, c_ = m.shape
        fig, ax = plt.subplots(figsize=(max(5, 0.5 * c_ + 4), max(2.5, 0.38 * r_ + 1.8)))
        cmap = plt.get_cmap("magma_r").copy(); cmap.set_bad("#eeeeee")
        im = ax.imshow(np.ma.masked_invalid(m.values.astype(float)), aspect="auto", cmap=cmap, vmin=0)
        for i in range(r_):
            for j in range(c_):
                v = m.values[i, j]
                if v == v:
                    ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7,
                            color="white" if v > 0.6 * np.nanmax(m.values) else "black")
        ax.set_xticks(range(c_)); ax.set_xticklabels(m.columns, rotation=50, ha="right", fontsize=8)
        ax.set_yticks(range(r_)); ax.set_yticklabels(m.index, fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.6).set_label("% cells expressing")
        ax.set_title(f"{gene} in {nm}: % expressing per group\n(compare within a row; grey = group absent)", fontsize=10)
        _save_ct(fig, out / "heatmap_groups.png", "Detection per group within this cell type, per study", key)
    cons = C["de_consistency"]
    if len(cons):
        top = cons[cons.consistency > 0].head(25)
        if len(top):
            fig, ax = plt.subplots(figsize=(7.5, max(2.5, 0.3 * len(top) + 1)))
            y = np.arange(len(top))
            ax.barh(y, top.n_studies_up, color="#c0392b", label="higher in gene+ cells (FDR < thr)")
            ax.barh(y, -top.n_studies_down, color="#2471a3", label="higher in gene- cells")
            for yi, x in zip(y, top.itertuples()):
                ax.text(max(x.n_studies_up, 0) + 0.1, yi, f"of {x.n_studies_tested}", va="center", fontsize=7)
            ax.axvline(0, color="#555", lw=0.7)
            ax.set_yticks(y); ax.set_yticklabels(top.gene, fontsize=8); ax.invert_yaxis()
            ax.set_xlabel("number of studies (right: up, left: down)")
            from matplotlib.ticker import MaxNLocator, FuncFormatter
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(int(v))}"))
            ax.legend(frameon=False, fontsize=7, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
            ax.set_title(f"Genes consistently DE in {gene}+ vs {gene}- {nm} cells (NEBULA)", fontsize=10, pad=22)
            for s_ in ["top", "right"]:
                ax.spines[s_].set_visible(False)
            _save_ct(fig, out / "de_consistency.png", "Genes DE in several studies (NEBULA, pooled mode)", key)
            de = C["de_long"]
            best = de.sort_values("fdr").drop_duplicates(["gene", "study"])
            m = best[best.gene.isin(top.gene)].pivot_table(index="gene", columns="study", values="log2FC").reindex(top.gene)
            q = best[best.gene.isin(top.gene)].pivot_table(index="gene", columns="study", values="fdr").reindex(top.gene)
            v = np.nanpercentile(np.abs(m.values), 98) if np.isfinite(m.values).any() else 1
            fig, ax = plt.subplots(figsize=(max(5, 0.6 * m.shape[1] + 4), max(3, 0.3 * m.shape[0] + 1.5)))
            cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#eeeeee")
            im = ax.imshow(np.ma.masked_invalid(m.values.astype(float)), cmap=cmap, vmin=-v, vmax=v, aspect="auto")
            for i in range(m.shape[0]):
                for j in range(m.shape[1]):
                    if q.values[i, j] == q.values[i, j] and q.values[i, j] < C["params"]["fdr"]:
                        ax.text(j, i, "*", ha="center", va="center", fontsize=9)
            ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=50, ha="right", fontsize=8)
            ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontsize=8)
            fig.colorbar(im, ax=ax, shrink=0.5).set_label("log2FC gene+ vs gene-")
            ax.set_title(f"log2FC across studies in {nm} (* FDR < {C['params']['fdr']:g}; grey = not tested)", fontsize=10)
            _save_ct(fig, out / "de_heatmap.png", "Fold changes of the most consistent genes across studies", key)


def export_celltype(C, out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    sheets = {"overview": C["overview"], "patients": C["patients"], "patients_all": C["donors"],
              "groups": C["groups"], "fisher_vs_rest": C["fisher_rest"], "fisher_between_groups": C["fisher_between"],
              "de_consistency": C["de_consistency"], "de_all_studies": C["de_long"]}
    for k, d in sheets.items():
        if isinstance(d, pd.DataFrame) and len(d):
            d.to_csv(out / f"{k}.csv", index=False)
    return sheets


CT_NOTES = {
    "overview": "Per study: cells of this type, cells expressing (≥1 count), % = expressing ÷ cells × 100. "
                "<i>rank_in_study</i> = position of this cell type when the study's cell types (with enough cells) "
                "are sorted by % expressing. For a broad class, counts are summed over its member cell types.",
    "patients": "Only patients with at least {min_pat_cells} cells of this type are counted (eligible). A patient "
                "expresses the gene in this cell type if ≥ {min_cells} of those cells express it and they are "
                "≥ {min_pct}% of them. 95% CI = Wilson interval.",
    "groups": "Within each study: % of this cell type's cells expressing the gene in each group, and eligible "
              "patients expressing / eligible patients per group. Compare groups within a study.",
    "fisher": "Taken from each study's DE step (per-study notebook, section 4). Odds ratio > 1: gene+ cells are "
              "relatively more common in this cell type (or in group 1) than in all other cells (or group 2). "
              "Fisher's test counts cells, so it ignores donor-to-donor variation.",
    "de": "Taken from each study's NEBULA results (pooled mode, reliable fits only). For each gene: number of studies "
          "where it is significantly higher (FDR < {fdr}) in gene+ cells, number where it is lower, and the median "
          "log2 fold change across studies. A gene DE in several independent studies is a more robust finding.",
}


def celltype_report_html(Cs, gene, details, out):
    css = ("body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1200px;margin:auto;"
           "padding:24px;background:#fff;color:#222}table.df{border-collapse:collapse;font-size:12px}"
           "table.df td,table.df th{border:1px solid #ddd;padding:3px 6px}.grid{display:flex;flex-wrap:wrap;gap:14px}"
           "figure{margin:0;flex:1 1 520px}figure img{width:100%}figcaption{font-size:12px;color:#555}"
           ".note{background:#f6f8fa;border-left:4px solid #4c72b0;padding:6px 10px;margin:6px 0;font-size:13px}"
           ".tbl{overflow-x:auto}h2{border-bottom:2px solid #ddd;margin-top:36px}")
    h = [f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{gene}: cell-type report</title><style>{css}"
         f"</style></head><body><h1>{gene}: cell-type report</h1>",
         "<p>Cell types: " + ", ".join(C["name"] for C in Cs) + ". Companion to the all-cells cross-study summary.</p>"]
    for C in Cs:
        p = C["params"]
        note = lambda k: f"<div class='note'>{CT_NOTES[k].format(**p)}</div>"
        figs = lambda f: "<div class='grid'>" + "".join(
            f"<figure><img src='data:image/png;base64,{img_b64(pp)}'><figcaption>{c}</figcaption></figure>"
            for k, c, pp in CT_FIGS if k == C["name"] and pp.name.startswith(f)) + "</div>"
        tbl = lambda d, n=40: f"<div class='tbl'>{d.head(n).to_html(index=False, classes='df', border=0)}</div>" if len(d) else "<p><i>Not available.</i></p>"
        h.append(f"<h2>{C['name']}" + (f" <small>(broad class: {', '.join(C['members'][:8])}"
                                       f"{'…' if len(C['members']) > 8 else ''})</small>" if C["kind"] == "broad" else "")
                 + "</h2>")
        ov = C["overview"]
        h.append(f"<h3>Studies</h3><p>{int(ov.present.sum())} of {len(ov)} studies have at least {p['min_ct_cells']} "
                 f"cells of this type.</p>" + note("overview")
                 + tbl(ov[["study", "tissue", "n_cells", "n_expressing", "pct_expressing", "rank_in_study",
                           "patients_expressing", "present"]]) + figs("bar_pct"))
        h.append("<h3>Patients</h3>" + note("patients") + tbl(C["patients"]) + figs("strip"))
        cols = [c for c in ["study", "group", "n_cells", "n_expressing", "pct_cells_expressing", "patients_expressing"]
                if c in C["groups"].columns]
        h.append("<h3>Groups</h3>" + note("groups") + (tbl(C["groups"][cols], 80) if len(C["groups"]) else "<p><i>No group information.</i></p>")
                 + figs("heatmap_groups"))
        if len(C["fisher_rest"]) or len(C["fisher_between"]):
            h.append("<h3>Cell proportions (Fisher's exact, from each study's DE step)</h3>" + note("fisher"))
            if len(C["fisher_rest"]):
                h.append(tbl(C["fisher_rest"][["study", "cell_type", "pct_pos_celltype", "pct_pos_other_cells",
                                               "odds_ratio", "or_ci95_low", "or_ci95_high", "fdr"]].round(4)))
            h.append(figs("fisher"))
            if len(C["fisher_between"]):
                h.append("<p><b>Between groups</b></p>" + tbl(C["fisher_between"][[
                    "study", "cell_type", "group_1", "group_2", "pct_pos_group_1", "pct_pos_group_2",
                    "odds_ratio_g1_vs_g2", "fdr"]].round(4)))
        else:
            h.append("<h3>Cell proportions and DE</h3><p><i>No DE results for this cell type. Run DE in the per-study "
                     "notebook (section 4) to add Fisher's and NEBULA results here.</i></p>")
        if len(C["de_consistency"]):
            h.append("<h3>Expression: genes DE across studies (NEBULA)</h3>" + note("de")
                     + tbl(C["de_consistency"].head(40)) + figs("de_"))
    h.append("<h2>References</h2><ol>" + "".join(f"<li>{c}</li>" for c in details.citation) + "</ol></body></html>")
    p = Path(out) / f"{gene}_celltype_report.html"
    p.write_text("\n".join(h))
    try:
        with pd.ExcelWriter(Path(out) / f"{gene}_celltype_report.xlsx") as xw:
            for C in Cs:
                for k, d in export_celltype(C, Path(out) / safe_name(C["name"])).items():
                    if isinstance(d, pd.DataFrame) and len(d):
                        d.to_excel(xw, sheet_name=f"{C['name'][:18]}_{k}"[:31], index=False)
    except ImportError:
        for C in Cs:
            export_celltype(C, Path(out) / safe_name(C["name"]))
    return p


# ----------------------------------------------------------------------------- calculations page & exports
CALC_CSS = """
.calcpage{background:#fff;color:#222;max-width:1100px;padding:20px 28px;border:1px solid #e5e5e5;
 font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;line-height:1.45}
.calcpage h3{margin-top:24px;border-bottom:1px solid #eee;padding-bottom:3px}
table.calc{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}
table.calc th{background:#f4f6f8;text-align:left} table.calc td,table.calc th{border:1px solid #ddd;padding:5px 8px;vertical-align:top}
"""


def calc_page(S):
    p, inv, gene = S["params"], S["inventory"], S["gene"]
    x = inv.iloc[0]
    rows = lambda items: ('<table class="calc"><tr><th>Item</th><th>How it is calculated</th><th>Example from your data'
                          '</th></tr>' + "".join(f"<tr><td><b>{a}</b></td><td>{b}</td><td>{c}</td></tr>"
                                                  for a, b, c in items) + "</table>")
    ctl = S["celltype_long"]
    e = ctl[(ctl.study == x.study) & ctl.enough_cells].sort_values("pct_expressing", ascending=False)
    e = e.iloc[0] if len(e) else ctl[ctl.study == x.study].iloc[0]
    b = ctl[(ctl.study == x.study) & (ctl.broad_class == e.broad_class)]
    cons = S["consistency_broad"]
    c0 = cons.iloc[0] if len(cons) else None
    pat = S["patients"]
    p0 = pat[pat.study == x.study].iloc[0] if len(pat) and x.study in set(pat.study) else None
    h = [f'<div class="calcpage"><h2>How the cross-study summary is calculated</h2>'
         f"<p>Built only from each study's own outputs. <b>Cells are never pooled across studies</b> and expression "
         f"levels are <b>not</b> compared between studies, because each study deposits differently normalized "
         f"data on CELLxGENE. The summary uses <b>proportions</b> (detection) and <b>patients</b>, which share "
         f"the same meaning in every study.</p>",
         "<h3>Inventory table</h3>", rows([
             ("pct_cells_expressing", "cells expressing ÷ all cells × 100 (from overall_summary.csv of that study)",
              f"{x.study}: {x.n_cells_expressing:,} ÷ {x.n_cells:,} × 100 = {x.pct_cells_expressing}%"),
             ("top_cell_type", f"cell type with the highest % expressing, among cell types with ≥ {p['min_ct_cells']} "
                               "cells in that study", f"{x.study}: {x.top_cell_type} ({x.top_cell_type_pct}%)"),
             ("patients_expressing", f"Re-computed here with one common rule for all studies: a patient expresses the "
                                     f"gene if ≥ {p['min_cells']} expressing cells AND ≥ {p['min_pct']}% of its cells. "
                                     "A donor with samples in several groups is counted once (cells summed).",
              f"{x.study}: {x.patients_expressing}"),
             ("pct_primary_cells", "share of the study's cells flagged by CELLxGENE as original (not re-used from "
                                   "another dataset)", f"{x.pct_primary_cells}%"),
         ]),
         "<h3>Study × cell type heatmaps</h3>", rows([
             ("Fine cell types", "Each cell: that study's % of cells of that CELLxGENE cell type expressing the gene "
                                 f"= n_expressing ÷ n_cells × 100. Grey if the study has < {p['min_ct_cells']} cells "
                                 "of that type.",
              f"{x.study} / {e.cell_type}: {int(e.n_expressing):,} ÷ {int(e.n_cells):,} × 100 = "
              f"{e.pct_expressing:.2f}%"),
             ("Broad classes", "Cell types are mapped to broad classes (editable CSV). Per study, counts are summed "
                               "over all cell types in the class first, then % = Σ n_expressing ÷ Σ n_cells × 100.",
              f"{x.study} / {e.broad_class}: {int(b.n_expressing.sum()):,} ÷ {int(b.n_cells.sum()):,} × 100 = "
              f"{100 * b.n_expressing.sum() / b.n_cells.sum():.2f}%"),
         ]),
         "<h3>Consistency tables</h3>", rows([
             ("n_studies_with_celltype", f"studies having ≥ {p['min_ct_cells']} cells of this type",
              f"{c0.broad_class}: {c0.n_studies_with_celltype}" if c0 is not None else ""),
             (f"n_studies_≥X%", "of those, studies where ≥ X% of these cells express the gene (X = 1, 5, "
                                f"{p['expr_thr']:g}, 25)",
              f"{c0.broad_class}: {c0.summary}" if c0 is not None else ""),
             ("n_studies_top_celltype", "studies in which this is the highest-% cell type (or class)",
              f"{c0.broad_class}: {c0.n_studies_top_celltype}" if c0 is not None else ""),
             ("median / min / max_pct", "median, minimum and maximum of the per-study % across studies that have it",
              f"{c0.broad_class}: median {c0.median_pct}%, range {c0.min_pct}–{c0.max_pct}%" if c0 is not None else ""),
         ]),
         "<h3>Patient plots</h3>", rows([
             ("% patients expressing", "patients expressing ÷ patients × 100, per study (common rule above)",
              f"{x.study}: {x.patients_expressing}"),
             ("95% CI", "Wilson score interval for a proportion; wide intervals mean few patients",
              f"{x.study}: {p0.ci95_low}–{p0.ci95_high}%" if p0 is not None else ""),
             ("Patient strip plot", "one dot per patient = that patient's expressing cells ÷ its cells × 100; "
                                    "red = meets the rule; black line = median of the study's patients", ""),
         ]),
         "<h3>Study details</h3>", rows([
             ("citation", "Full reference looked up from the study DOI via Crossref (cached in "
                          "<code>citations_cache.json</code>); otherwise built from CELLxGENE's publisher information "
                          "(authors, journal, year). A 'citation_override' in the editable sheet always wins.", ""),
             ("samples used", "From this study's run: the groups you selected, donors (and samples, if the authors "
                              "provide a sample column) actually included, and their cell numbers "
                              "(<code>samples_used.csv</code>).", ""),
             ("how samples were obtained", "CELLxGENE has no standard field for this. Shown: CELLxGENE tissue type "
                                           "(tissue / organoid / cell culture) and the authors' own sample columns "
                                           "(e.g. sample_type, procedure, site), plus the authors' collection "
                                           "description and data/protocol links. Your entry in "
                                           "<code>study_details_editable.csv</code> replaces the automatic text.", ""),
             ("technology", "CELLxGENE assay labels (share of cells) and suspension type (single-cell vs "
                            "single-nucleus).", ""),
         ]),
         "<h3>Groups</h3>", rows([
             ("Study × group heatmap", "within each study, % cells expressing per group (e.g. disease) = Σ expressing "
                                       "÷ Σ cells of that group. Compare groups within a study (row), not across rows.",
              ""),
             ("subjects_expressing", "per group, subjects meeting the rule ÷ subjects in that group", ""),
         ]),
         "</div>"]
    return "\n".join(h)


def img_b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()


def export_all(S, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    sheets = {"inventory": S["inventory"],
              "study_details": S["details"].drop(columns=["author_description"]),
              "samples_used": S["samples_used"], "consistency_broad": S["consistency_broad"],
              "consistency_fine": S["consistency_fine"], "heatmap_broad_pct": S["broad_pct"].round(3),
              "heatmap_fine_pct": S["fine_pct"].round(3), "patients_per_study": S["patients"],
              "patients_all": S["donors"], "groups_per_study": S["groups"], "celltype_long": S["celltype_long"]}
    for name, df in sheets.items():
        if isinstance(df, pd.DataFrame) and not df.empty:
            df.to_csv(out / f"{name}.csv", index=name.startswith("heatmap"))
    try:
        with pd.ExcelWriter(out / f"{S['gene']}_cross_study_summary.xlsx") as xw:
            for name, df in sheets.items():
                if isinstance(df, pd.DataFrame) and not df.empty:
                    df.to_excel(xw, sheet_name=name[:31], index=name.startswith("heatmap"))
    except ImportError:
        pass
    (out / "HOW_IT_IS_CALCULATED.html").write_text(
        f"<!DOCTYPE html><html><head><meta charset='utf-8'><style>body{{background:#fff}}{CALC_CSS}</style></head>"
        f"<body>{calc_page(S)}</body></html>")
    inv = S["inventory"]
    cols = ["ref", "study", "tissue", "assay", "groups", "n_cells", "pct_cells_expressing", "top_cell_type",
            "top_cell_type_pct", "patients_expressing", "pct_primary_cells"]
    html = [f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{S['gene']} cross-study summary</title><style>"
            "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1200px;margin:auto;padding:24px}"
            "table.df{border-collapse:collapse;font-size:12px}table.df td,table.df th{border:1px solid #ddd;padding:3px 6px}"
            ".grid{display:flex;flex-wrap:wrap;gap:14px}figure{margin:0;flex:1 1 520px}figure img{width:100%}"
            "figcaption{font-size:12px;color:#555}.warn{background:#fff6e5;border-left:4px solid #e69500;padding:6px 10px;"
            "margin:4px 0}.tbl{overflow-x:auto}</style></head><body>",
            f"<h1>{S['gene']}: cross-study summary ({len(inv)} studies)</h1>"]
    for w in S["warnings"]:
        html.append(f"<div class='warn'>{w}</div>")
    html.append(f"<style>{DETAILS_CSS}</style>")
    html.append(f"<h2>Inventory</h2><div class='tbl'>{inv[[c for c in cols if c in inv.columns]].to_html(index=False, classes='df', border=0)}</div>")
    if not S["consistency_broad"].empty:
        html.append("<h2>Consistency (broad classes)</h2><div class='tbl'>"
                    + S["consistency_broad"].to_html(index=False, classes="df", border=0) + "</div>")
    for tab in dict.fromkeys(t for t, _, _ in FIGS):
        html.append(f"<h2>{tab}</h2><div class='grid'>" + "".join(
            f"<figure><img src='data:image/png;base64,{img_b64(p)}'><figcaption>{c}</figcaption></figure>"
            for t, c, p in FIGS if t == tab) + "</div>")
    html.append("<h2>Study details</h2>" + details_html(S["details"]))
    html.append("<h2>References</h2><ol>" + "".join(f"<li>{c}</li>" for c in S["details"].citation) + "</ol>")
    html.append("</body></html>")
    (out / f"{S['gene']}_cross_study_summary.html").write_text("\n".join(html))


def make_all_plots(S, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    FIGS.clear()
    g, p = S["gene"], S["params"]
    plot_study_overall(S["inventory"], out / "bar_overall_pct_by_study.png", g, "Overview")
    plot_heatmap(S["broad_pct"], out / "heatmap_study_x_broad_class.png", f"% cells expressing {g}: study × broad cell class",
                 "Study × broad cell class (% cells expressing)", "Cell types")
    plot_heatmap(S["fine_pct"], out / "heatmap_study_x_celltype.png", f"% cells expressing {g}: study × cell type",
                 "Study × CELLxGENE cell type (most widely shared cell types shown; full matrix in CSV)", "Cell types")
    plot_consistency(S["consistency_broad"], "broad_class", p["expr_thr"], out / "consistency_broad_class.png",
                     f"In how many studies is {g} expressed (≥{p['expr_thr']:g}% of cells)? Broad classes",
                     "Consistency across studies: broad cell classes", "Consistency")
    plot_consistency(S["consistency_fine"], "cell_type", p["expr_thr"], out / "consistency_celltype.png",
                     f"In how many studies is {g} expressed (≥{p['expr_thr']:g}% of cells)? Cell types",
                     "Consistency across studies: CELLxGENE cell types (top 40)", "Consistency")
    rule = f"≥{p['min_cells']} cells & ≥{p['min_pct']}%"
    plot_donor_strip(S["donors"], S["order"], out / "strip_patient_pct_by_study.png", g, "Patients")
    plot_groups(S["groups"], out / "heatmap_study_x_group.png", g, "Groups")
