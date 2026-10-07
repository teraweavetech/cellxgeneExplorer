"""Browser app. Start with:  cxg-explorer app"""
import json
import shutil
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(page_title="CELLxGENE Gene Explorer", layout="wide")

from cxg_explorer import agent as A  # noqa: E402
from cxg_explorer import interactive as IX  # noqa: E402
from cxg_explorer import notebook_agent as NA  # noqa: E402
from cxg_explorer import summary as SM  # noqa: E402
from cxg_explorer import theme as T  # noqa: E402
from cxg_explorer import __version__  # noqa: E402

ss = st.session_state
T.inject()
T.reset_panels()


def chart(fig, **kw):
    st.plotly_chart(T.style_fig(fig), use_container_width=True, **kw)


# ----------------------------------------------------------------------------- helpers
def html_block(body, css="", height=900):
    components.html(f"<style>body{{background:#fff;font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif}}"
                    f"{css}</style>{body}", height=height, scrolling=True)


def show_images(items, cols=2):
    items = list(items)
    for i in range(0, len(items), cols):
        cs = st.columns(cols)
        for c, (cap, p) in zip(cs, items[i:i + cols]):
            c.image(str(p), caption=cap, use_container_width=True)


def zip_download(folder, label, key):
    folder = Path(folder)
    if folder.exists():
        tmp = Path(tempfile.mkdtemp()) / folder.name
        arch = shutil.make_archive(str(tmp), "zip", folder)
        st.download_button(label, Path(arch).read_bytes(), file_name=f"{folder.name}.zip", key=key)


def file_download(path, label, key, mime="text/html"):
    path = Path(path)
    if path.exists():
        st.download_button(label, path.read_bytes(), file_name=path.name, mime=mime, key=key)


def section_view(prefix):
    secs = list(dict.fromkeys([s for s, *_ in A.FIGS + A.TABLES] + [s for s, _ in A.NOTES]))
    secs = [s for s in secs if s.startswith(prefix)]
    if not secs:
        st.caption("Not run for this analysis.")
    for s in secs:
        st.subheader(s)
        for sec, txt in A.NOTES:
            if sec == s:
                st.markdown(txt, unsafe_allow_html=True)
        for sec, cap, path, prev in A.TABLES:
            if sec == s:
                st.markdown(f"**{cap}** ({path.name})")
                st.dataframe(pd.read_csv(path), use_container_width=True, height=min(420, 38 + 35 * len(prev)))
        show_images([(cap, p) for sec, cap, p in A.FIGS if sec == s])


# ----------------------------------------------------------------------------- single study
def open_dataset(path, meta, cx):
    D = A.CxgDataset(path)
    NA.STATE["cx"] = cx or {}
    ss.D, ss.meta, ss.cx = D, meta, cx or {}
    ss.res = ss.de = ss.report_body = None
    coll = (cx or {}).get("collection")
    ss.study_name = NA.default_study_name(coll) if coll else NA.default_study_name(fallback=Path(path).stem.split("__")[0])


def page_single(out_dir, data_dir):
    if not ss.get("res"):
        T.page_header("cellxgeneExplorer", subtitle="Find a paper on CELLxGENE, download its data, and analyse one gene "
                                               "at the cellular, cell-type and patient level.")
        st.write("")
    with st.expander("1. Dataset", expanded=ss.get("D") is None):
        t1, t2 = st.tabs(["Search CELLxGENE", "Open a local H5AD file"])
        with t1:
            q = st.text_area("Paper title as listed on CELLxGENE, a DOI, or a CELLxGENE collection link", height=70)
            if st.button("Search CELLxGENE", type="primary"):
                with st.spinner("Searching CELLxGENE Discover"):
                    try:
                        ss.hits = A.find_collections(q.strip(), Path(data_dir))
                    except Exception as e:
                        st.error(f"Search failed: {e}. Check internet access to api.cellxgene.cziscience.com.")
            if ss.get("hits"):
                labels = [f"{c.get('name')}  ({len(c.get('datasets') or [])} datasets; {why}, score {min(s, 1):.2f})"
                          for s, c, why in ss.hits]
                i = st.selectbox("Collection", range(len(labels)), format_func=lambda k: labels[k])
                if ss.hits[0][0] < 0.45:
                    st.warning("No strong match. Try the paper's DOI or the CELLxGENE collection link.")
                if st.button("Show datasets"):
                    c = ss.hits[i][1]
                    if not c.get("datasets") or "assets" not in (c["datasets"][0] or {}):
                        c = A.get_json(f"{A.CURATION_API}/collections/{c['collection_id']}")
                    ss.coll = c
                    ss.datasets = sorted([d for d in c.get("datasets", []) if d], key=lambda d: -(d.get("cell_count") or 0))
            if ss.get("datasets"):
                def dl(d):
                    h5 = next((a for a in d.get("assets", []) if a.get("filetype") == "H5AD"), {})
                    return (f"{d.get('title')} | {d.get('cell_count', 0):,} cells | {A.labels_of(d.get('assay'))} | "
                            f"H5AD {A.human_size(h5.get('filesize'))}")
                j = st.selectbox("Dataset", range(len(ss.datasets)), format_func=lambda k: dl(ss.datasets[k]))
                st.markdown(f"[Open this collection on CELLxGENE]({A.COLLECTION_URL}{ss.coll['collection_id']})")
                if st.button("Download and open", type="primary"):
                    d, c = ss.datasets[j], ss.coll
                    with st.spinner("Downloading from CELLxGENE (large files take several minutes)"):
                        path = A.download_h5ad(d, Path(data_dir))
                        cx = {"collection": c, "dataset": d}
                        Path(str(path) + ".cellxgene.json").write_text(json.dumps(cx, default=str))
                        open_dataset(path, {"Collection": c.get("name"), "DOI": c.get("doi") or "-",
                                            "Collection URL": f"{A.COLLECTION_URL}{c['collection_id']}",
                                            "Dataset": d.get("title"), "Dataset ID": d.get("dataset_id")}, cx)
                    st.rerun()
        with t2:
            p = st.text_input("Path to an H5AD file downloaded from CELLxGENE")
            if st.button("Open file"):
                pth = Path(p.strip())
                if not pth.exists():
                    st.error(f"File not found: {pth}")
                else:
                    side = Path(str(pth) + ".cellxgene.json")
                    with st.spinner("Opening"):
                        open_dataset(pth, {"Source": f"local file {pth.name}"},
                                     json.loads(side.read_text()) if side.exists() else {})
                    st.rerun()

    D = ss.get("D")
    if D is None:
        st.info("Start by searching for a paper or opening a local H5AD file.")
        return
    obs = D.obs
    if not ss.get("res"):
        st.markdown(f"**Dataset:** {D.title}")
        T.metric_cards([{"label": "Cells", "value": f"{D.n_obs:,}"},
                        {"label": "Cell types", "value": obs["cell_type"].nunique() if "cell_type" in obs.columns else "n/a"},
                        {"label": "Donors", "value": obs["donor_id"].nunique() if "donor_id" in obs.columns else "n/a"},
                        {"label": "Raw counts", "value": "yes" if D.has_raw else "no"}])

    with st.expander("2. Analysis settings", expanded=ss.get("res") is None):
        study = st.text_input("Study folder (results go to <results folder>/<study folder>/<GENE>__<dataset>)",
                              value=ss.get("study_name", ""))
        gcols = NA.group_columns(obs)
        gcol = st.selectbox("Group by", gcols) if gcols else None
        groups = None
        if gcol:
            vc = obs[gcol].astype(str).value_counts()
            groups = T.dropdown_checklist("Groups to analyse (all selected by default; keep one or two to focus)",
                                          list(vc.index), key=f"grp_{gcol}_{D.title}",
                                          fmt=lambda g: f"{g} ({vc[g]:,} cells)")
            if len(groups) == 2:
                st.caption("Two groups selected: direct comparison (Mann-Whitney U on patients).")
            elif len(groups) > 2:
                st.caption("Several groups selected: multi-group comparison (Kruskal-Wallis on patients).")
        gene = st.text_input("Gene of interest (official symbol or Ensembl ID)", value=ss.get("gene_q", ""))
        ctc = "cell_type" if "cell_type" in obs.columns else obs.columns[0]
        cvc = obs[ctc].astype(str).value_counts()
        pops = T.dropdown_checklist("Cell populations for the cell-specific level (optional)", list(cvc.index),
                                    key=f"pops_{D.title}", default=[], fmt=lambda c: f"{c} ({cvc[c]:,} cells)")
        a, b = st.columns(2)
        min_cells = a.number_input("Patient counts as expressing with at least this many expressing cells", 1,
                                   10000, 3)
        min_pct = b.number_input("... and at least this % of its cells", 0.0, 100.0, 1.0, 0.5)
        go = st.button("Run analysis", type="primary")
        if go:
            if gcol and not groups:
                st.error("Select at least one group.")
                return
            vid, sym, sugg = D.resolve_gene(gene)
            if vid is None:
                st.error(f"'{gene}' is not in this dataset." + (f" Did you mean: {', '.join(sugg)}?" if sugg else
                                                                " Use the official symbol (e.g. CD274 for PD-L1)."))
                return
            ss.gene_q = gene
            NA.OUT_DIR = Path(out_dir)
            with st.spinner(f"Analysing {sym}"):
                ss.res = NA.run_analysis(D, vid, sym, pops, int(min_cells), float(min_pct), ss.meta,
                                         study=study.strip() or None, group_col=gcol, groups=groups)
            ss.de = ss.report_body = None
            st.rerun()

    res = ss.get("res")
    if not res:
        return
    out, sym = res["folder"], res["symbol"]
    ov = pd.read_csv(out / "1_cellular_level/overall_summary.csv").iloc[0]
    df = res["df"]
    sub = [f"{int(ov.n_cells_total):,} cells"]
    if res.get("donor_col") and res["donor_col"] in df.columns:
        sub.append(f"{df[res['donor_col']].nunique()} donors")
    sub.append(f"{df[res['ct_col']].nunique()} cell types")
    if A.GROUP_COL in df.columns and df[A.GROUP_COL].nunique() > 1:
        sub.append(f"grouped by {A.GROUP_COL}")
    st.write("")
    T.page_header(sym, crumb=f"Single study / {D.title}", subtitle=" · ".join(sub))
    cards = [{"label": "Cells analysed", "value": f"{int(ov.n_cells_total):,}"},
             {"label": f"Cells expressing {sym}", "value": f"{int(ov.n_cells_expressing):,}",
              "pill": f"{ov.pct_expressing:.1f}%"}]
    hp = out / "3_sample_level/subjects_expressing_headline.csv"
    if hp.exists():
        h = pd.read_csv(hp)
        allr = h[h.group == "all subjects"].iloc[0]
        cards.append({"label": "Patients expressing",
                      "value": f"{int(allr.n_subjects_expressing)} / {int(allr.n_subjects_total)}",
                      "pill": f"{allr.pct_expressing:.0f}%"})
        grp = h[h.group != "all subjects"].group.astype(str).tolist()
        cards.append({"label": "Groups", "chips": grp or ["one group"]})
    T.metric_cards(cards)
    st.caption(f"Saved in {out}")

    tabs = st.tabs(["Cellular", "Cell-specific", "Sample", "Differential expression", "How it's calculated",
                    "Methods", "Descriptive report", "Download"])
    with tabs[0]:
        cs = pd.read_csv(out / "1_cellular_level/celltype_summary.csv")
        st.caption("Hover for details, drag to zoom, double-click to reset.")
        a, b = st.columns([3, 2])
        with a:
            with T.panel():
                if "emb_1" in df.columns:
                    chart(IX.umap_expression(df, sym, res["ct_col"], title=f"{sym} on the UMAP"))
                else:
                    st.caption("This dataset has no embedding.")
        with b:
            with T.panel():
                rows = cs.sort_values("pct_expressing", ascending=False)
                T.bars_panel("% of cells expressing, by cell type",
                             [(r[cs.columns[0]], r.pct_expressing, r.n_expressing, r.n_cells)
                              for _, r in rows.iterrows()])
                if len(rows) > 8:
                    st.caption(f"Top 8 of {len(rows)} cell types; all are in the chart and table below.")
        a, b = st.columns(2)
        with a:
            with T.panel():
                if "emb_1" in df.columns:
                    chart(IX.umap_celltype(df, res["ct_col"]))
        with b:
            with T.panel():
                chart(IX.bar_celltype(cs, sym, name_col=cs.columns[0]))
        section_view("1.")
    with tabs[1]:
        section_view("2.")
    with tabs[2]:
        section_view("3.")
    with tabs[3]:
        de_panel(res)
    with tabs[4]:
        html_block(A.calc_guide_html(out, sym), A.CALC_CSS, height=1400)
    with tabs[5]:
        st.markdown(A.methods_markdown())
    with tabs[6]:
        report_panel(res)
    with tabs[7]:
        file_download(out / "report.html", "Data report (HTML)", "dl_rep")
        file_download(out / "HOW_IT_IS_CALCULATED.html", "How it's calculated (HTML)", "dl_calc")
        file_download(out / "METHODS.md", "Methods (Markdown)", "dl_meth", "text/markdown")
        zip_download(out, "All results (zip)", "dl_zip")


def de_panel(res):
    st.markdown("Within each cell type, compares cells with at least one count of the gene (gene+) against cells "
                "without it (gene-): **expression** of all other genes with NEBULA, and **cell proportions** with "
                "Fisher's exact test.")
    df, ctc = res["df"], res["ct_col"]
    a, b, c = st.columns([3, 1, 1])
    rs = a.text_input("Rscript", value=ss.get("rscript", A.find_rscript() or ""))
    ss.rscript = rs
    if b.button("Check R / nebula"):
        ok, msg = A.r_status(rs or None)
        (st.success if ok else st.warning)(msg.replace("<b>", "").replace("</b>", ""))
    if c.button("Install nebula"):
        with st.spinner("Installing nebula in R (several minutes)"):
            ok, log = A.r_install_nebula(rs)
        (st.success if ok else st.error)("nebula installed." if ok else "Installation failed; see the log below.")
        st.code(log[-2000:])
    cnt = df.groupby(ctc).detected.agg(["sum", "size"])
    mc = st.number_input("Minimum gene+ and gene- cells per cell type", 3, 10000, 20)
    elig = [k for k, r in cnt.iterrows() if r["sum"] >= mc and r["size"] - r["sum"] >= mc]
    cts = T.dropdown_checklist("Cell types (those with enough gene+ and gene- cells are ticked)",
                               list(cnt.sort_values("sum", ascending=False).index),
                               key=f"de_cts_{res['folder']}", default=elig,
                               fmt=lambda k: f"{k} ({int(cnt.loc[k, 'sum']):,} + / "
                                             f"{int(cnt.loc[k, 'size'] - cnt.loc[k, 'sum']):,} -)")
    mode = st.selectbox("Mode", ["Pooled and within each group", "Pooled groups only (group as covariate)",
                                 "Within each group only"])
    pooled, within = mode != "Within each group only", mode != "Pooled groups only (group as covariate)"
    y1, y2, y3, y4 = st.columns(4)
    mx = y1.number_input("Max cells per class", 50, 1_000_000, 5000)
    det = y2.number_input("Min detection %", 0.0, 100.0, 10.0)
    fdr = y3.number_input("FDR", 0.0001, 1.0, 0.05)
    cores = y4.number_input("R cores", 1, 128, 1)
    run = st.button("Run DE", type="primary")
    if run:
        scopes = [s for s, v in [("pooled", pooled), ("within", within)] if v]
        gid = res["meta"]["Gene"].split("(")[-1].rstrip(")")
        prog = st.empty()
        r = A.run_de(ss.D, df, res["symbol"], gid, ctc, res["donor_col"], A.GROUP_COL, res["folder"] / "4_DE",
                     celltypes=cts or None, scopes=scopes or ["pooled"], min_cells=int(mc), max_cells=int(mx),
                     min_det=det / 100, fdr_thr=fdr, rscript=rs or None, ncore=int(cores),
                     progress=lambda m: prog.info(m))
        prog.empty()
        nv = ""
        if "summary" in r and len(r["summary"]) and "nebula_version" in r["summary"].columns \
                and r["summary"]["nebula_version"].notna().any():
            nv = str(r["summary"]["nebula_version"].dropna().iloc[0])
        A.RUN_INFO["de"] = {"scopes": scopes, "min_cells": int(mc), "min_donors": 3, "min_cells_per_donor": 3,
                            "max_cells": int(mx), "min_det": det / 100, "fdr": fdr, "nebula_version": nv}
        out = res["folder"]
        (out / "METHODS.md").write_text("# Methods\n\n" + A.methods_markdown())
        A.write_calc_page(out, res["symbol"])
        A.write_report(out, res["meta"])
        res["digest"] = A.build_digest(out)
        ss.de, ss.de_fdr = r, fdr
    r = ss.get("de")
    if not r:
        return
    msg = r.get("nebula_msg")
    if msg:
        hint = NA.nebula_fix_hint(msg).replace("<b>", "").replace("</b>", "").replace("<code>", "`").replace("</code>", "`")
        st.warning(f"Expression (NEBULA) did not run: {msg}. {hint} The cell-proportion results (Fisher) are complete.")
    else:
        st.success(f"Expression (NEBULA) ran for {len(r.get('jobs', []))} cell type and mode combination(s); "
                   "cell proportions (Fisher) ran for all cell types.")
    st.subheader("A. Expression: NEBULA")
    if "summary" in r and len(r["summary"]):
        st.dataframe(r["summary"], use_container_width=True)
    hm = [(cap, p) for s, cap, p in A.FIGS if s == A.DE_SECTION and p.name == "heatmap_top_de_genes.png"]
    show_images(hm, cols=1)
    if r.get("jobs"):
        k = st.selectbox("Cell type and mode", range(len(r["jobs"])),
                         format_func=lambda i: f"{r['jobs'][i]['celltype']} [{r['jobs'][i]['scope']}]")
        j = r["jobs"][k]
        chart(IX.volcano(j["table"], ss.get("de_fdr", 0.05), f"{j['celltype']} [{j['scope']}]"))
        st.dataframe(j["table"], use_container_width=True, height=380)
    st.subheader("B. Cell proportions: Fisher's exact test")
    if len(r["fisher_vs_rest"]):
        st.markdown("**Each cell type vs all other cells**")
        st.dataframe(r["fisher_vs_rest"].round(5), use_container_width=True)
    if len(r["fisher_between"]):
        st.markdown("**Between groups, within each cell type**")
        st.dataframe(r["fisher_between"].round(5), use_container_width=True)
    show_images([(cap, p) for s, cap, p in A.FIGS if s == A.DE_SECTION and p.name.startswith("fisher")], cols=1)


def report_panel(res):
    st.markdown("Review the results first. When you are ready, a local **gpt-oss:20b** model (via Ollama) writes a "
                "descriptive report from the computed summary statistics only.")
    a, b = st.columns(2)
    model = a.text_input("Model", A.LLM_MODEL)
    url = b.text_input("Server", A.LLM_URL)
    if st.button("Generate descriptive report", type="primary"):
        ok, msg = A.check_llm(model, url)
        if not ok:
            st.warning(msg.replace("<code>", "`").replace("</code>", "`"))
        else:
            box = st.empty()
            body, css = A.generate_descriptive_report(res, model, url, progress=lambda m: box.info(f"{model}: {m}"))
            box.empty()
            ss.report_body, ss.report_css = body, css
    if ss.get("report_body"):
        html_block(ss.report_body, ss.report_css.replace("body{", ".x{"), height=1400)
        file_download(res["folder"] / "DESCRIPTIVE_REPORT.html", "Download report (HTML)", "dl_desc")


# ----------------------------------------------------------------------------- cross-study summary
def page_summary(out_dir):
    T.page_header("Cross-study summary", subtitle="Summarise one gene across all studies analysed in the "
                  "single-study section. Proportions and patients are compared; cells are never pooled across studies.")
    st.write("")
    root = Path(st.text_input("Results folder", value=out_dir))
    if st.button("Scan", type="primary"):
        ss.runs = SM.discover_runs(root) if root.exists() else []
        ss.S = None
        if not ss.runs:
            st.warning("No results found in this folder.")
    runs = ss.get("runs")
    if not runs:
        return
    genes = pd.Series([r["gene"].upper() for r in runs]).value_counts()
    gene = st.selectbox("Gene", list(genes.index), format_func=lambda g: f"{g} ({genes[g]} runs)")
    keys, labels = SM.study_keys(runs, gene)
    st.markdown(f"**{len(keys)} studies** for {gene}. Untick to exclude; choose the run where a study was "
                "analysed more than once.")
    chosen = {}
    for k, rs in sorted(keys.items(), key=lambda kv: labels[kv[0]]):
        c1, c2 = st.columns([1, 2])
        inc = c1.checkbox(labels[k], True, key=f"inc_{labels[k]}")
        opts = [("all groups" if not r["subset"] else "subset: " + ", ".join(map(str, r["groups"]))) +
                f" [{r['label_run']}]" for r in rs]
        i = c2.selectbox("run", range(len(rs)), format_func=lambda x: opts[x], key=f"run_{labels[k]}",
                         label_visibility="collapsed", disabled=len(rs) == 1)
        if inc:
            chosen[k] = rs[i]
    with st.expander("Settings"):
        a, b, c, d = st.columns(4)
        mcell = a.number_input("Patient: min expressing cells", 1, 1000, 3)
        mpct = b.number_input("Patient: min %", 0.0, 100.0, 1.0, 0.5)
        mct = c.number_input("Min cells per cell type", 1, 10000, 20)
        thr = d.number_input("'Expressed' if at least %", 0.1, 100.0, 10.0)
        cite = st.checkbox("Look up full citations online (Crossref)", True)
    if st.button("Build summary", type="primary"):
        if not chosen:
            st.error("Select at least one study.")
            return
        out = root / "_cross_study" / SM.safe_name(gene)
        with st.spinner(f"Summarising {len(chosen)} studies"):
            S = SM.build_summary(chosen, labels, gene, int(mcell), float(mpct), int(mct), float(thr),
                                 mapping_path=root / "_cross_study" / "celltype_to_broad_class.csv",
                                 online_citations=cite)
            SM.make_all_plots(S, out)
            SM.export_all(S, out)
        ss.S, ss.S_out, ss.C = S, out, None
    S = ss.get("S")
    if S is None:
        return
    out, g, inv = ss.S_out, S["gene"], S["inventory"]
    for w in S["warnings"]:
        st.warning(w.replace("<b>", "").replace("</b>", "").replace("<code>", "`").replace("</code>", "`"))
    figs = lambda t: show_images([(c, p) for tt, c, p in SM.FIGS if tt == t])
    tabs = st.tabs(["Overview", "Study details", "Cell types", "Consistency", "Patients", "Groups",
                    "How it's calculated", "Cell-type report", "Download"])
    with tabs[0]:
        st.dataframe(inv, use_container_width=True)
        chart(IX.bar_by_study(inv, "pct_cells_expressing", g, f"Overall detection of {g} per study"))
    with tabs[1]:
        st.caption(f"Edit {S['details_editable']} to add how samples were obtained, notes or citation corrections, "
                   "then build again.")
        html_block(SM.details_html(S["details"]) + "<h3>References</h3><ol>" +
                   "".join(f"<li>{c}</li>" for c in S["details"].citation) + "</ol>", SM.DETAILS_CSS, 1100)
        if len(S["samples_used"]):
            st.markdown("**Samples used (all studies)**")
            st.dataframe(S["samples_used"], use_container_width=True)
    with tabs[2]:
        chart(IX.heatmap(S["broad_pct"], f"% cells expressing {g}: study x broad class"))
        chart(IX.heatmap(S["fine_pct"], f"% cells expressing {g}: study x cell type"))
        st.caption(f"Broad classes come from {S['params']['mapping']} (editable).")
    with tabs[3]:
        st.dataframe(S["consistency_broad"], use_container_width=True)
        st.dataframe(S["consistency_fine"], use_container_width=True)
        figs("Consistency")
    with tabs[4]:
        st.dataframe(S["patients"], use_container_width=True)
        if len(S["donors"]):
            chart(IX.strip_patients(S["donors"], g))
    with tabs[5]:
        if len(S["groups"]):
            st.dataframe(S["groups"], use_container_width=True)
            figs("Groups")
        else:
            st.caption("No group information.")
    with tabs[6]:
        html_block(SM.calc_page(S), SM.CALC_CSS, 1300)
    with tabs[7]:
        celltype_panel(S, out)
    with tabs[8]:
        file_download(out / f"{g}_cross_study_summary.html", "All-cells report (HTML)", "dl_s_html")
        file_download(out / f"{g}_cross_study_summary.xlsx", "All-cells tables (Excel)", "dl_s_xlsx",
                      "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        zip_download(out, "Everything (zip)", "dl_s_zip")


def celltype_panel(S, out):
    opts = SM.celltype_choices(S)
    labels = {v: lab for lab, v in opts}
    pick = T.dropdown_checklist("Cell types of interest (broad classes are marked [class])",
                                [v for _, v in opts], key=f"ctpick_{S['gene']}", default=[],
                                fmt=lambda v: labels[v])
    a, b = st.columns(2)
    mpat = a.number_input("Patient: min cells of the type", 1, 10000, 10)
    fdr = b.number_input("DE FDR", 0.0001, 1.0, 0.05)
    if st.button("Build cell-type report", type="primary"):
        if not pick:
            st.error("Select at least one cell type.")
            return
        p = S["params"]
        SM.CT_FIGS.clear()
        rep = out / "celltype_report"
        Cs = []
        with st.spinner("Building the cell-type report"):
            for v in pick:
                C = SM.build_celltype_summary(S, tuple(v), p["min_cells"], p["min_pct"], p["min_ct_cells"],
                                              int(mpat), float(fdr))
                SM.make_celltype_plots(C, S["gene"], rep / SM.safe_name(C["name"]))
                Cs.append(C)
            ss.C, ss.C_path = Cs, SM.celltype_report_html(Cs, S["gene"], S["details"], rep)
    Cs = ss.get("C")
    if not Cs:
        return
    file_download(ss.C_path, "Cell-type report (HTML)", "dl_ct")
    for C, t in zip(Cs, st.tabs([C["name"][:30] for C in Cs])):
        with t:
            p = C["params"]
            imgs = lambda f: show_images([(c, pp) for k, c, pp in SM.CT_FIGS if k == C["name"] and pp.name.startswith(f)])
            if C["kind"] == "broad":
                st.caption("Broad class made of: " + ", ".join(C["members"]))
            st.markdown("**Studies**")
            st.caption(SM.CT_NOTES["overview"].replace("<i>", "").replace("</i>", ""))
            st.dataframe(C["overview"], use_container_width=True)
            imgs("bar_pct")
            st.markdown("**Patients**")
            st.caption(SM.CT_NOTES["patients"].format(**p))
            st.dataframe(C["patients"], use_container_width=True)
            if len(C["donors"]):
                chart(IX.strip_patients(C["donors"], S["gene"]))
            st.markdown("**Groups**")
            if len(C["groups"]):
                st.dataframe(C["groups"], use_container_width=True)
            imgs("heatmap_groups")
            if len(C["fisher_rest"]) or len(C["fisher_between"]):
                st.markdown("**Cell proportions (Fisher's exact)**")
                st.caption(SM.CT_NOTES["fisher"])
                st.dataframe(C["fisher_rest"], use_container_width=True)
                if len(C["fisher_between"]):
                    st.dataframe(C["fisher_between"], use_container_width=True)
            else:
                st.caption("No DE results for this cell type. Run DE in the single-study section to add them.")
            if len(C["de_consistency"]):
                st.markdown("**Expression: genes DE across studies (NEBULA)**")
                st.caption(SM.CT_NOTES["de"].format(**p))
                st.dataframe(C["de_consistency"].head(60), use_container_width=True)
                imgs("de_")


# ----------------------------------------------------------------------------- help
def page_help():
    T.page_header("Help")
    st.write("")
    st.markdown("""
**Single study.** Search CELLxGENE by paper title, DOI or collection link (or open a local H5AD file), choose the
groups, gene and cell populations, and run the analysis. Results are organised as cellular, cell-specific and
patient levels, with optional differential expression (NEBULA in R, Fisher's exact test) and an optional
descriptive report written by a local gpt-oss:20b model.

**Cross-study summary.** After analysing several studies for the same gene, scan the results folder and build a
summary across studies, then a report for the cell types you choose.

**Every number is explained** in the "How it's calculated" tabs, with worked examples from your data.

**Requirements for optional parts**
- NEBULA differential expression: R 4.4 or newer and the nebula package (`cxg-explorer setup-r`).
- Descriptive report: [Ollama](https://ollama.com) with `ollama pull gpt-oss:20b`.

**Try it without downloading anything:** `cxg-explorer example` writes a small synthetic dataset; open it under
"Open a local H5AD file" and analyse the gene `EXAMPLE1`.
""")


# ----------------------------------------------------------------------------- layout
T.sidebar_brand()
page = st.sidebar.radio("Section", ["Single study", "Cross-study summary", "Help"], label_visibility="collapsed")
st.sidebar.write("")
out_dir = st.sidebar.text_input("Results folder", value=ss.get("out_dir", "cellxgene_results"))
data_dir = st.sidebar.text_input("Download folder", value=ss.get("data_dir", "cellxgene_data"))
ss.out_dir, ss.data_dir = out_dir, data_dir
Path(data_dir).mkdir(parents=True, exist_ok=True)
st.sidebar.caption("Data: CELLxGENE Discover only. Annotations are used as curated by CELLxGENE.")
T.sidebar_credit(__version__)

if page == "Single study":
    page_single(out_dir, data_dir)
elif page == "Cross-study summary":
    page_summary(out_dir)
else:
    page_help()
