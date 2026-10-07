# cellxgeneExplorer
# Copyright (c) 2026 TeraWeave Technologies. Licensed under the PolyForm Noncommercial License 1.0.0 (see LICENSE.md).
# Required Notice: Copyright (c) 2026 TeraWeave Technologies
"""Interactive Jupyter panel for the cross-study summary. Use: `from cxg_explorer import summary_panel; summary_panel()`."""
from .summary import *  # noqa: F401,F403
from . import interactive as IX

import ipywidgets as W
from IPython.display import display, HTML, clear_output

STATE = {}
L = lambda w="95%": W.Layout(width=w)
root_in = W.Text(value="cellxgene_results", description="Results folder", layout=L("70%"),
                 style={"description_width": "110px"})
scan_btn = W.Button(description="Scan", button_style="primary", layout=L("140px"))
gene_dd = W.Dropdown(description="Gene", disabled=True, layout=L("320px"), style={"description_width": "110px"})
studies_box = W.VBox()
min_cells_in = W.BoundedIntText(value=3, min=1, max=1000, description="patient: min cells", layout=L("230px"),
                                style={"description_width": "130px"})
min_pct_in = W.BoundedFloatText(value=1.0, min=0, max=100, step=0.5, description="patient: min %", layout=L("230px"),
                                style={"description_width": "130px"})
min_ct_in = W.BoundedIntText(value=20, min=1, max=10000, description="min cells / cell type", layout=L("230px"),
                             style={"description_width": "130px"})
thr_in = W.BoundedFloatText(value=10, min=0.1, max=100, step=1, description="'expressed' if ≥ %", layout=L("230px"),
                            style={"description_width": "130px"})
cite_cb = W.Checkbox(value=True, indent=False,
                     description="Look up full citations online (Crossref, from each study's DOI)")
ct_pick = W.SelectMultiple(description="Cell types", rows=10, layout=L(), disabled=True,
                           style={"description_width": "110px"})
min_pat_in = W.BoundedIntText(value=10, min=1, max=10000, description="patient: min cells of type", layout=L("260px"),
                              style={"description_width": "160px"})
ct_fdr_in = W.BoundedFloatText(value=0.05, min=0.0001, max=1, step=0.01, description="DE FDR", layout=L("180px"),
                               style={"description_width": "60px"})
ct_btn = W.Button(description="Build cell-type report", button_style="info", disabled=True,
                  layout=L("230px"))
build_btn = W.Button(description="Build summary", button_style="success", disabled=True, layout=L("200px"))
status = W.HTML()
results = W.Output()


def set_status(msg, kind="info"):
    color = {"info": "#1f5f9e", "ok": "#2e7d32", "warn": "#b26a00", "err": "#c62828"}[kind]
    status.value = f'<div style="padding:6px 10px;border-left:4px solid {color};background:#f7f9fb">{msg}</div>'


def guarded(fn):
    def w(_=None):
        try:
            fn()
        except Exception as e:
            import traceback
            set_status(f"{type(e).__name__}: {e}", "err")
            traceback.print_exc()
    return w


@guarded
def on_scan():
    root = Path(root_in.value.strip())
    if not root.exists():
        return set_status(f"Folder not found: {root}", "warn")
    runs = discover_runs(root)
    if not runs:
        return set_status("No agent results found in this folder.", "warn")
    STATE["runs"] = runs
    genes = pd.Series([r["gene"].upper() for r in runs]).value_counts()
    gene_dd.options = [(f"{g}  ({n} runs)", g) for g, n in genes.items()]
    gene_dd.disabled = False
    gene_dd.value = genes.index[0]
    on_gene()
    set_status(f"Found {len(runs)} run(s) for {len(genes)} gene(s). Check the study list, then press "
               "<b>Build summary</b>.", "ok")


def on_gene(_=None):
    if not STATE.get("runs") or not gene_dd.value:
        return
    keys, labels = study_keys(STATE["runs"], gene_dd.value)
    STATE["keys"], STATE["labels"] = keys, labels
    rows, STATE["choices"] = [], {}
    for k, rs in sorted(keys.items(), key=lambda kv: labels[kv[0]]):
        cb = W.Checkbox(value=True, indent=False, layout=L("30px"))
        opts = [((("all groups" if not r["subset"] else "subset: " + ", ".join(map(str, r["groups"])))
                  + f"  [{r['label_run']}]")[:110], i) for i, r in enumerate(rs)]
        dd = W.Dropdown(options=opts, value=0, layout=L("560px"), disabled=len(rs) == 1)
        STATE["choices"][k] = (cb, dd)
        rows.append(W.HBox([cb, W.HTML(f"<div style='width:260px;padding-top:4px'><b>{labels[k]}</b></div>"), dd]))
    studies_box.children = [W.HTML(f"<b>{len(keys)} studies</b> for {gene_dd.value} (untick to exclude):")] + rows
    build_btn.disabled = False


gene_dd.observe(on_gene, "value")


def _grid(items, width=520):
    return HTML('<div style="display:flex;flex-wrap:wrap;gap:14px">' + "".join(
        f'<figure style="margin:0;width:{width}px"><img src="data:image/png;base64,{img_b64(p)}" style="width:100%">'
        f'<figcaption style="font-size:11px;color:#666">{c}</figcaption></figure>' for c, p in items) + "</div>")


@guarded
def on_build():
    selected = {k: STATE["keys"][k][dd.value] for k, (cb, dd) in STATE["choices"].items() if cb.value}
    if len(selected) < 1:
        return set_status("Select at least one study.", "warn")
    gene = gene_dd.value
    root = Path(root_in.value.strip())
    out = root / "_cross_study" / safe_name(gene)
    set_status(f"Summarising {len(selected)} studies ...")
    S = build_summary(selected, STATE["labels"], gene, min_cells_in.value, min_pct_in.value, min_ct_in.value,
                      thr_in.value, mapping_path=root / "_cross_study" / "celltype_to_broad_class.csv",
                      online_citations=cite_cb.value)
    make_all_plots(S, out)
    export_all(S, out)
    STATE["S"], STATE["out"] = S, out
    ct_pick.options = celltype_choices(S)
    ct_pick.value = ()
    ct_pick.disabled = ct_btn.disabled = False
    show(S, out)
    n = len(selected)
    set_status(f"Summary of {n} {'study' if n == 1 else 'studies'} saved in <code>{out}</code>", "ok")


def show(S, out):
    inv, g = S["inventory"], S["gene"]
    tabs, names = [], []

    def tab(name, fn):
        o = W.Output()
        with o:
            fn()
        tabs.append(o); names.append(name)

    def figs(t):
        it = [(c, p) for tt, c, p in FIGS if tt == t]
        if it:
            display(_grid(it))

    def overview():
        for w in S["warnings"]:
            display(HTML(f"<div style='background:#fff6e5;border-left:4px solid #e69500;padding:6px 10px;margin:4px 0'>"
                         f"{w}</div>"))
        cols = ["ref", "study", "tissue", "assay", "organism", "groups", "n_cells", "n_cells_expressing",
                "pct_cells_expressing", "top_cell_type", "top_cell_type_pct", "patients_expressing", "pct_primary_cells"]
        display(HTML(f"<h3>{len(inv)} studies</h3>"))
        display(inv[[c for c in cols if c in inv.columns]])
        show_ix(IX.bar_by_study(inv, "pct_cells_expressing", g, f"Overall detection of {g} per study"), "overview_by_study")
        figs("Overview")

    def celltypes():
        display(HTML(f"<p>Broad classes come from <code>{S['params']['mapping']}</code>; edit that file and press "
                     "Build again to change assignments.</p>"))
        display(HTML("<b>% cells expressing: study × broad class</b>"))
        display(S["broad_pct"].round(2))
        show_ix(IX.heatmap(S["broad_pct"], f"% cells expressing {g}: study x broad class"), "heatmap_broad_class")
        show_ix(IX.heatmap(S["fine_pct"], f"% cells expressing {g}: study x cell type"), "heatmap_cell_type")
        display(HTML("<p><small>Static versions (used in the reports):</small></p>"))
        figs("Cell types")

    def consistency():
        display(HTML("<b>Broad classes</b>"))
        display(S["consistency_broad"])
        display(HTML("<b>CELLxGENE cell types</b> (top 40; all in CSV)"))
        display(S["consistency_fine"].head(40))
        figs("Consistency")

    def patients():
        p = S["params"]
        display(HTML(f"<p>Common rule for every study: a patient expresses {g} with ≥ {p['min_cells']} expressing "
                     f"cells and ≥ {p['min_pct']}% of its cells.</p>"))
        if not S["patients"].empty:
            display(S["patients"])
        if len(S["donors"]):
            show_ix(IX.strip_patients(S["donors"], g), "patients_by_study")
        figs("Patients")

    def groups():
        if S["groups"].empty:
            display(HTML("<i>No group information found.</i>"))
            return
        cols = [c for c in ["study", "group", "n_cells", "n_expressing", "pct_cells_expressing", "subjects_expressing"]
                if c in S["groups"].columns]
        display(S["groups"][cols])
        figs("Groups")

    def details():
        display(HTML(f"<style>{DETAILS_CSS}</style><p>Citation, samples used in this analysis, how the samples were "
                     "obtained, and technology for each study. To fill in or correct the sampling method, notes or a "
                     f"citation, edit <code>{S['details_editable']}</code> and press <b>Build summary</b> again; your "
                     "entries replace the automatic text.</p>" + details_html(S["details"])
                     + "<h3>References</h3><ol>" + "".join(f"<li>{c}</li>" for c in S["details"].citation) + "</ol>"))
        if not S["samples_used"].empty:
            display(HTML("<h3>Samples used (all studies)</h3><p><small>One row per donor (× sample, if the authors "
                         "provide a sample column). Full table in <code>samples_used.csv</code> and the Excel file."
                         "</small></p>"))
            display(S["samples_used"].head(60))

    def show_ix(fig, name):
        IX.save(fig, out / "interactive" / f"{name}.html")
        display(fig)

    def calc():
        display(HTML(f"<style>{CALC_CSS}</style>" + calc_page(S)))

    for n, f in [("Overview", overview), ("Study details", details), ("Cell types", celltypes), ("Consistency", consistency),
                 ("Patients", patients), ("Groups", groups), ("How it's calculated", calc)]:
        tab(n, f)
    STATE["ct_out"] = W.Output()
    with STATE["ct_out"]:
        display(HTML("<p style='padding:12px'>Choose one or more cell types in section 3 of the panel and press "
                     "<b>Build cell-type report</b>.</p>"))
    tabs.append(STATE["ct_out"]); names.append("Cell-type report")
    t = W.Tab(children=tabs)
    for i, n in enumerate(names):
        t.set_title(i, n)
    STATE["tab"] = t
    with results:
        clear_output()
        display(HTML(f"<h2>{g} across {len(inv)} studies</h2><p>Files: <code>{out}</code> "
                     f"(<code>{g}_cross_study_summary.xlsx</code>, <code>{g}_cross_study_summary.html</code>, CSVs, PNGs)</p>"), t)


@guarded
def on_ct():
    S, out = STATE.get("S"), STATE.get("out")
    if S is None:
        return set_status("Build the summary first.", "warn")
    if not ct_pick.value:
        return set_status("Select at least one cell type.", "warn")
    set_status(f"Building the cell-type report for {len(ct_pick.value)} cell type(s) ...")
    CT_FIGS.clear()
    rep_dir = out / "celltype_report"
    Cs = []
    for v in ct_pick.value:
        C = build_celltype_summary(S, tuple(v), min_cells_in.value, min_pct_in.value, min_ct_in.value,
                                   min_pat_in.value, ct_fdr_in.value)
        make_celltype_plots(C, S["gene"], rep_dir / safe_name(C["name"]))
        Cs.append(C)
    rp = celltype_report_html(Cs, S["gene"], S["details"], rep_dir)
    sub, names = [], []
    for C in Cs:
        o = W.Output()
        p = C["params"]
        with o:
            ov = C["overview"]
            figs = lambda f: _grid([(c, pp) for k, c, pp in CT_FIGS if k == C["name"] and pp.name.startswith(f)])
            if C["kind"] == "broad":
                display(HTML(f"<p><b>Broad class</b> made of: {', '.join(C['members'])}</p>"))
            display(HTML(f"<h3>Studies</h3><p>{int(ov.present.sum())} of {len(ov)} studies have ≥ {p['min_ct_cells']} "
                         f"cells of this type.</p><div style='{NOTE_CSS}'>{CT_NOTES['overview']}</div>"))
            display(ov[["study", "tissue", "n_cells", "n_expressing", "pct_expressing", "rank_in_study",
                        "patients_expressing", "present"]])
            display(figs("bar_pct"))
            display(HTML(f"<h3>Patients</h3><div style='{NOTE_CSS}'>{CT_NOTES['patients'].format(**p)}</div>"))
            if len(C["patients"]):
                display(C["patients"])
            display(figs("strip"))
            display(HTML(f"<h3>Groups</h3><div style='{NOTE_CSS}'>{CT_NOTES['groups']}</div>"))
            if len(C["groups"]):
                display(C["groups"][[c for c in ["study", "group", "n_cells", "n_expressing", "pct_cells_expressing",
                                                 "patients_expressing"] if c in C["groups"].columns]])
            display(figs("heatmap_groups"))
            if len(C["fisher_rest"]) or len(C["fisher_between"]):
                display(HTML(f"<h3>Cell proportions (Fisher's exact)</h3><div style='{NOTE_CSS}'>{CT_NOTES['fisher']}</div>"))
                if len(C["fisher_rest"]):
                    display(C["fisher_rest"][["study", "cell_type", "pct_pos_celltype", "pct_pos_other_cells",
                                              "odds_ratio", "or_ci95_low", "or_ci95_high", "fdr"]].round(4))
                display(figs("fisher"))
                if len(C["fisher_between"]):
                    display(HTML("<b>Between groups</b>"))
                    display(C["fisher_between"][["study", "cell_type", "group_1", "group_2", "pct_pos_group_1",
                                                 "pct_pos_group_2", "odds_ratio_g1_vs_g2", "fdr"]].round(4))
            else:
                display(HTML("<h3>Cell proportions and DE</h3><p><i>No DE results for this cell type. Run DE in the "
                             "per-study notebook (section 4) to add Fisher's and NEBULA results here.</i></p>"))
            if len(C["de_consistency"]):
                display(HTML(f"<h3>Expression: genes DE across studies (NEBULA)</h3>"
                             f"<div style='{NOTE_CSS}'>{CT_NOTES['de'].format(**p)}</div>"))
                display(C["de_consistency"].head(40))
                display(figs("de_"))
        sub.append(o); names.append(C["name"][:30])
    st = W.Tab(children=sub)
    for i, n in enumerate(names):
        st.set_title(i, n)
    with STATE["ct_out"]:
        clear_output()
        display(HTML(f"<p>Cell-type report: <code>{rp}</code> (with Excel and CSVs in <code>{rep_dir}</code>). "
                     f"All-cells report: <code>{out / (S['gene'] + '_cross_study_summary.html')}</code>.</p>"), st)
    STATE["tab"].selected_index = len(STATE["tab"].children) - 1
    set_status(f"Cell-type report written: <code>{rp}</code>", "ok")


NOTE_CSS = "background:#f6f8fa;border-left:4px solid #4c72b0;padding:6px 10px;margin:6px 0;font-size:13px"

scan_btn.on_click(on_scan)
ct_btn.on_click(on_ct)
build_btn.on_click(on_build)
PANEL = W.VBox([
    W.HTML("<h3>1. Results</h3>"), W.HBox([root_in, scan_btn]), gene_dd, studies_box,
    W.HTML("<h3>2. Settings</h3>"),
    W.HBox([min_cells_in, min_pct_in]), W.HBox([min_ct_in, thr_in]),
    W.HTML("<small><b>patient rule</b>: applied identically to every study. <b>min cells / cell type</b>: a study's "
           "% for a cell type is shown only if it has at least this many cells of that type. "
           "<b>'expressed' if ≥ %</b>: threshold used in the consistency counts.</small>"),
    cite_cb, build_btn,
    W.HTML("<h3>3. Cell-type report</h3><small>After building the summary, the cell types found across your studies "
           "are listed here (broad classes first, then CELLxGENE cell types, with the number of studies that contain "
           "each). Select the ones you are interested in (Ctrl/Cmd-click for several).</small>"),
    ct_pick, W.HBox([min_pat_in, ct_fdr_in]), ct_btn, status])


def launch(results_folder="cellxgene_results"):
    """Show the cross-study panel for the explorer's results folder."""
    root_in.value = str(results_folder)
    display(PANEL, results)
    set_status("Check the results folder and press <b>Scan</b>.")
