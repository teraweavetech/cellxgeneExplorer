"""Interactive (plotly) figures used by the notebooks and the browser app.

Every function returns a plotly Figure: hover for details, drag to zoom, double-click to reset.
Static PNG versions of the same plots are written by the analysis engines for reports.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

MAX_POINTS = 150_000          # embedding plots are randomly subsampled above this, for responsiveness
TEMPLATE = "simple_white"


def _sample(df, n=MAX_POINTS, seed=0):
    return (df.sample(n, random_state=seed), True) if len(df) > n else (df, False)


def _note(fig, sampled, n_total):
    if sampled:
        fig.add_annotation(text=f"showing a random {MAX_POINTS:,} of {n_total:,} cells", xref="paper", yref="paper",
                           x=0, y=-0.08, showarrow=False, font=dict(size=10, color="#666"))
    return fig


def _hover_cols(df, cols):
    return [c for c in cols if c in df.columns]


def umap_celltype(df, ct_col, title="Cell types"):
    d, sampled = _sample(df)
    fig = go.Figure()
    for i, (ct, x) in enumerate(sorted(d.groupby(ct_col), key=lambda kv: -len(kv[1]))):
        extra = _hover_cols(x, ["donor_id", "disease"])
        fig.add_trace(go.Scattergl(
            x=x["emb_1"], y=x["emb_2"], mode="markers", name=str(ct)[:40],
            marker=dict(size=3, color=px.colors.qualitative.Alphabet[i % 26] if i >= 10 else px.colors.qualitative.Plotly[i]),
            customdata=x[extra].values if extra else None,
            hovertemplate=f"<b>{ct}</b><br>" + "".join(f"{c}: %{{customdata[{k}]}}<br>" for k, c in enumerate(extra))
                          + "<extra></extra>"))
    fig.update_layout(template=TEMPLATE, title=title, height=620, legend=dict(itemsizing="constant", font=dict(size=10)),
                      xaxis=dict(title="UMAP 1", showticklabels=False), yaxis=dict(title="UMAP 2", showticklabels=False))
    return _note(fig, sampled, len(df))


def umap_expression(df, gene, ct_col, title=None):
    d, sampled = _sample(df)
    neg, pos = d[~d["detected"]], d[d["detected"]].sort_values("expression")
    extra = _hover_cols(d, [ct_col, "donor_id", "disease"])
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=neg["emb_1"], y=neg["emb_2"], mode="markers", name=f"{gene} not detected",
                               marker=dict(size=3, color="#d3d3d3"), customdata=neg[extra].values,
                               hovertemplate="".join(f"{c}: %{{customdata[{k}]}}<br>" for k, c in enumerate(extra))
                                             + "not detected<extra></extra>"))
    vmax = float(np.percentile(pos["expression"], 99)) if len(pos) else 1.0
    fig.add_trace(go.Scattergl(x=pos["emb_1"], y=pos["emb_2"], mode="markers", name=f"{gene} expressing",
                               marker=dict(size=4, color=pos["expression"], colorscale="Viridis", cmin=0,
                                           cmax=max(vmax, 1e-6), colorbar=dict(title="expression")),
                               customdata=np.column_stack([pos[extra].values, pos["expression"].round(3)]),
                               hovertemplate="".join(f"{c}: %{{customdata[{k}]}}<br>" for k, c in enumerate(extra))
                                             + f"expression: %{{customdata[{len(extra)}]}}<extra></extra>"))
    fig.update_layout(template=TEMPLATE, title=title or f"{gene} expression", height=620,
                      legend=dict(orientation="h", y=1.04),
                      xaxis=dict(title="UMAP 1", showticklabels=False), yaxis=dict(title="UMAP 2", showticklabels=False))
    return _note(fig, sampled, len(df))


def bar_celltype(summary, gene, name_col=None, title=None):
    s = summary.reset_index() if name_col is None else summary.copy()
    name_col = name_col or s.columns[0]
    s = s.sort_values("pct_expressing")
    fig = go.Figure(go.Bar(
        x=s["pct_expressing"], y=s[name_col].astype(str), orientation="h", marker_color="#4c72b0",
        customdata=np.column_stack([s["n_expressing"], s["n_cells"], s.get("mean_expr_expressing_cells", s["n_cells"] * 0)]),
        hovertemplate="<b>%{y}</b><br>% expressing: %{x:.2f}<br>expressing / cells: %{customdata[0]:,} / "
                      "%{customdata[1]:,}<br>mean in expressing cells: %{customdata[2]:.3f}<extra></extra>"))
    fig.update_layout(template=TEMPLATE, title=title or f"% of cells expressing {gene}",
                      height=max(320, 26 * len(s) + 120), xaxis_title=f"% cells expressing {gene}",
                      margin=dict(l=10))
    return fig


def volcano(t, fdr_thr=0.05, title="NEBULA"):
    d = t.dropna(subset=["fdr"]).copy()
    d["neglog10p"] = -np.log10(d["p_value"].clip(lower=1e-300))
    d["class"] = np.where(d["fdr"] < fdr_thr, np.where(d["log2FC"] > 0, "higher in gene+ cells",
                                                        "higher in gene- cells"), "not significant")
    fig = px.scatter(d, x="log2FC", y="neglog10p", color="class", template=TEMPLATE, title=title,
                     color_discrete_map={"higher in gene+ cells": "#c0392b", "higher in gene- cells": "#2471a3",
                                         "not significant": "#bbbbbb"},
                     hover_data={"gene": True, "log2FC": ":.3f", "fdr": ":.2e", "pct_detected_pos": True,
                                 "pct_detected_neg": True, "neglog10p": False, "class": False}, height=560)
    fig.add_vline(x=0, line_color="#999", line_width=1)
    fig.update_traces(marker=dict(size=6))
    fig.update_layout(xaxis_title="log2 fold change (gene+ vs gene- cells)", yaxis_title="-log10 p")
    return fig


def heatmap(mat, title, colorbar="% cells expressing", diverging=False, max_cols=80):
    m = mat.copy()
    if m.shape[1] > max_cols:
        m = m[m.notna().sum().sort_values(ascending=False).index[:max_cols]]
    z = m.values.astype(float)
    vmax = np.nanmax(np.abs(z)) if np.isfinite(z).any() else 1
    fig = go.Figure(go.Heatmap(
        z=z, x=[str(c) for c in m.columns], y=[str(i) for i in m.index],
        colorscale="RdBu_r" if diverging else "Magma_r", zmid=0 if diverging else None,
        zmin=-vmax if diverging else 0, zmax=vmax if diverging else None, hoverongaps=False,
        colorbar=dict(title=colorbar), hovertemplate="%{y}<br>%{x}<br>value: %{z:.2f}<extra></extra>"))
    fig.update_layout(template=TEMPLATE, title=title, height=max(350, 28 * m.shape[0] + 200),
                      xaxis=dict(tickangle=-50), margin=dict(b=160))
    return fig


def strip_patients(donors, gene, title=None):
    d = donors.copy()
    d["status"] = np.where(d["expressing"], "expressing patient", "not expressing")
    fig = px.strip(d, x="study", y="pct_expressing", color="status", template=TEMPLATE,
                   color_discrete_map={"expressing patient": "#c0392b", "not expressing": "#999999"},
                   hover_data=[c for c in ["donor", "n_cells", "n_expressing"] if c in d.columns], height=520,
                   title=title or f"Patient-level {gene} detection by study (one point = one patient)")
    fig.update_layout(yaxis_title=f"% of a patient's cells expressing {gene}", xaxis_title="")
    return fig


def bar_by_study(inv, value, gene, title, color="tissue"):
    d = inv.sort_values(value)
    fig = px.bar(d, x=value, y="study", color=color if color in d.columns else None, orientation="h",
                 template=TEMPLATE, title=title, height=max(320, 28 * len(d) + 120),
                 hover_data=[c for c in ["n_cells", "n_cells_expressing", "n_expressing", "assay"] if c in d.columns])
    fig.update_layout(xaxis_title=f"% cells expressing {gene}", yaxis_title="")
    return fig


def save(fig, path):
    """Write a standalone interactive HTML file (opens in any browser)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(path), include_plotlyjs="cdn")
    return path
