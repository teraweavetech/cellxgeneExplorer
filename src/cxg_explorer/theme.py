"""Look and feel of the browser app (clinical style: white page, dark teal sidebar, IBM Plex type),
the "Made by" credit, and dropdown checklists used instead of multi-select boxes."""
import html
import re
from pathlib import Path

import streamlit as st

TEAL = "#0F766E"
TEAL_DARK = "#0B5C55"
INK = "#10212B"
MUTED = "#55646E"
LINE = "#DDE3E7"
GROUND = "#F4F6F8"
SIDEBAR = "#0B2530"
EXPRESSION_SCALE = [[0.0, "#B7EBE4"], [0.35, "#5CCBBE"], [0.65, "#2BB3A3"], [0.85, "#0F766E"], [1.0, "#0B4F4A"]]
CATEGORICAL = ["#0F766E", "#1D4ED8", "#B45309", "#7C3AED", "#BE185D", "#0369A1", "#4D7C0F", "#9F1239",
               "#475569", "#A16207"]

CSS = f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600;700&display=swap');
.stApp {{ background: {GROUND}; font-family: 'IBM Plex Sans', system-ui, sans-serif; color: {INK}; }}
.stApp h1, .stApp h2, .stApp h3, .stApp h4 {{ font-family: 'IBM Plex Sans', system-ui, sans-serif; color: {INK};
  letter-spacing: -0.01em; }}
.block-container {{ padding-top: 0.8rem; max-width: 1400px; }}

/* sidebar */
[data-testid="stSidebar"] {{ background: {SIDEBAR}; }}
[data-testid="stSidebar"] p, [data-testid="stSidebar"] label, [data-testid="stSidebar"] span,
[data-testid="stSidebar"] div {{ color: #E6F0F2; }}
[data-testid="stSidebar"] input {{ background: #0F2F3A !important; color: #E6F0F2 !important;
  font-family: 'IBM Plex Mono', monospace; font-size: 13px; }}
[data-testid="stSidebar"] [data-testid="stWidgetLabel"] p {{ color: #9DB4BC; font-size: 12px; text-transform: uppercase;
  letter-spacing: 0.06em; }}
[data-testid="stSidebar"] [role="radiogroup"] {{ gap: 4px; }}
[data-testid="stSidebar"] [role="radiogroup"] label {{ padding: 10px 14px; border-radius: 8px; width: 100%;
  margin: 0; cursor: pointer; }}
[data-testid="stSidebar"] [data-testid="stRadioOption"] > div > div:first-child:not([data-testid="stMarkdownContainer"])
  {{ display: none; }}
[data-testid="stSidebar"] [role="radiogroup"] label p {{ font-size: 15px; color: #C9D8DC; }}
[data-testid="stSidebar"] [data-testid="stRadioOption"][data-selected="true"] {{ background: #123844; }}
[data-testid="stSidebar"] [data-testid="stRadioOption"][data-selected="true"] p {{ color: #FFFFFF; font-weight: 600; }}
[data-testid="stSidebar"] [role="radiogroup"] label:has(input:checked) {{ background: #123844; }}
[data-testid="stSidebar"] [role="radiogroup"] label:has(input:checked) p {{ color: #FFFFFF; font-weight: 600; }}
[data-testid="stSidebar"] [role="radiogroup"] label:hover {{ background: #10303B; }}
[data-testid="stSidebar"] [data-testid="stTextInputRootElement"], [data-testid="stSidebar"] [data-baseweb="input"]
  {{ background: #0F2F3A; border: 1px solid #2B4C57 !important; border-radius: 6px; }}
[data-testid="stSidebar"] a {{ color: #5EEAD4; }}

[data-testid="stSidebar"] [data-testid="stTextInputRootElement"] {{ border-color: #2B4C57 !important; }}
[data-testid="stAppDeployButton"] {{ display: none; }}
[data-testid="stHeader"] {{ background: transparent; }}

/* buttons */
.stApp button[kind="primary"], .stApp [data-testid="stBaseButton-primary"] {{ background: {TEAL}; border: none;
  color: #FFFFFF; font-weight: 600; border-radius: 8px; }}
.stApp button[kind="primary"]:hover, .stApp [data-testid="stBaseButton-primary"]:hover {{ background: {TEAL_DARK};
  color: #FFFFFF; }}
.stApp button[kind="secondary"], .stApp [data-testid="stBaseButton-secondary"] {{ background: #FFFFFF;
  border: 1px solid #C8D1D7; color: {INK}; border-radius: 8px; font-weight: 500; }}
.stApp button[kind="secondary"]:hover, .stApp [data-testid="stBaseButton-secondary"]:hover {{ border-color: {TEAL};
  color: {TEAL}; }}

/* tabs */
.stTabs [data-baseweb="tab-list"] {{ gap: 4px; border-bottom: 1px solid {LINE}; }}
.stTabs [data-baseweb="tab"] {{ padding: 10px 14px; }}
.stTabs [data-baseweb="tab"] p {{ font-size: 14px; color: {MUTED}; }}
.stTabs [aria-selected="true"] p {{ color: {INK}; font-weight: 600; }}
.stTabs [data-baseweb="tab-highlight"] {{ background: {TEAL}; height: 3px; }}

/* panels and expanders */
[class*="st-key-panel"] {{ background: #FFFFFF; border: 1px solid {LINE}; border-radius: 12px; padding: 18px 20px; }}
[data-testid="stExpander"] details {{ background: #FFFFFF; border: 1px solid {LINE}; border-radius: 12px; }}
[data-testid="stExpander"] summary p {{ font-weight: 600; font-size: 15px; }}
[data-testid="stPopover"] button {{ justify-content: space-between; }}

/* custom blocks */
.gx-crumb {{ font-size: 13px; color: {MUTED}; margin-bottom: 4px; }}
.gx-title {{ font-size: 10px; font-weight: 400; letter-spacing: -0.02em; margin: 0; line-height: 1.15; color: {INK}; }}
.gx-sub {{ font-size: 15px; color: {MUTED}; margin-top: 6px; }}
.gx-cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 16px; margin: 18px 0 8px; }}
.gx-card {{ background: #FFFFFF; border: 1px solid {LINE}; border-radius: 12px; padding: 16px 20px; }}
.gx-card .l {{ font-size: 13px; color: {MUTED}; }}
.gx-card .v {{ font-family: 'IBM Plex Mono', monospace; font-size: 28px; font-weight: 500; color: {INK};
  display: flex; align-items: baseline; gap: 10px; flex-wrap: wrap; margin-top: 4px; }}
.gx-pill {{ font-family: 'IBM Plex Sans', sans-serif; font-size: 14px; font-weight: 600; color: {TEAL};
  background: #E3F5F2; padding: 3px 8px; border-radius: 999px; }}
.gx-chips {{ display: flex; gap: 8px; flex-wrap: wrap; margin-top: 8px; }}
.gx-chip {{ font-size: 14px; padding: 5px 10px; border-radius: 999px; border: 1px solid #C8D1D7; color: {INK}; }}
.gx-bars h3 {{ font-size: 18px; margin: 0 0 14px; }}
.gx-bar {{ margin-bottom: 14px; }}
.gx-bar .top {{ display: flex; justify-content: space-between; gap: 12px; font-size: 14px; }}
.gx-bar .top span:last-child {{ font-family: 'IBM Plex Mono', monospace; }}
.gx-bar .track {{ height: 10px; background: #EEF2F4; border-radius: 5px; margin: 6px 0 4px; }}
.gx-bar .fill {{ height: 10px; background: {TEAL}; border-radius: 5px; }}
.gx-bar .n {{ font-size: 12px; color: {MUTED}; }}
.gx-brand {{ display: flex; align-items: center; gap: 12px; margin: 4px 0 18px; }}
.gx-brand .t {{ font-weight: 700; font-size: 18px; color: #FFFFFF; }}
.gx-brand .s {{ font-size: 12px; color: #9DB4BC; }}
.gx-credit {{ border-top: 1px solid #23434E; margin-top: 22px; padding-top: 16px; font-size: 13px; line-height: 1.7; }}
.gx-credit .m {{ color: #9DB4BC; }}
</style>
"""

BRAND_SVG = ('<svg width="34" height="34" viewBox="0 0 34 34" fill="none" stroke="#5EEAD4" stroke-width="2" '
             'aria-hidden="true"><circle cx="17" cy="17" r="14"/><circle cx="12" cy="14" r="3"/>'
             '<circle cx="21" cy="12" r="2"/><circle cx="20" cy="22" r="3.5"/></svg>')


def inject():
    st.markdown(CSS, unsafe_allow_html=True)


def _md(s):
    st.markdown(s, unsafe_allow_html=True)


def esc(x):
    return html.escape(str(x))


# ----------------------------------------------------------------------------- credit ("Made by")
def _placeholder(s):
    return not s or ("[" in s and "]" in s)


def credits():
    """Author name and repository link, read from CITATION.cff (or the installed package metadata)."""
    name, url = None, None
    cff = Path(__file__).resolve().parents[2] / "CITATION.cff"
    if cff.exists():
        t = cff.read_text(encoding="utf-8", errors="ignore")
        m = re.search(r"^authors:\s*\n((?:[ \t].*\n?|-.*\n?)*)", t, re.M)
        if m:
            names = []
            for e in re.split(r"^\s*-\s+", m.group(1), flags=re.M):
                kv = dict(re.findall(r'^\s*([\w-]+):\s*"?([^"\n]*?)"?\s*$', e, re.M))
                nm = kv.get("name") or " ".join(x for x in [kv.get("given-names"), kv.get("family-names")] if x)
                if nm and not _placeholder(nm):
                    names.append(nm)
            name = ", ".join(names) or None
        u = re.search(r'^repository-code:\s*"?([^"\n]+)"?', t, re.M)
        url = u.group(1).strip() if u and not _placeholder(u.group(1)) else None
    if not name or not url:
        try:
            from importlib.metadata import metadata
            md = metadata("cellxgene-gene-explorer")
            if not name:
                a = md.get("Author") or (md.get("Author-email") or "").split("<")[0].strip().strip('"')
                name = None if _placeholder(a) else a or None
            if not url:
                for v in md.get_all("Project-URL") or []:
                    if "github.com" in v and not _placeholder(v):
                        url = v.split(",", 1)[-1].strip()
                        break
        except Exception:
            pass
    return name, url


def sidebar_brand():
    st.sidebar.markdown(f'<div class="gx-brand"><div><div class="t">What do you want to Analyze</div>', unsafe_allow_html=True)


def sidebar_credit(version):
    name, url = credits()
    rows = []
    if name:
        rows.append(f"Made by <strong>{esc(name)}</strong>")
    if url:
        rows.append(f'<a href="{esc(url)}" target="_blank">{esc(re.sub(r"^https?://", "", url))}</a>')
    rows.append(f'<span class="m">v{esc(version)} · PolyForm Noncommercial</span>')
    st.sidebar.markdown('<div class="gx-credit">' + "<br>".join(rows) + "</div>", unsafe_allow_html=True)


# ----------------------------------------------------------------------------- blocks
def page_header(title, crumb=None, subtitle=None):
    _md((f'<div class="gx-crumb">{esc(crumb)}</div>' if crumb else "")
        + f'<h1 class="gx-title">{esc(title)}</h1>'
        + (f'<div class="gx-sub">{esc(subtitle)}</div>' if subtitle else ""))


def metric_cards(items):
    """items: dicts with label, and value (+ optional pill) or chips (list of strings)."""
    out = []
    for it in items:
        if it.get("chips") is not None:
            body = '<div class="gx-chips">' + "".join(f'<span class="gx-chip">{esc(c)}</span>'
                                                      for c in it["chips"]) + "</div>"
        else:
            pill = f'<span class="gx-pill">{esc(it["pill"])}</span>' if it.get("pill") else ""
            body = f'<div class="v"><span>{esc(it["value"])}</span>{pill}</div>'
        out.append(f'<div class="gx-card"><div class="l">{esc(it["label"])}</div>{body}</div>')
    _md('<div class="gx-cards">' + "".join(out) + "</div>")


def bars_panel(title, rows, max_rows=8):
    """rows: list of (name, pct, n_expressing, n_cells), highest first."""
    rows = list(rows)[:max_rows]
    top = max([r[1] for r in rows] + [1e-9])
    parts = [f'<div class="gx-bars"><h3>{esc(title)}</h3>']
    for name, pct, k, n in rows:
        parts.append(f'<div class="gx-bar"><div class="top"><span>{esc(name)}</span><span>{pct:.1f}%</span></div>'
                     f'<div class="track"><div class="fill" style="width:{100 * pct / top:.1f}%"></div></div>'
                     f'<div class="n">{int(k):,} of {int(n):,} cells</div></div>')
    _md("".join(parts) + "</div>")


_panel_n = [0]


def panel():
    """A white, rounded box (falls back to a plain container on older Streamlit)."""
    _panel_n[0] += 1
    try:
        return st.container(border=False, key=f"panel_{_panel_n[0]}")
    except TypeError:
        return st.container(border=True)


def reset_panels():
    _panel_n[0] = 0


def style_fig(fig):
    """Match plotly figures to the app: fonts, white panels, teal expression scale."""
    fig.update_layout(font=dict(family="IBM Plex Sans, system-ui, sans-serif", color=INK, size=13),
                      paper_bgcolor="#FFFFFF", plot_bgcolor="#FAFBFC", colorway=CATEGORICAL,
                      title=dict(font=dict(size=16)), margin=dict(t=56, l=10, r=10, b=40))
    for tr in fig.data:
        mk = getattr(tr, "marker", None)
        if tr.type in ("scattergl", "scatter") and mk is not None and mk.colorscale is not None \
                and mk.colorbar is not None and mk.colorbar.title.text:
            tr.update(marker=dict(colorscale=EXPRESSION_SCALE))
    return fig


# ----------------------------------------------------------------------------- dropdown checklists
def _opt_key(key, o):
    return f"{key}__{re.sub(r'[^A-Za-z0-9_]+', '_', str(o))}"


def dropdown_checklist(label, options, key, default=None, fmt=str, help=None):
    """A dropdown menu with a tick box per option (instead of a multi-select with removable chips).
    Returns the ticked options in their original order."""
    options = list(options)
    default = options if default is None else list(default)
    ss = st.session_state
    for o in options:
        k = _opt_key(key, o)
        if k not in ss:
            ss[k] = o in default
    sel = [o for o in options if ss.get(_opt_key(key, o))]
    if not options:
        summary = "none available"
    elif len(sel) == len(options):
        summary = "all selected" if len(options) > 1 else fmt(sel[0])
    elif not sel:
        summary = "none selected"
    elif len(sel) == 1:
        summary = fmt(sel[0])
    else:
        summary = f"{len(sel)} of {len(options)} selected"

    def set_all(v):
        for o in options:
            ss[_opt_key(key, o)] = v

    st.markdown(f'<div style="font-size:14px;margin-bottom:4px">{esc(label)}</div>', unsafe_allow_html=True)
    title = summary if len(summary) <= 70 else summary[:67] + "..."
    if hasattr(st, "popover"):
        try:
            box = st.popover(title, help=help, width="stretch")
        except TypeError:
            box = st.popover(title, help=help, use_container_width=True)
    else:
        box = st.expander(title)
    with box:
        if len(options) > 1:
            a, b = st.columns(2)
            a.button("Select all", key=f"{key}__all", on_click=set_all, args=(True,), use_container_width=True)
            b.button("Clear", key=f"{key}__none", on_click=set_all, args=(False,), use_container_width=True)
        for o in options:
            st.checkbox(fmt(o), key=_opt_key(key, o))
    return [o for o in options if ss.get(_opt_key(key, o))]
