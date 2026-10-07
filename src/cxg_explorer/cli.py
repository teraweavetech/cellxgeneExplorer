"""Command line: `cxg-explorer app | setup-r | example | version`."""
import argparse
import subprocess
import sys
from pathlib import Path


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cxg-explorer", description="CELLxGENE gene explorer")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("app", help="open the browser app")
    p.add_argument("--port", type=int, default=8501)
    p = sub.add_parser("setup-r", help="install the R package nebula (needed for NEBULA differential expression)")
    p.add_argument("--rscript", default=None, help="path to Rscript (auto-detected if omitted)")
    p = sub.add_parser("example", help="write a small synthetic dataset for trying the tools")
    p.add_argument("--out", default="example_data/example_pbmc.h5ad")
    sub.add_parser("version")
    a = ap.parse_args(argv)
    if a.cmd == "app":
        app = Path(__file__).with_name("app.py")
        theme = ["--theme.base", "light", "--theme.primaryColor", "#0F766E", "--theme.backgroundColor", "#F4F6F8",
                 "--theme.secondaryBackgroundColor", "#EEF2F4", "--theme.textColor", "#10212B", "--client.toolbarMode", "minimal"]
        sys.exit(subprocess.call([sys.executable, "-m", "streamlit", "run", str(app), "--server.port", str(a.port),
                                  "--browser.gatherUsageStats", "false"] + theme))
    if a.cmd == "setup-r":
        from .agent import find_rscript, r_install_nebula, r_status
        rs = find_rscript(a.rscript)
        if not rs:
            sys.exit("Rscript not found. Install R (https://cran.r-project.org) or pass --rscript PATH.")
        print(f"Using {rs}. Installing nebula (this can take several minutes) ...")
        ok, log = r_install_nebula(rs)
        print(log[-1500:])
        print(r_status(rs)[1])
        sys.exit(0 if ok else 1)
    if a.cmd == "example":
        from .example import make_example
        print(f"Example dataset written to {make_example(a.out)}  (gene to try: EXAMPLE1)")
    if a.cmd == "version":
        from . import __version__
        print(__version__)


if __name__ == "__main__":
    main()
