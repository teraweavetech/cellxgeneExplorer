"""Apply the new app style to an existing cellxgene-gene-explorer folder.

Usage (from anywhere):
    python apply_style_update.py PATH/TO/cellxgene-gene-explorer

What it does:
  - updates src/cxg_explorer/app.py and src/cxg_explorer/cli.py, KEEPING the copyright header lines
    at the top of each file exactly as they are in your copy (your name and year);
  - adds src/cxg_explorer/theme.py with the same header as your app.py;
  - saves your previous app.py and cli.py as app.py.bak and cli.py.bak.
No other file is touched.
"""
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent / "files"


def header_of(path):
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    out = []
    for ln in lines:
        if ln.startswith("#"):
            out.append(ln)
        else:
            break
    return "".join(out)


def main():
    repo = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    pkg = repo / "src" / "cxg_explorer"
    if not (pkg / "app.py").exists():
        sys.exit(f"Could not find src/cxg_explorer/app.py under {repo}. Pass the cellxgene-gene-explorer folder.")
    app_header = header_of(pkg / "app.py")
    for name in ["app.py", "cli.py"]:
        target = pkg / name
        head = header_of(target)
        shutil.copy2(target, target.with_name(name + ".bak"))
        target.write_text(head + (HERE / name).read_text(encoding="utf-8"), encoding="utf-8")
        print(f"updated  {target.relative_to(repo)}  (header kept, backup: {name}.bak)")
    theme = pkg / "theme.py"
    theme.write_text(app_header + (HERE / "theme.py").read_text(encoding="utf-8"), encoding="utf-8")
    print(f"added    {theme.relative_to(repo)}")
    print("Done. Start the app with:  cxg-explorer app")


if __name__ == "__main__":
    main()
