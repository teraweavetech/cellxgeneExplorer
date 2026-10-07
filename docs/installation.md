# Installation

## Python package

```bash
git clone https://github.com/teraweavetech/cellxgeneExplorer.git
cd cellxgeneExplorer
pip install -e .
```

Check it works with `cxg-explorer version`.

Conda users: `conda env create -f environment.yml && conda activate cxg-explorer`. The environment includes JupyterLab and R.

## Internet access

The explorer needs to reach `api.cellxgene.cziscience.com` and `datasets.cellxgene.cziscience.com`. Citation lookup also uses `api.crossref.org`, but it can be switched off. If you are behind a proxy, set `HTTPS_PROXY` before starting Python.

## R and nebula (optional, for differential expression)

1. Install R 4.4 or newer from https://cran.r-project.org.
2. Run `cxg-explorer setup-r` (or `cxg-explorer setup-r --rscript /path/to/Rscript`), or press **Install nebula** in the app or notebook. This installs nebula and upgrades `parallelly`, `future`, `future.apply`, `globals` and `doFuture` if they are too old. Older versions of these packages satisfy nebula's declared requirements but fail at run time.
3. Press **Check R / nebula**. It runs a small test fit and should report success.

You can also run the installer directly in R: `Rscript src/cxg_explorer/r/install_nebula.R`.

Typical Rscript locations: Windows `C:\Program Files\R\R-4.4.x\bin\Rscript.exe`; macOS `/Library/Frameworks/R.framework/Resources/bin/Rscript` or `/opt/homebrew/bin/Rscript`; Linux `/usr/bin/Rscript`.

## gpt-oss:20b (optional, for the descriptive report)

Install [Ollama](https://ollama.com), then run `ollama pull gpt-oss:20b`. Ollama serves it at `http://localhost:11434`, the default in the app. Any OpenAI-compatible server for gpt-oss-20b also works; enter its URL ending in `/v1`.
