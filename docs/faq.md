# FAQ and troubleshooting

**The paper isn't found.** Use the title exactly as on CELLxGENE, or paste the DOI or the collection link. Not every paper's data is on CELLxGENE.

**"Gene not found".** Use the official symbol as stored by CELLxGENE (for example CD274 for PD-L1). Close matches are suggested.

**Only Fisher results appear in the DE tab.** NEBULA did not run; the banner above the results says why. Usually Rscript isn't found (enter its path), nebula isn't installed (press **Install nebula**), or the dataset has no raw counts.

**NEBULA fails with `unused argument (max = ncore)` or a `with(plan(...))` error.** Your `parallelly` or `future` packages are too old. Run `cxg-explorer setup-r`, which upgrades them.

**A patient count says 0 of N although cells express the gene.** The patient rule needs both a minimum number of expressing cells and a minimum percentage. The sensitivity table (`subjects_expressing_overview.csv`) shows the answer under several thresholds.

**Why are expression levels not compared between studies?** Each study deposits differently normalized data, so the same value means different things in different studies. The summary compares detection proportions and patients instead.

**The same donors may appear in two studies.** CELLxGENE atlases often re-use cells from other datasets. The summary warns when a study contains re-used cells, or when donor IDs repeat across studies.

**Large datasets.** Expression is read from disk one gene at a time, so memory use stays low. Downloads can be several GB, so set the download folder on a drive with space; files are reused on later runs.
