# Usage

## Browser app: `cxg-explorer app`

The sidebar has the results and download folders and three sections.

**Single study**
1. *Dataset*: search CELLxGENE by paper title, DOI or collection link, choose the collection and dataset, then **Download and open**. Or open an H5AD you already have.
2. *Analysis settings*: study folder name (defaults to first author and year), grouping column and groups, gene, optional cell populations, and the rule for calling a patient "expressing".
3. *Results*: tabs for Cellular, Cell-specific, Sample, Differential expression, How it's calculated, Methods, Descriptive report and Download (a zip of all results).

**Cross-study summary**: choose the results folder, **Scan**, pick the gene and studies, then **Build summary**. The *Cell-type report* tab builds reports for the cell types you select.

## Jupyter

`notebooks/01_gene_explorer.ipynb` and `notebooks/02_cross_study_summary.ipynb` open the same workflows as panels with buttons and drop-downs. Results appear in tabs below the panel; interactive plots are also saved as standalone HTML files in `interactive/`.

## Tips

- Use the paper title exactly as shown on CELLxGENE. If it isn't found, use the DOI or collection link.
- Gene names must be official symbols (CD274, not PD-L1). Close matches are suggested.
- Analysing a subset of groups writes to a separate folder (`..._normal_vs_COVID-19`), so runs don't overwrite each other.
- For the cross-study summary, analyse every study with the same gene. Re-running a study reuses the downloaded file.
- Edit `_cross_study/celltype_to_broad_class.csv` to change how cell types are grouped into broad classes, and `_cross_study/study_details_editable.csv` to add sampling details, notes or corrected citations. Then build again.
