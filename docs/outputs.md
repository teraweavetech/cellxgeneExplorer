# Output files

## Single study: `<results>/<study>/<GENE>__<dataset>/`

| Folder or file | Contents |
|---|---|
| `1_cellular_level/` | `overall_summary.csv`, `celltype_summary.csv`, `summary_by_<disease/tissue/sex...>.csv`, `per_cell_expression.csv`; UMAPs (cell types, expression, expressing vs not), bar, violin and dot plots |
| `2_cell_specific_level/<population>/` | population summary, breakdowns by covariate, UMAPs, and a `sample_level/` folder with per-patient results for that population |
| `3_sample_level/` | `subject_summary.csv`, `subjects_expressing_headline.csv`, `subjects_expressing_overview.csv` and `_by_threshold.csv` (number and proportion of patients), patient x cell type tables, `celltype_n_subjects_expressing.csv`, plots |
| `4_DE/` | `fisher_celltype_vs_rest.csv`, `fisher_between_groups.csv` (+ heatmap), `nebula_summary.csv`, `<cell type>/<mode>/nebula_de_results.csv` with a volcano plot (PNG and interactive HTML), `heatmap_top_de_genes.png` |
| `interactive/` | interactive HTML plots (open in any browser) |
| `report.html` | all tables and figures in one page |
| `HOW_IT_IS_CALCULATED.html` | every column explained with worked examples, plus automatic consistency checks |
| `METHODS.md` | methods text for this run |
| `DESCRIPTIVE_REPORT.html/.md` | narrative report (only if generated) |
| `study_info.json`, `samples_used.csv` | dataset metadata, citation information and samples used (read by the cross-study summary) |

## Cross-study: `<results>/_cross_study/<GENE>/`

| File | Contents |
|---|---|
| `<GENE>_cross_study_summary.html` / `.xlsx` | all-cells report and all tables |
| `inventory.csv`, `study_details.csv`, `samples_used.csv` | one row per study; citations, technology, samples, sampling information |
| `heatmap_broad_pct.csv`, `heatmap_fine_pct.csv` | % cells expressing per study x cell class / cell type |
| `consistency_broad.csv`, `consistency_fine.csv` | in how many studies the gene is expressed in each cell type |
| `patients_per_study.csv`, `patients_all.csv`, `groups_per_study.csv` | patient-level and group-level summaries |
| `celltype_report/` | `<GENE>_celltype_report.html` / `.xlsx` and one folder per selected cell type |
| `HOW_IT_IS_CALCULATED.html` | calculation guide for the summary |

Shared editable files in `_cross_study/`: `celltype_to_broad_class.csv`, `study_details_editable.csv`, and `citations_cache.json`.
