# Example data

Create a small synthetic dataset in the CELLxGENE format (no download needed):

```bash
cxg-explorer example --out example_data/example_pbmc.h5ad
```

Open it in the app or notebook with "Open a local H5AD file" and analyse the gene **EXAMPLE1**. It is enriched in classical monocytes, more so in COVID-19 than in normal donors, and genes G1-G5 are co-expressed with it (G6-G7 anti-correlated), so the differential expression step has a known answer to recover.
