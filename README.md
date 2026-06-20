# Scanpy Nextflow Pipeline

A production-grade Nextflow pipeline for single-cell RNA-seq analysis using Scanpy.

## Features

- 10x Genomics data ingestion
- Quality control and filtering
- Normalization and HVG selection
- Dimensionality reduction (PCA, UMAP)
- Leiden clustering
- Publication-quality visualizations

## Installation

```bash
conda activate scanpy-latest
nextflow run main.nf --input_mtx ./data/filtered_gene_bc_matrices/hg19/ --outdir results
```

## Results

- `adata_processed.h5ad` - Processed single-cell object
- `umap.png` - UMAP visualization with clusters
- `report.txt` - Analysis summary

## Author

Your Name
