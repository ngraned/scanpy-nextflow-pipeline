import sys
import scanpy as sc
import warnings
import shutil
import os

warnings.filterwarnings('ignore')

matrix_dir = sys.argv[1]

print("=" * 60)
print("SCANPY PBMC ANALYSIS")
print("=" * 60)

print(f"\nLoading 10x data from: {matrix_dir}")
adata = sc.read_10x_mtx(matrix_dir)
print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

print("Quality control...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(f"After QC: {adata.n_obs} cells")

print("Normalization...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print("HVG selection...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

print("PCA...")
sc.pp.pca(adata)

print("Neighbors & UMAP...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)

print("Leiden clustering...")
sc.tl.leiden(adata, resolution=0.7)
n_clusters = len(adata.obs['leiden'].unique())

print("Saving results...")
adata.write_h5ad('adata_processed.h5ad')

# Save UMAP to current directory (not figures/)
sc.pl.umap(adata, color='leiden', save=False, show=False)
import matplotlib.pyplot as plt
plt.savefig('umap.png', dpi=100, bbox_inches='tight')
plt.close()

with open('report.txt', 'w') as f:
    f.write("PBMC ANALYSIS REPORT\n")
    f.write("=" * 60 + "\n")
    f.write(f"Cells: {adata.n_obs}\n")
    f.write(f"Genes: {adata.n_vars}\n")
    f.write(f"Clusters: {n_clusters}\n")

print("=" * 60)
print("COMPLETE!")
print("=" * 60)
