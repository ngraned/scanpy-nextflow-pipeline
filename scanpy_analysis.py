import sys
import scanpy as sc
import warnings
import pandas as pd
import matplotlib.pyplot as plt

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

print("Finding marker genes...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

print("Saving results...")
adata.write_h5ad('adata_processed.h5ad')

# Save UMAP
print("Generating UMAP visualization...")
sc.pl.umap(adata, color='leiden', save=False, show=False)
plt.savefig('umap.png', dpi=100, bbox_inches='tight')
plt.close()

# Save marker genes visualization (without using figures/ directory)
print("Generating marker genes plot...")
sc.set_figure_params(figsize=(12, 8))
sc.pl.rank_genes_groups(adata, n_genes=5, save=False, show=False)
plt.savefig('marker_genes_plot.png', dpi=100, bbox_inches='tight')
plt.close()

# Extract marker genes using built-in function
result = sc.get.rank_genes_groups_df(adata, group=None)

# Get top 5 genes per cluster
top_genes_list = []
for cluster in sorted(adata.obs['leiden'].unique()):
    cluster_result = result[result['group'] == str(cluster)].head(5)
    genes = ', '.join(cluster_result['names'].tolist())
    top_genes_list.append({'Cluster': cluster, 'Top_Markers': genes})

markers_df = pd.DataFrame(top_genes_list)

# Write comprehensive report
with open('report.txt', 'w') as f:
    f.write("PBMC ANALYSIS REPORT\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"Cells: {adata.n_obs}\n")
    f.write(f"Genes: {adata.n_vars}\n")
    f.write(f"Clusters: {n_clusters}\n\n")
    f.write("TOP 5 MARKER GENES PER CLUSTER\n")
    f.write("-" * 60 + "\n")
    for _, row in markers_df.iterrows():
        f.write(f"Cluster {row['Cluster']}: {row['Top_Markers']}\n")

# Save marker genes as CSV
markers_df.to_csv('marker_genes.csv', index=False)

print("=" * 60)
print("COMPLETE!")
print("=" * 60)
print("\nMarker genes identified per cluster:")
print(markers_df.to_string(index=False))
