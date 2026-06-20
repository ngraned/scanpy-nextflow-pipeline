import scanpy as sc
import warnings
import os
import sys

warnings.filterwarnings('ignore')

def main(input_file):
    print("=" * 50)
    print("SCANPY ANALYSIS PIPELINE")
    print("=" * 50)
    
    os.makedirs('figures', exist_ok=True)
    
    print("Loading data...")
    adata = sc.read_h5ad(input_file)
    print(f"Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
    
    print("Quality control...")
    print(f"Before QC: {adata.n_obs} cells")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"After QC: {adata.n_obs} cells")
    
    print("Normalization...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    
    print("Dimensionality reduction...")
    sc.pp.pca(adata)
    print("Computing neighbors...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    print("Computing UMAP...")
    sc.tl.umap(adata)
    
    print("Clustering with Leiden...")
    sc.tl.leiden(adata, resolution=0.7, flavor='igraph', n_iterations=2, directed=False)
    n_clusters = len(adata.obs['leiden'].unique())
    print(f"Found {n_clusters} clusters")
    
    print("Generating visualizations...")
    sc.settings.figdir = 'figures'
    sc.pl.umap(adata, color='leiden', save='umap.png', show=False)
    
    print("Saving results...")
    adata.write_h5ad('adata_processed.h5ad')
    
    with open('clustering_report.txt', 'w') as f:
        f.write("SCANPY ANALYSIS REPORT\n")
        f.write("=" * 50 + "\n")
        f.write(f"Cells: {adata.n_obs}\n")
        f.write(f"Genes: {adata.n_vars}\n")
        f.write(f"Clusters: {n_clusters}\n")
    
    print("=" * 50)
    print("PIPELINE COMPLETE!")
    print("=" * 50)

if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else "pbmc_data.h5ad"
    main(input_file)
