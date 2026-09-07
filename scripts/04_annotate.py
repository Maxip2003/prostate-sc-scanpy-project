import scanpy as sc
import pandas as pd

from utils import save_fig

adata = sc.read_h5ad("data/03_clustered.h5ad")

# Build a neighborhood graph using the first 15 PCs chosen from the elbow plot
sc.pp.neighbors(adata, n_pcs=15)

# Compute UMAP coordinates for visualization
sc.tl.umap(adata)

# Cluster cells into groups based on their neighborhood graph
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)

# Visualize clusters alongside the cell types already annotated in this dataset
ax = sc.pl.umap(adata, color=["leiden", "cell_type"], show=False)
save_fig(ax, "umap_clusters.png")

adata.write("data/04_annotated.h5ad")
print("Saved to data/04_annotated.h5ad")

crosstab = pd.crosstab(adata.obs["leiden"], adata.obs["cell_type"])
print(crosstab)

# Sanity check: KLK3 should be high in luminal cells
ax = sc.pl.violin(adata, ["KLK3"], groupby="cell_type", rotation=90, show=False)
save_fig(ax, "violin_klk3_check.png")

