import matplotlib.pyplot as plt
import scanpy as sc

from utils import save_fig

adata = sc.read_h5ad("data/02_cleaned.h5ad")

# Start from raw counts, not the pre-processed X from the original file
adata = adata.raw.to_adata()

# Keep a copy of raw counts before transforming anything
adata.layers["counts"] = adata.X.copy()

# Normalize: make every cell's total count comparable, then log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find the genes that vary the most across cells — these carry the real signal
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Save the full normalized data before subsetting to HVGs only
adata.raw = adata
adata = adata[:, adata.var.highly_variable].copy()

# Scale so every gene contributes equally to the PCA
sc.pp.scale(adata, max_value=10)

# Reduce dimensions
sc.tl.pca(adata, svd_solver="arpack")

# Elbow plot: helps decide how many PCs actually carry signal
sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False)
save_fig(plt.gcf(), "pca_variance_ratio_elbow.png")

adata.write("data/03_clustered.h5ad")
print("Saved to data/03_clustered.h5ad")