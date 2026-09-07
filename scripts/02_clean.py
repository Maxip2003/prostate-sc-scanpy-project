import scanpy as sc

from utils import save_fig

adata = sc.read_h5ad("data/01_loaded.h5ad")

# Recompute QC metrics on raw counts (X in this dataset holds processed values)
adata.var["mt"] = adata.var["feature_name"].str.startswith("MT-")
adata_raw = adata.raw.to_adata()
adata_raw.var["mt"] = adata.var["mt"].values
sc.pp.calculate_qc_metrics(adata_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# Copy the QC columns back onto the main object, keeping everything else as-is
adata.obs["total_counts"] = adata_raw.obs["total_counts"].values
adata.obs["n_genes_by_counts"] = adata_raw.obs["n_genes_by_counts"].values
adata.obs["pct_counts_mt"] = adata_raw.obs["pct_counts_mt"].values

# Plot distributions to choose thresholds visually
save_fig(sc.pl.violin(adata, ["n_genes_by_counts"], jitter=0.3, show=False), "violin_qc_before.png")
save_fig(sc.pl.violin(adata, ["pct_counts_mt"], jitter=0.3, show=False), "violin_qc_2_before.png")

# Apply the QC thresholds decided from the plots above
n_before = adata.n_obs

adata = adata[adata.obs["n_genes_by_counts"] >= 500, :]
adata = adata[adata.obs["pct_counts_mt"] <= 20, :]

n_after = adata.n_obs
print(f"Cells before filtering: {n_before}")
print(f"Cells after filtering: {n_after}")
print(f"Removed: {n_before - n_after} cells ({100 * (n_before - n_after) / n_before:.1f}%)")

# Save the cleaned dataset for the next step
adata.write("data/02_cleaned.h5ad")
print("Saved to data/02_cleaned.h5ad")