import scanpy as sc
import yaml
from utils import save_fig

adata = sc.read_h5ad("data/01_loaded.h5ad")

# Load pipeline parameters from a single config file instead of hardcoding values,
# so switching datasets only means editing the YAML, not the code
with open("config/params.yaml") as f:
    params = yaml.safe_load(f)

min_genes = params["qc"]["min_genes"]
max_pct_mt = params["qc"]["max_pct_mt"]
max_genes = params["qc"].get("max_genes")  # optional, None if not set — keeps old configs working

# Recompute QC metrics on raw counts (X in this dataset holds processed values)
adata.var["mt"] = adata.var["feature_name"].str.startswith("MT-")
adata_raw = adata.raw.to_adata()
adata_raw.var["mt"] = adata.var["mt"].values
sc.pp.calculate_qc_metrics(adata_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# Copy the QC columns back onto the main object, keeping everything else as-is
adata.obs["total_counts"] = adata_raw.obs["total_counts"].values
adata.obs["n_genes_by_counts"] = adata_raw.obs["n_genes_by_counts"].values
adata.obs["pct_counts_mt"] = adata_raw.obs["pct_counts_mt"].values

del adata_raw  # free this copy, we already extracted what we needed

# Print exact percentiles to pick thresholds precisely, complementing the plots below
qc_summary = adata.obs[["n_genes_by_counts", "pct_counts_mt"]].describe(
    percentiles=[.01, .05, .25, .5, .75, .95, .99]
)
print(qc_summary)
qc_summary.to_csv("figures/qc_summary.csv")
print("Saved to figures/qc_summary.csv")

# Plot distributions to choose thresholds visually
save_fig(sc.pl.violin(adata, ["n_genes_by_counts"], jitter=0.3, show=False), "violin_qc_before.png")
save_fig(sc.pl.violin(adata, ["pct_counts_mt"], jitter=0.3, show=False), "violin_qc_2_before.png")

# Apply the QC thresholds decided from the plots and percentiles above
n_before = adata.n_obs

adata = adata[adata.obs["n_genes_by_counts"] >= min_genes, :]
if max_genes is not None:
    adata = adata[adata.obs["n_genes_by_counts"] <= max_genes, :]
adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mt, :]

n_after = adata.n_obs
print(f"Cells before filtering: {n_before}")
print(f"Cells after filtering: {n_after}")
print(f"Removed: {n_before - n_after} cells ({100 * (n_before - n_after) / n_before:.1f}%)")

# Save the cleaned dataset for the next step
adata.write("data/02_cleaned.h5ad")
print("Saved to data/02_cleaned.h5ad")