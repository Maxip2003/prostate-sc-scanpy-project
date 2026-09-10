import argparse
import scanpy as sc
import yaml
from utils import save_fig

parser = argparse.ArgumentParser()
parser.add_argument("--config", default="config/params.yaml")
args = parser.parse_args()

with open(args.config) as f:
    params = yaml.safe_load(f)

output_dir = params["dataset"]["output_dir"]
min_genes = params["qc"]["min_genes"]
max_pct_mt = params["qc"]["max_pct_mt"]
max_genes = params["qc"].get("max_genes")

adata = sc.read_h5ad(f"{output_dir}/01_loaded.h5ad")

adata.var["mt"] = adata.var["feature_name"].str.startswith("MT-")
adata_raw = adata.raw.to_adata()
adata_raw.var["mt"] = adata.var["mt"].values
sc.pp.calculate_qc_metrics(adata_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

adata.obs["total_counts"] = adata_raw.obs["total_counts"].values
adata.obs["n_genes_by_counts"] = adata_raw.obs["n_genes_by_counts"].values
adata.obs["pct_counts_mt"] = adata_raw.obs["pct_counts_mt"].values

del adata_raw

qc_summary = adata.obs[["n_genes_by_counts", "pct_counts_mt"]].describe(
    percentiles=[.01, .05, .25, .5, .75, .95, .99]
)
print(qc_summary)
qc_summary.to_csv("figures/qc_summary.csv")

save_fig(sc.pl.violin(adata, ["n_genes_by_counts"], jitter=0.3, show=False), "violin_qc_before.png")
save_fig(sc.pl.violin(adata, ["pct_counts_mt"], jitter=0.3, show=False), "violin_qc_2_before.png")

n_before = adata.n_obs
adata = adata[adata.obs["n_genes_by_counts"] >= min_genes, :]
if max_genes is not None:
    adata = adata[adata.obs["n_genes_by_counts"] <= max_genes, :]
adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mt, :]

n_after = adata.n_obs
print(f"Cells before filtering: {n_before}")
print(f"Cells after filtering: {n_after}")
print(f"Removed: {n_before - n_after} cells ({100 * (n_before - n_after) / n_before:.1f}%)")

adata.write(f"{output_dir}/02_cleaned.h5ad")
print(f"Saved to {output_dir}/02_cleaned.h5ad")