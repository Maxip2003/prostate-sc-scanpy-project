import argparse
import scanpy as sc
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("--config", default="config/params.yaml")
args = parser.parse_args()

with open(args.config) as f:
    params = yaml.safe_load(f)

output_dir = params["dataset"]["output_dir"]

adata = sc.read_h5ad(params["dataset"]["input_path"])

# Use gene symbols instead of Ensembl IDs as var_names, so downstream
# scripts can reference genes of interest (KLK3, PCA3...) directly.
# .raw has its own var_names independent of adata.var, so rename both.
raw = adata.raw.to_adata()
raw.var_names = raw.var["feature_name"].values
raw.var_names_make_unique()
raw.var.index.name = None
adata.raw = raw

adata.var_names = adata.var["feature_name"].values
adata.var_names_make_unique()
adata.var.index.name = None

print(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

adata.write(f"{output_dir}/01_loaded.h5ad")
print(f"Saved to {output_dir}/01_loaded.h5ad")