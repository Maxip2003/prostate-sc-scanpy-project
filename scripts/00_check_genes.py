import argparse
import scanpy as sc
import yaml

# Accept a config file from the command line, defaulting to the real dataset —
# lets this script run against either the real or the test config without edits
parser = argparse.ArgumentParser()
parser.add_argument("--config", default="config/params.yaml")
args = parser.parse_args()

with open(args.config) as f:
    params = yaml.safe_load(f)

# This script always reads the ORIGINAL CELLxGENE file, never a processed one —
# var_names here are still Ensembl IDs, so we can look genes up directly by code
adata = sc.read_h5ad(params["dataset"]["input_path"])

with open("config/genes.yaml") as f:
    ids = yaml.safe_load(f)["genes"]

genes = list(ids.keys())

# Sanity check 1: confirm the genes are even listed in this dataset's annotation
found = adata.var[adata.var["feature_name"].isin(genes)]
print(found)

# Sanity check 2: confirm they're actually detected in cells, not just annotated
print("\n--- detection per cell ---")
for name, ensg in ids.items():
    values = adata.raw[:, ensg].X
    values = values.toarray() if hasattr(values, "toarray") else values
    n = (values > 0).sum()
    pct = 100 * n / adata.n_obs
    print(f"{name}: {n} of {adata.n_obs} cells ({pct:.1f}%)")