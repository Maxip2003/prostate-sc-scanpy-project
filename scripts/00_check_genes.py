import scanpy as sc
import yaml

with open("config/params.yaml") as f:
    params = yaml.safe_load(f)

adata = sc.read_h5ad(params["dataset"]["input_path"])

with open("config/genes.yaml") as f:
    ids = yaml.safe_load(f)["genes"]

genes = list(ids.keys())
found = adata.var[adata.var["feature_name"].isin(genes)]
print(found)

print("\n--- detection per cell ---")
for name, ensg in ids.items():
    values = adata.raw[:, ensg].X
    values = values.toarray() if hasattr(values, "toarray") else values
    n = (values > 0).sum()
    pct = 100 * n / adata.n_obs
    print(f"{name}: {n} of {adata.n_obs} cells ({pct:.1f}%)")