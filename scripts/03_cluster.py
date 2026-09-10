import argparse
import matplotlib.pyplot as plt
import scanpy as sc
import yaml
from utils import save_fig

parser = argparse.ArgumentParser()
parser.add_argument("--config", default="config/params.yaml")
args = parser.parse_args()

with open(args.config) as f:
    params = yaml.safe_load(f)

output_dir = params["dataset"]["output_dir"]

adata = sc.read_h5ad(f"{output_dir}/02_cleaned.h5ad")

adata = adata.raw.to_adata()
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")

sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False)
save_fig(plt.gcf(), "pca_variance_ratio_elbow.png")

adata.write(f"{output_dir}/03_clustered.h5ad")
print(f"Saved to {output_dir}/03_clustered.h5ad")