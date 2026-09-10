import argparse
import scanpy as sc
import pandas as pd
import yaml
from utils import save_fig

parser = argparse.ArgumentParser()
parser.add_argument("--config", default="config/params.yaml")
args = parser.parse_args()

with open(args.config) as f:
    params = yaml.safe_load(f)

output_dir = params["dataset"]["output_dir"]
n_pcs = params["clustering"]["n_pcs"]
resolution = params["clustering"]["leiden_resolution"]

adata = sc.read_h5ad(f"{output_dir}/03_clustered.h5ad")

sc.pp.neighbors(adata, n_pcs=n_pcs)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2, directed=False)

ax = sc.pl.umap(adata, color=["leiden", "cell_type"], show=False)
save_fig(ax, "umap_clusters.png")

adata.write(f"{output_dir}/04_annotated.h5ad")
print(f"Saved to {output_dir}/04_annotated.h5ad")

crosstab = pd.crosstab(adata.obs["leiden"], adata.obs["cell_type"])
print(crosstab)

positive_control_gene = params["sanity_check"]["positive_control_gene"]
ax = sc.pl.violin(adata, [positive_control_gene], groupby="cell_type", rotation=90, show=False)
save_fig(ax, f"violin_{positive_control_gene.lower()}_check.png")