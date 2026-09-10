import argparse
import subprocess
from datetime import datetime

import duckdb
import pandas as pd
import scanpy as sc
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("--config", default="config/params.yaml")
args = parser.parse_args()

with open(args.config) as f:
    params = yaml.safe_load(f)

with open("config/genes.yaml") as f:
    gene_ids = yaml.safe_load(f)["genes"]

d = params["dataset"]
q = params["qc"]
c = params["clustering"]
output_dir = d["output_dir"]

adata = sc.read_h5ad(f"{output_dir}/04_annotated.h5ad")

con = duckdb.connect(d["db_path"])

con.execute("DELETE FROM gene_expression_summary")
con.execute("DELETE FROM clusters")
con.execute("DELETE FROM qc_metrics")
con.execute("DELETE FROM datasets")
con.execute("DELETE FROM runs")

n_cells_raw = sc.read_h5ad(f"{output_dir}/01_loaded.h5ad", backed="r").n_obs

dataset_row = pd.DataFrame([{
    "id": 1, "accession": d["accession"], "n_cells_raw": n_cells_raw,
    "n_cells_after_qc": adata.n_obs, "download_date": d["download_date"],
    "source_url": d["source_url"],
}])
con.execute("INSERT INTO datasets SELECT * FROM dataset_row")

qc_row = pd.DataFrame([{
    "id": 1, "dataset_id": 1, "min_genes_threshold": q["min_genes"],
    "max_pct_mt_threshold": q["max_pct_mt"],
    "cells_removed": n_cells_raw - adata.n_obs,
    "pct_removed": round(100 * (n_cells_raw - adata.n_obs) / n_cells_raw, 1),
}])
con.execute("INSERT INTO qc_metrics SELECT * FROM qc_row")

crosstab = pd.crosstab(adata.obs["leiden"], adata.obs["cell_type"])

cluster_rows = []
cluster_id_map = {}
next_id = 1
for leiden_label, row in crosstab.iterrows():
    n_cells = int(row.sum())
    dominant_type = row.idxmax()
    dominant_pct = round(100 * row.max() / n_cells, 1)
    cluster_rows.append({
        "id": next_id, "dataset_id": 1, "leiden_label": leiden_label,
        "n_cells": n_cells, "dominant_cell_type": dominant_type,
        "dominant_cell_type_pct": dominant_pct,
    })
    cluster_id_map[leiden_label] = next_id
    next_id += 1

clusters_df = pd.DataFrame(cluster_rows)
con.execute("INSERT INTO clusters SELECT * FROM clusters_df")

expr_rows = []
next_expr_id = 1
for gene, ensembl_id in gene_ids.items():
    match = adata.raw.var.loc[adata.raw.var["feature_name"] == gene]
    if match.empty:
        print(f"Warning: {gene} not found, skipping")
        continue
    row_key = match.index[0]
    values = adata.raw[:, row_key].X
    values = values.toarray().flatten() if hasattr(values, "toarray") else values.flatten()

    for leiden_label in adata.obs["leiden"].cat.categories:
        cell_mask = (adata.obs["leiden"] == leiden_label).values
        gene_values = values[cell_mask]
        pct_expr = round(100 * (gene_values > 0).sum() / len(gene_values), 2)
        mean_expr = round(float(gene_values.mean()), 4)
        expr_rows.append({
            "id": next_expr_id, "cluster_id": cluster_id_map[leiden_label],
            "gene_symbol": gene, "ensembl_id": ensembl_id,
            "pct_cells_expressing": pct_expr, "mean_expression": mean_expr,
        })
        next_expr_id += 1

expr_df = pd.DataFrame(expr_rows)
con.execute("INSERT INTO gene_expression_summary SELECT * FROM expr_df")

git_commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()

run_row = pd.DataFrame([{
    "id": 1, "timestamp": datetime.now().isoformat(timespec="seconds"),
    "git_commit": git_commit, "scanpy_version": sc.__version__,
    "n_pcs_used": c["n_pcs"], "leiden_resolution": c["leiden_resolution"],
}])
con.execute("INSERT INTO runs SELECT * FROM run_row")

con.close()
print(f"Data loaded into {d['db_path']}")