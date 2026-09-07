import scanpy as sc
import pandas as pd
import duckdb
import yaml
from datetime import datetime

# --- Load pipeline outputs ---
adata = sc.read_h5ad("data/04_annotated.h5ad")

with open("config/genes.yaml") as f:
    gene_names = yaml.safe_load(f)["genes"]

con = duckdb.connect("databases/prostate.duckdb")

# --- 1. datasets ---
dataset_row = pd.DataFrame([{
    "id": 1,
    "accession": "4d0a653a-291d-44f6-966d-c3a7f1f3bd09",
    "n_cells_raw": 268,
    "n_cells_after_qc": adata.n_obs,
    "download_date": "19/08/2026",
    "source_url": "https://cellxgene.cziscience.com/collections/bdac7a53-fe34-4f04-8c46-f9bd5297c099"
}])
con.execute("INSERT INTO datasets SELECT * FROM dataset_row")

# --- 2. qc_metrics ---
qc_row = pd.DataFrame([{
    "id": 1,
    "dataset_id": 1,
    "min_genes_threshold": 500,
    "max_pct_mt_threshold": 20.0,
    "cells_removed": 8,
    "pct_removed": 3.0,
}])
con.execute("INSERT INTO qc_metrics SELECT * FROM qc_row")

# --- 3. clusters ---
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

# --- 4. gene_expression_summary ---
expr_rows = []
next_expr_id = 1
for gene in gene_names:
    match = adata.raw.var.loc[adata.raw.var["feature_name"] == gene]
    if match.empty:
        print(f"Warning: {gene} not found, skipping")
        continue
    gene_id = match.index[0]

    values = adata.raw[:, gene_id].X
    values = values.toarray().flatten() if hasattr(values, "toarray") else values.flatten()

    for leiden_label in adata.obs["leiden"].cat.categories:
        cell_mask = (adata.obs["leiden"] == leiden_label).values
        gene_values = values[cell_mask]
        pct_expr = round(100 * (gene_values > 0).sum() / len(gene_values), 2)
        mean_expr = round(float(gene_values.mean()), 4)

        expr_rows.append({
            "id": next_expr_id, "cluster_id": cluster_id_map[leiden_label],
            "gene_symbol": gene, "ensembl_id": gene_id,
            "pct_cells_expressing": pct_expr, "mean_expression": mean_expr,
        })
        next_expr_id += 1

expr_df = pd.DataFrame(expr_rows)
con.execute("INSERT INTO gene_expression_summary SELECT * FROM expr_df")

# --- 5. runs ---
run_row = pd.DataFrame([{
    "id": 1,
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "git_commit": None,               # TODO: fill with `git rev-parse --short HEAD`
    "scanpy_version": sc.__version__,
    "n_pcs_used": 15,
    "leiden_resolution": 1.0,
}])
con.execute("INSERT INTO runs SELECT * FROM run_row")

con.close()
print("Data loaded into databases/prostate.duckdb")