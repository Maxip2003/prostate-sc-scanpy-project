import subprocess
from datetime import datetime

import duckdb
import pandas as pd
import scanpy as sc
import yaml
from datetime import datetime

# --- Load pipeline outputs ---
adata = sc.read_h5ad("data/04_annotated.h5ad")

with open("config/genes.yaml") as f:
    gene_ids = yaml.safe_load(f)["genes"]  # {"KLK3": "ENSG...", "PCA3": "ENSG...", ...}

con = duckdb.connect("databases/prostate.duckdb")

# --- 0. Clear existing rows so this script can be re-run safely ---
con.execute("DELETE FROM gene_expression_summary")
con.execute("DELETE FROM clusters")
con.execute("DELETE FROM qc_metrics")
con.execute("DELETE FROM datasets")
con.execute("DELETE FROM runs")

# Load dataset metadata, QC thresholds, and clustering parameters from config instead of hardcoding them — this is what changes when we switch datasets
with open("config/params.yaml") as f:
    params = yaml.safe_load(f)
d = params["dataset"]
q = params["qc"]
c = params["clustering"]

# --- 1. datasets ---
dataset_row = pd.DataFrame([{
    "id": 1,
    "accession": d["accession"],
    "n_cells_raw": d["n_cells_raw"],
    "n_cells_after_qc": adata.n_obs,
    "download_date": d["download_date"],
    "source_url": d["source_url"],
}])
con.execute("INSERT INTO datasets SELECT * FROM dataset_row")

# --- 2. qc_metrics ---
qc_row = pd.DataFrame([{
    "id": 1,
    "dataset_id": 1,
    "min_genes_threshold": q["min_genes"],
    "max_pct_mt_threshold": q["max_pct_mt"],
    # cells_removed / pct_removed are derived, not stored in config,
    # so they stay correct automatically for any dataset size
    "cells_removed": d["n_cells_raw"] - adata.n_obs,
    "pct_removed": round(100 * (d["n_cells_raw"] - adata.n_obs) / d["n_cells_raw"], 1),
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
for gene, ensembl_id in gene_ids.items():
    # adata.raw.var_names mixes Ensembl codes and gene symbols depending on the
    # gene, so we look up the row by symbol to get expression values, but we
    # always store the Ensembl ID from config/genes.yaml (verified on Ensembl),
    # never whatever happens to be in the file's index.
    match = adata.raw.var.loc[adata.raw.var["feature_name"] == gene]
    if match.empty:
        print(f"Warning: {gene} not found, skipping")
        continue
    row_key = match.index[0]  # used only to locate the column, not stored

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

# --- 5. runs ---
git_commit = subprocess.check_output(
    ["git", "rev-parse", "--short", "HEAD"], text=True
).strip()

run_row = pd.DataFrame([{
    "id": 1,
    "timestamp": datetime.now().isoformat(timespec="seconds"),
    "git_commit": git_commit,
    "scanpy_version": sc.__version__,
    "n_pcs_used": c["n_pcs"],
    "leiden_resolution": c["leiden_resolution"],
}])
con.execute("INSERT INTO runs SELECT * FROM run_row")

con.close()
print("Data loaded into databases/prostate.duckdb")