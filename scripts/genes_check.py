import sys
import scanpy as sc

if len(sys.argv) < 2:
    print("Usage: python scripts/genes_check.py <GENE_NAME>")
    sys.exit(1)

gene = sys.argv[1]

adata = sc.read_h5ad("data/04_annotated.h5ad")

# Look up this gene's actual index value inside THIS file,
# instead of relying on a manually copied Ensembl ID
match = adata.raw.var.loc[adata.raw.var["feature_name"] == gene]
if match.empty:
    print(f"Gene '{gene}' not found in this dataset.")
    sys.exit(1)

gene_id = match.index[0]

sc.pl.violin(adata, [gene_id], groupby="cell_type", rotation=90, save=f"_{gene}_check.png")
print(f"Saved figures/violin_{gene}_check.png")