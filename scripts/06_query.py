import argparse
import duckdb

parser = argparse.ArgumentParser(description="Query gene expression by cell cluster")
parser.add_argument("--gene", required=True, help="Gene symbol, e.g. PCA3")
args = parser.parse_args()

con = duckdb.connect("databases/prostate.duckdb")

result = con.execute("""
    SELECT
        c.leiden_label,
        c.dominant_cell_type,
        c.n_cells AS cluster_size,
        g.pct_cells_expressing,
        g.mean_expression
    FROM gene_expression_summary g
    JOIN clusters c ON g.cluster_id = c.id
    WHERE g.gene_symbol = ?
    ORDER BY g.pct_cells_expressing DESC
""", [args.gene]).fetchdf()

con.close()

if result.empty:
    print(f"No data found for gene '{args.gene}'")
else:
    print(f"\nExpression of {args.gene} by cluster:\n")
    print(result.to_string(index=False))