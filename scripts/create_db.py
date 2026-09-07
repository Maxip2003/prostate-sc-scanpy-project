import duckdb

con = duckdb.connect("databases/prostate.duckdb")

with open("scripts/schema.sql") as f:
    con.execute(f.read())

print("Database created with all tables")
con.close()