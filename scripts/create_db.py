import argparse
import duckdb
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("--config", default="config/params.yaml")
args = parser.parse_args()

with open(args.config) as f:
    params = yaml.safe_load(f)

db_path = params["dataset"]["db_path"]

con = duckdb.connect(db_path)

with open("scripts/schema.sql") as f:
    con.execute(f.read())

print(f"Database created with all tables at {db_path}")
con.close()