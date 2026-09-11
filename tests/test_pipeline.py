import subprocess
from pathlib import Path
import duckdb

TEST_CONFIG = "config/params.test.yaml"
TEST_DB = Path("databases/test_prostate.duckdb")

def run_step(script):
    cmd = ["python", f"scripts/{script}", "--config", TEST_CONFIG]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"{script} failed:\n{result.stderr}"

def test_pipeline_runs_end_to_end():
    # Start from a clean slate — remove any leftover test database from a previous run
    TEST_DB.unlink(missing_ok=True)

    run_step("create_db.py")

    for script in ["01_load.py", "02_clean.py", "03_cluster.py", "04_annotate.py", "05_load_to_db.py"]:
        run_step(script)

    con = duckdb.connect(str(TEST_DB))
    n_datasets = con.execute("SELECT COUNT(*) FROM datasets").fetchone()[0]
    n_clusters = con.execute("SELECT COUNT(*) FROM clusters").fetchone()[0]
    n_expr = con.execute("SELECT COUNT(*) FROM gene_expression_summary").fetchone()[0]
    con.close()

    assert n_datasets == 1
    assert n_clusters > 0
    assert n_expr > 0