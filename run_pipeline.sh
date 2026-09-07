#!/bin/bash
set -e #stop if any step fails

echo "Running pipeline..."
python scripts/01_load.py
python scripts/02_clean.py
python scripts/03_cluster.py
python scripts/04_annotate.py
python scripts/05_load_to_db.py
echo "Pipeline finished."
