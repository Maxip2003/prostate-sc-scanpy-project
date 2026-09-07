#!/bin/bash
set -e

echo "Running load + clean steps..."
python scripts/01_load.py
python scripts/02_clean.py

echo ""
echo "Check figures/violin_qc_before.png and figures/violin_qc_2_before.png"
echo "Then set the right thresholds in config/params.yaml and run ./run_pipeline.sh"
