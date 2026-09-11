
## Quickstart

```bash
# Environment
conda create -n prostata python=3.11 -y
conda activate prostata
pip install -r requirements.txt

# Run the full pipeline on the small test dataset (already included in the repo)
python scripts/create_db.py --config config/params.test.yaml
python scripts/01_load.py --config config/params.test.yaml
python scripts/02_clean.py --config config/params.test.yaml
python scripts/03_cluster.py --config config/params.test.yaml
python scripts/04_annotate.py --config config/params.test.yaml
python scripts/05_load_to_db.py --config config/params.test.yaml

# Query the results
python scripts/06_query.py --gene PCA3 --config config/params.test.yaml
```

To reproduce the full result on the 68,322-cell dataset: download it from [CELLxGENE](https://cellxgene.cziscience.com/e/881e0e6b-a185-45f1-9b56-81606e91e7ff.cxg/) (accession `68b23fda-7191-46a5-8870-819feca3e66e`), place it at the path set in `config/params.yaml`, and run the same six commands without `--config` (it defaults to `config/params.yaml`).

## Database schema

Five tables in DuckDB: `datasets`, `qc_metrics`, `clusters`, `gene_expression_summary`, and `runs` — the last one logs the exact git commit, package versions, and parameters used for each run, so any result can be traced back to the code that produced it. Foreign keys enforce referential integrity between tables.

## Reproducibility

Every parameter (QC thresholds, number of PCA components, positive-control gene, dataset paths) lives in a single YAML config, never hardcoded in the scripts — switching datasets or re-running with different thresholds never requires touching code. `tests/test_pipeline.py` runs the entire pipeline end-to-end against the small test dataset and verifies the database is populated correctly; this runs automatically on every push via GitHub Actions.

## Limitations

- This is an exploratory analysis, not a publication-ready finding. The PCA3–malignancy association is correlational, observed within a single dataset, and has not been statistically tested or corrected for patient-of-origin (cells from the same donor are not independent samples, and some clusters may be dominated by a small number of patients out of the cohort's 24).
- The doublet-suspicious upper threshold on genes-per-cell was chosen visually/by percentile, not with a dedicated doublet-detection tool.
- Leiden clustering resolution and the number of PCA components were chosen from the elbow plot and are not exhaustively optimised.

## Data source

Apostolov, E., Roden, D. L., Holliday, H., Cazet, A., Harvey, K., Zhang, H., Wu, S. Z., van der Leij, S., Jieun Kim, H., Selth, L. A., Bartonicek, N., Al-Eryani, G., Reeves, J. L., He, M., Lundeberg, J., Potter, A. J., Kench, J. G., Stricker, P. D., Joshua, A. M., Horvath, L. G., … Swarbrick, A. (2026). *Single-Cell and Spatial Transcriptomic Profiling Reveals Epithelial Functional States and Fibroblast Phenotypes in Hormone Therapy-Naïve Localized Prostate Cancer*. Cancer Research, 86(8), 1836–1853. https://doi.org/10.1158/0008-5472.CAN-25-1202

Data accessed via [CELLxGENE Discover](https://cellxgene.cziscience.com/collections/bdac7a53-fe34-4f04-8c46-f9bd5297c099).

## Author

Pedro Carrillo Alarcón — MSc Bioinformatics — [github.com/Maxip2003](https://github.com/Maxip2003)