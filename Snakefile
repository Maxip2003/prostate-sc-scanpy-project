import yaml

# Which params file to use — override at the command line with:
# snakemake --cores 1 --config config_path=config/params.test.yaml
CONFIG_PATH = config.get("config_path", "config/params.yaml")

with open(CONFIG_PATH) as f:
    params = yaml.safe_load(f)

OUTPUT_DIR = params["dataset"]["output_dir"]
DB_PATH = params["dataset"]["db_path"]

rule all:
    input:
        f"{OUTPUT_DIR}/.pipeline_complete"

rule create_db:
    output:
        DB_PATH
    shell:
        "python scripts/create_db.py --config {CONFIG_PATH}"

rule load:
    input:
        CONFIG_PATH
    output:
        f"{OUTPUT_DIR}/01_loaded.h5ad"
    shell:
        "python scripts/01_load.py --config {CONFIG_PATH}"

rule clean:
    input:
        f"{OUTPUT_DIR}/01_loaded.h5ad"
    output:
        f"{OUTPUT_DIR}/02_cleaned.h5ad"
    shell:
        "python scripts/02_clean.py --config {CONFIG_PATH}"

rule cluster:
    input:
        f"{OUTPUT_DIR}/02_cleaned.h5ad"
    output:
        f"{OUTPUT_DIR}/03_clustered.h5ad"
    shell:
        "python scripts/03_cluster.py --config {CONFIG_PATH}"

rule annotate:
    input:
        f"{OUTPUT_DIR}/03_clustered.h5ad"
    output:
        f"{OUTPUT_DIR}/04_annotated.h5ad"
    shell:
        "python scripts/04_annotate.py --config {CONFIG_PATH}"

rule load_to_db:
    input:
        h5ad = f"{OUTPUT_DIR}/04_annotated.h5ad",
        db = DB_PATH
    output:
        touch(f"{OUTPUT_DIR}/.pipeline_complete")
    shell:
        "python scripts/05_load_to_db.py --config {CONFIG_PATH}"