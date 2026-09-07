CREATE TABLE datasets (
    id INTEGER PRIMARY KEY,
    accession TEXT,
    n_cells_raw INTEGER,
    n_cells_after_qc INTEGER,
    download_date TEXT,
    source_url TEXT
);

CREATE TABLE qc_metrics (
    id INTEGER PRIMARY KEY,
    dataset_id INTEGER REFERENCES datasets(id),
    min_genes_threshold INTEGER,
    max_pct_mt_threshold REAL,
    cells_removed INTEGER,
    pct_removed REAL
);

CREATE TABLE clusters (
    id INTEGER PRIMARY KEY,
    dataset_id INTEGER REFERENCES datasets(id),
    leiden_label TEXT,
    n_cells INTEGER,
    dominant_cell_type TEXT,
    dominant_cell_type_pct REAL
);

CREATE TABLE gene_expression_summary (
    id INTEGER PRIMARY KEY,
    cluster_id INTEGER REFERENCES clusters(id),
    gene_symbol TEXT,
    ensembl_id TEXT,
    pct_cells_expressing REAL,
    mean_expression REAL
);

CREATE TABLE runs (
    id INTEGER PRIMARY KEY,
    timestamp TEXT,
    git_commit TEXT,
    scanpy_version TEXT,
    n_pcs_used INTEGER,
    leiden_resolution REAL
);