![smoke-test](https://github.com/crabbyPatty10/Single-Cell-Reanalysis/actions/workflows/smoke.yml/badge.svg)

# Single-Cell Reanalysis

## Smoke test (end-to-end)
This repo includes a minimal, end-to-end smoke test to confirm the pipeline runs start -> finish.

### 1) Create environment (Windows / conda)
    conda create -n scr_smoke -c conda-forge python=3.11 scanpy anndata python-igraph leidenalg pyyaml matplotlib pandas numpy -y
    conda activate scr_smoke

### 2) Generate toy data
    python src/pipeline/make_toy_data.py

### 3) Run pipeline
    python src/pipeline/run.py --config configs/smoke.yaml

### Outputs
- outputs/processed_smoke.h5ad
- reports/smoke/qc_summary.csv
- reports/smoke/qc_violin.png
- reports/smoke/qc_counts_vs_mito.png
- reports/smoke/umap_leiden.png

## Full analysis run (discovery + validation figures)
From the repo root, with the analysis environment activated:

    conda activate scr_smoke
    python scripts/run_full.py --config configs/full.yaml

This regenerates:
- results/tables/state_summary.tsv
- results/tables/paired_subject_deltas.tsv
- results/figures/fig1*.png
- results/figures/figV*.png (if validation tables exist)

## Data required (not included in repo)

You can run the **smoke test** end-to-end with the toy dataset (no private data).
To reproduce the **full discovery analysis** (and optional validation figures), you must provide the following inputs locally.

### Discovery inputs (required for full run)

Expected paths:

- `data/discovery/taurus_v3/myeloid_final.h5ad`
- `data/discovery/taurus_v3/samples.tsv`
- `src/gene_sets/hallmark_selected.gmt` (already in repo)

`samples.tsv` must be whitespace- or tab-delimited and include columns:

- `sample_id`
- `subject_id`
- `timepoint`
- `response` (optional but recommended)

Note: `scripts/run_full.py` will automatically create:

- `data/discovery/taurus_v3/paired_for_run_discovery.tsv`

from `samples.tsv` if it does not already exist.

### Validation outputs (optional)

If you already generated validation tables for **GSE298464**, place them under:

- `results/tables/validation_gse298464/`

`run_full.py` will generate validation figures **only if** the required TSV tables exist in that folder.

## Environment setup (recommended)

This repo uses a pinned conda-forge environment.

1) Create the environment:

    conda env create -f environment.yml

2) Activate it:

    conda activate scr_smoke

3) Quick import check:

    python -c "import numpy, pandas, scanpy, anndata; print('ok')"


