![smoke-test](https://github.com/crabbyPatty10/Single-Cell-Reanalysis/actions/workflows/smoke.yml/badge.svg)

# Single-Cell Reanalysis

## Abstract

This repo contains a reproducible single-cell reanalysis pipeline that quantifies **within-subject (post − baseline)** changes in myeloid immune programs and links them to clinical outcome. In the discovery cohort (n=34), subjects reaching **Remission** show larger paired decreases in **ROS**, **OXPHOS**, and **inflammatory response** signals than **Non_Remission**, consistent with immune-state normalization. The project includes an end-to-end smoke test (no private data), a one-command full runner, and documented Results/Methods with embedded figures.

- Results: `docs/Results.md`
- Methods: `docs/Methods.md`
- Data acquisition: `docs/Data_Acquisition.md`

## Quickstart

Smoke test (no real data required):

    conda env create -f environment.yml
    conda activate scr_smoke
    python src/pipeline/make_toy_data.py
    python src/pipeline/run.py --config configs/smoke.yaml

Full run (requires discovery inputs; see `docs/Data_Acquisition.md`):

    conda activate scr_smoke
    python scripts/run_full.py --config configs/full.yaml

---

## Docker (recommended for portability)

Build and run the smoke test with no local setup required:

    docker build -t scr .
    docker run --rm scr

To run the full analysis, mount your local data directory into the container:

    docker run --rm \
      -v $(pwd)/data:/app/data \
      -v $(pwd)/results:/app/results \
      scr python scripts/run_full.py --config configs/full.yaml

To run tests:

    docker run --rm scr python -m pytest tests/

To open an interactive shell inside the container:

    docker run --rm -it scr bash

---

## Environment setup (recommended)

This repo uses a pinned conda-forge environment.

1) Create the environment:

    conda env create -f environment.yml

2) Activate it:

    conda activate scr_smoke

3) Quick import check:

    python -c "import numpy, pandas, scanpy, anndata; print('ok')"

### Alternative: manual environment creation

If you prefer not to use `environment.yml`, you can create the environment manually:

    conda create -n scr_smoke -c conda-forge python=3.11 scanpy anndata python-igraph leidenalg pyyaml matplotlib pandas numpy -y
    conda activate scr_smoke

---

## Smoke test (end-to-end)

This repo includes a minimal, end-to-end smoke test to confirm the pipeline runs start -> finish.

1) Generate toy data:

    python src/pipeline/make_toy_data.py

2) Run pipeline:

    python src/pipeline/run.py --config configs/smoke.yaml

Outputs (not committed to git):
- `outputs/processed_smoke.h5ad`
- `reports/smoke/qc_summary.csv`
- `reports/smoke/qc_violin.png`
- `reports/smoke/qc_counts_vs_mito.png`
- `reports/smoke/umap_leiden.png`

---

## Full analysis run (discovery + validation figures)

From the repo root, with the analysis environment activated:

    conda activate scr_smoke
    python scripts/run_full.py --config configs/full.yaml

This regenerates:
- `results/tables/state_summary.tsv`
- `results/tables/paired_subject_deltas.tsv`
- `results/figures/fig1*.png`
- `results/figures/figV*.png` (if validation tables exist)

---

## Data required (not included in repo)

See `docs/Data_Acquisition.md` for detailed setup instructions.

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
