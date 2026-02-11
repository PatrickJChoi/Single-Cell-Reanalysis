# Data acquisition (not included in repo)

This repository is designed to be reproducible without distributing large or restricted datasets.
You can always run the **smoke test** using the toy dataset.

For the **full discovery + optional validation figures**, provide the inputs below in the expected paths.

---

## 1) Discovery dataset (TAURUS / local discovery)

### Required files

Place these files here:

- `data/discovery/taurus_v3/myeloid_final.h5ad`
- `data/discovery/taurus_v3/samples.tsv`

The gene set file used by the discovery pipeline is already in the repo:

- `src/gene_sets/hallmark_selected.gmt`

### `samples.tsv` format

`samples.tsv` must be whitespace- or tab-delimited and include:

- `sample_id`
- `subject_id`
- `timepoint`
- `response` (optional but recommended)

### Paired file creation

The full runner will generate the paired table automatically (if missing):

- `data/discovery/taurus_v3/paired_for_run_discovery.tsv`

Generated columns:

- `sample_id`, `subject_id`, `timepoint`, `response`

---

## 2) Validation dataset (GSE298464)

This repo includes scripts to build and analyze a validation dataset from **GSE298464**.

By default, `scripts/run_full.py` **does not** download or rebuild validation data automatically.
It only creates validation figures if precomputed validation tables exist.

### Expected validation tables

If you have already generated the validation tables, place them under:

- `results/tables/validation_gse298464/`

When these TSVs exist, `scripts/run_full.py` will generate:

- `results/figures/figV1_gse298464_mito_ros_by_group.png`
- `results/figures/figV2_gse298464_mito_ros_delta.png`

---

## 3) What you can run without any real data

### Smoke test (toy dataset)

1) Create toy data:

    conda activate scr_smoke  
    python src/pipeline/make_toy_data.py

2) Run the smoke pipeline:

    python src/pipeline/run.py --config configs/smoke.yaml

Outputs:

- `outputs/processed_smoke.h5ad`
- `reports/smoke/` (QC + UMAP plots)

---

## 4) Full run (requires discovery inputs)

    conda activate scr_smoke
    python scripts/run_full.py --config configs/full.yaml

Outputs:

- `results/tables/state_summary.tsv`
- `results/tables/paired_subject_deltas.tsv`
- `results/figures/fig1*.png`

Optional (if validation tables exist):

- `results/figures/figV1*.png`
- `results/figures/figV2*.png`
