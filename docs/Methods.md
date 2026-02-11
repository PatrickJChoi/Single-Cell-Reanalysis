\# Methods (high-level)



This document describes the analysis flow at a level suitable for readers who want to understand \*what was done\* and \*why\*, without diving into implementation details.



---



\## Overview



1\) Load preprocessed discovery AnnData (`.h5ad`) and sample metadata  

2\) Compute QC metrics and standardize expression (where applicable)  

3\) Embed and cluster myeloid cells (PCA → neighbors → UMAP → Leiden)  

4\) Score biologically meaningful programs (e.g., mitochondrial/ROS, OXPHOS, inflammation)  

5\) Compute \*\*paired within-subject deltas\*\* (post − baseline)  

6\) Summarize group-level effects and generate figures  

7\) Validate patterns in an external cohort (GSE298464) when available



---



\## Inputs



\### Discovery (required for full run)



\- `data/discovery/taurus\_v3/myeloid\_final.h5ad`

\- `data/discovery/taurus\_v3/samples.tsv`

\- Gene sets: `src/gene\_sets/hallmark\_selected.gmt`



\### Validation (optional)



Validation uses tables/outputs under:

\- `results/tables/validation\_gse298464/`



Raw validation construction is supported by scripts but not run automatically.



---



\## Key analysis choices



\### Clustering and visualization (smoke test and/or exploratory steps)



\- Highly variable genes (HVGs) selected for dimensionality reduction

\- PCA used for embedding

\- Neighborhood graph constructed from PCA space

\- UMAP computed for visualization

\- Leiden clustering used for community detection



\### Biological scoring



\- Gene-set–based scores are computed to represent biological programs (e.g., ROS/mitochondrial stress, OXPHOS, inflammation).

\- Scores are used to interpret clusters/states and to quantify changes between timepoints.



\### Paired delta design (core result)



For each subject:

\- Identify baseline and post samples

\- Compute delta per metric:  

&nbsp; \*\*delta = post − baseline\*\*

\- Aggregate deltas at the subject level and compare across groups (e.g., response categories)



This design emphasizes \*\*within-person change\*\*, reducing confounding from baseline differences between individuals.



---



\## Outputs



Discovery tables:

\- `results/tables/state\_summary.tsv`

\- `results/tables/paired\_subject\_deltas.tsv`



Figures:

\- `results/figures/fig1\*.png`

\- Validation (if present): `results/figures/figV\*.png`



---



\## Reproducibility



\- Pinned environment: `environment.yml`

\- Smoke test: `configs/smoke.yaml` (toy data, end-to-end)

\- Full run: `configs/full.yaml` + `scripts/run\_full.py`

\- CI: GitHub Actions runs the smoke test on pushes/PRs



