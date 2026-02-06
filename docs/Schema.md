# Schema

This document defines the required columns for the two key tables:
1) `samples.tsv` (one row per sample/biopsy)
2) `state_summary.tsv` (one row per sample × celltype; initially celltype = myeloid)

All TSVs are tab-separated with a header row.

---

## 1) samples.tsv

### Required columns
- sample_id
  - Unique ID for a biopsy/sample used in `adata.obs` to group cells
- subject_id
  - Patient/donor identifier
- cohort
  - Discovery vs validation label (e.g., "TAURUS_GSE282122" or "GSE298464")
- timepoint
  - Pre/post (or discrete timepoint label if multiple post timepoints exist)
- treatment
  - anti-TNF drug name if available (e.g., adalimumab, infliximab)
- response
  - Binary or categorical response label (e.g., responder/nonresponder OR remission/nonremission)
- disease
  - CD/UC/IBD if available
- tissue_site
  - Colon/ileum/rectum/etc if available
- batch
  - Library/batch identifier if available

### Optional (include if available)
- collection_day
  - numeric day relative to treatment start
- dataset_accession
  - GEO accession or source identifier
- notes
  - free text for quirks (missing label, ambiguous timepoint, etc.)

---

## 2) state_summary.tsv

One row per (sample_id × celltype). For this project, we will focus on celltype = "myeloid".

### Required identifier columns
- sample_id
- subject_id
- cohort
- timepoint
- response
- celltype
  - expected values: "myeloid" (future-proof if we add other cell types later)
- n_cells
  - number of cells included in this row (e.g., number of myeloid cells in sample)

### Module score columns (means)
Mean module score across cells in that (sample_id × celltype):
- mean_HALLMARK_OXIDATIVE_PHOSPHORYLATION
- mean_HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
- mean_HALLMARK_TNFA_SIGNALING_VIA_NFKB
- mean_HALLMARK_INFLAMMATORY_RESPONSE
- mean_HALLMARK_INTERFERON_ALPHA_RESPONSE
- mean_HALLMARK_INTERFERON_GAMMA_RESPONSE
- mean_ribosomal_score
- mean_S_score
- mean_G2M_score

### State-high fraction columns
Fraction of cells in the top 20% of the module score distribution within-sample (q=0.80), computed within myeloid cells per sample:
- frac_statehigh_HALLMARK_OXIDATIVE_PHOSPHORYLATION
- frac_statehigh_HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
- frac_statehigh_HALLMARK_TNFA_SIGNALING_VIA_NFKB
- frac_statehigh_HALLMARK_INFLAMMATORY_RESPONSE
- frac_statehigh_HALLMARK_INTERFERON_ALPHA_RESPONSE
- frac_statehigh_HALLMARK_INTERFERON_GAMMA_RESPONSE

### Optional QC/summary columns (include if available)
- pct_mito_mean
  - mean percent mitochondrial reads per cell (sample × celltype)
- total_counts_mean
  - mean total UMI counts per cell
- n_genes_by_counts_mean
  - mean number of genes detected per cell
- batch
  - include if the same sample_id appears across batches (rare)

---

## Notes / conventions
- Use consistent string values across tables (e.g., response labels).
- `state_summary.tsv` should be generated automatically from the AnnData object + `samples.tsv`.

