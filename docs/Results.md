# Results

This page summarizes the key discovery results and the external validation.

> Note: Figures shown here are copied into `docs/assets/figures/` for easy viewing in GitHub.

---

## Discovery (primary)

### Figure 1 — Paired deltas by clinical response (Remission vs Non_Remission)

![Fig 1: paired deltas by response](assets/figures/fig1_response_deltas.png)

**Takeaways (discovery):**
- Subjects who reached **Remission** show larger decreases (post − baseline) in **ROS**, **OXPHOS**, and **inflammatory response** signals compared with **Non_Remission**.
- The **MITO/ROS-high fraction** also decreases more in Remission on average.
- This supports an “immune-state normalization” interpretation in the Remission group.

---

### Figure S1 — MITO/ROS-high change by response (single-metric view)

![Fig S1: mito/ROS delta](assets/figures/fig1_mito_ros_delta.png)

**Takeaways (discovery):**
- Focused view of the key state: the **MITO/ROS-high fraction** decreases more in **Remission** than **Non_Remission** on average.
- The spread/outliers highlight patient heterogeneity, motivating the multi-metric view in Figure 1.

---

### Figure 1b — Relationship to inflammation

![Fig 1b: mito/ROS vs inflammation](assets/figures/fig1b_mito_ros_vs_inflammation.png)

**Takeaways (discovery):**
- The mito/ROS signal is compared against an inflammation-related signal to contextualize biology.
- The plot helps separate “general inflammation” from a more specific mitochondrial/oxidative pattern.

---

### Figure 1c–1d — Mean delta summaries

![Fig 1c: mean ROS delta](assets/figures/fig1c_delta_mean_ROS.png)

![Fig 1d: mean OXPHOS delta](assets/figures/fig1d_delta_mean_OXPHOS.png)

**Takeaways (discovery):**
- These panels provide compact summaries of mean changes in ROS- and OXPHOS-related measures.
- Together they support the primary mito/ROS change pattern.

---

## Validation (external)

### Figure V1 — Validation cohort group comparison

![Fig V1: validation group comparison](assets/figures/figV1_gse298464_mito_ros_by_group.png)

**Takeaways (validation):**
- The key mito/ROS pattern is evaluated in an independent dataset (GSE298464).
- This is a first-pass validation intended to check directional consistency.

---

### Figure V2 — Validation paired deltas

> Note: In GSE298464 validation, response labels for the paired-delta plot (Fig V2) are inferred (the paired-deltas table lacks an explicit `response` column); see `scripts/make_fig1_validation_gse298464.py` for the mapping logic.
> See [`scripts/make_fig1_validation_gse298464.py`](../scripts/make_fig1_validation_gse298464.py) for the mapping logic.

![Fig V2: validation deltas](assets/figures/figV2_gse298464_mito_ros_delta.png)

**Takeaways (validation):**
- The same paired-delta logic is applied to the validation cohort where pairing is available.
- Agreement in direction strengthens confidence that the discovery signal is not dataset-specific.

---

## Where to find the underlying tables

Discovery:
- `results/tables/state_summary.tsv`
- `results/tables/paired_subject_deltas.tsv`
- `results/tables/paired_subject_deltas_with_metadata.tsv`

Validation (if generated):
- `results/tables/validation_gse298464/`
