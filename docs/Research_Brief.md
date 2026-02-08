# Research Brief

Personal intellectual question:
"Why do some IBD patients fail anti-TNF therapy even when inflammation decreases--could a persistent mitochondrial/oxidative-stress state in gut myeloid cells predict nonresponse across cohorts?"

Formal research question:
"In longitudinal gut biopsy scRNA-seq from Crohn's/UC patients treated with adalimumab, does donor-level abundance of a predefined myeloid mitochondrial/oxidative-stress state persist in nonresponders after treatment, independent of IFN/TNF states, and replicate in an external anti-TNF cohort?"

---

## 2026-02-08 -- TAURUS discovery run + Figure 1

### Core question
Does a portable myeloid mitochondrial/oxidative-stress state persist post-anti-TNF in IBD nonresponders (beyond inflammatory changes) and replicate across cohorts?

### Discovery dataset
TAURUS (processed myeloid-only AnnData: `myeloid_final.h5ad`) with paired sample list for longitudinal subset.

### State definition
- Module scoring: Hallmark OXPHOS, ROS, TNFA/NFKB, Inflammatory response, IFN-alpha, IFN-gamma.
- State-high rule: top 20% within each biopsy sample (`sample_id`) among myeloid cells.
- Combined state: MITO_ROS-high = (OXPHOS_high AND ROS_high).

### Outputs (generated locally; not committed)
- Tables:
  - `results/tables/state_summary.tsv` (per biopsy sample)
  - `results/tables/paired_subject_deltas.tsv` (paired subject Post-Pre deltas)
- Figures:
  - `results/figures/fig1_mito_ros_delta.png` (delta MITO_ROS by response)
  - `results/figures/fig1b_mito_ros_vs_inflammation.png` (delta MITO_ROS vs delta fraction inflamed)

### Key numbers (TAURUS paired subjects)
- Paired subjects: 32 (Non_Remission n=18, Remission n=14)
- Median delta(Post-Pre) percent MITO_ROS-high:
  - Non_Remission: -1.15
  - Remission: -2.51

### Beyond inflammation quick check
- Correlation between delta MITO_ROS and delta fraction inflamed is small (~0.12).
- Simple regression: delta MITO_ROS ~ nonremission + delta fraction inflamed shows delta fraction inflamed is not strongly associated in this small sample.
- Interpretation: preliminary support that MITO_ROS change is not explained solely by inflammation change, but sample size is small.

### Next steps
1) Decide final primary endpoint (MITO_ROS vs ROS vs OXPHOS) and finalize Figure 1 style.
2) Add robustness checks (disease/site stratification; continuous scores vs percent-high).
3) Select validation cohort and build validation runner to output the same TSV schema.

### Interpretation (quick)
In the TAURUS discovery cohort (32 paired subjects), the MITO_ROS-high fraction generally decreased from pre to post in both responders and nonresponders. Changes in this state did not track strongly with changes in the fraction of inflamed biopsies, suggesting the MITO/ROS program is not simply a proxy for inflammation status in this small sample. Mean ROS and mean OXPHOS module scores also shifted downward post-treatment, with larger decreases in the remission group.
