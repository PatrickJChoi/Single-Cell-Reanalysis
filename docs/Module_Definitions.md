# Module Definitions

This project defines “cell states” using portable gene-set (module) scores computed in myeloid cells.

## Source
- Gene sets: MSigDB Hallmark (gene symbols)
- We will use these exact set names (6 Hallmark modules) for scoring:
  - HALLMARK_OXIDATIVE_PHOSPHORYLATION
  - HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
  - HALLMARK_TNFA_SIGNALING_VIA_NFKB
  - HALLMARK_INFLAMMATORY_RESPONSE
  - HALLMARK_INTERFERON_ALPHA_RESPONSE
  - HALLMARK_INTERFERON_GAMMA_RESPONSE

## Primary modules (focus)
- HALLMARK_OXIDATIVE_PHOSPHORYLATION (OXPHOS)
- HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY (ROS)

## Comparator biology modules
- HALLMARK_TNFA_SIGNALING_VIA_NFKB (TNF/NFκB)
- HALLMARK_INFLAMMATORY_RESPONSE (Inflammation)
- HALLMARK_INTERFERON_ALPHA_RESPONSE (IFN-α)
- HALLMARK_INTERFERON_GAMMA_RESPONSE (IFN-γ)

## Control scores (not Hallmark)
### Ribosomal score
- Definition: genes with prefixes RPL* and RPS* (gene symbols)

### Cell cycle scores
- S score + G2M score
- Method: Scanpy cell cycle scoring using canonical S and G2M gene lists

## State-high definition (portable rule)
Within myeloid cells **per sample**:
- For each module score, label cells as **state-high** if their score is in the **top 20%** within that sample.
- Equivalent quantile threshold: q = 0.80 within each (sample_id × myeloid) group.

This state-high rule must be applied identically in discovery and validation cohorts.

