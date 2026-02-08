from __future__ import annotations

import sys
from pathlib import Path

# Ensure repo root is on PYTHONPATH so `import src...` works
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import argparse
import pandas as pd
import scanpy as sc

from src.state_scoring import load_gmt, score_modules, assign_state_high

HALLMARKS = [
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]

DEFAULT_META_COLS = [
    "sample_id",        # biopsy/sample id (CID...)
    "Patient",          # subject label (UC1/CD10...)
    "Disease",
    "Site",
    "Treatment",        # Pre/Post
    "Remission_status", # Remission / Non_Remission
    "Inflammation",     # Inflamed / Non_Inflamed
    "Age",
    "Gender",
    "Ethnicity",
]

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Discovery pipeline: score myeloid states and write summary TSVs."
    )
    ap.add_argument("--h5ad", required=True, help="Path to processed myeloid .h5ad")
    ap.add_argument("--paired", required=False, help="Path to paired_sample_list.csv (tab-separated) for TAURUS")
    ap.add_argument("--gmt", required=True, help="Path to hallmark_selected.gmt")
    ap.add_argument("--outdir", default="results/tables", help="Output directory for TSV tables")
    return ap.parse_args()

def safe_read_paired(paired_path: Path) -> pd.DataFrame:
    """
    TAURUS paired_sample_list.csv is tab-separated with columns: sample_id, Category
    """
    df = pd.read_csv(paired_path, sep="\t")
    required = {"sample_id", "Category"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"paired file missing required columns: {missing}. Found: {list(df.columns)}")
    return df

def build_state_summary(adata) -> pd.DataFrame:
    """
    Build one row per biopsy/sample_id:
    - metadata
    - n_myeloid_cells
    - mean scores for each module
    - pct_high_MITO_ROS
    """
    # Metadata (one row per sample)
    meta_cols = [c for c in DEFAULT_META_COLS if c in adata.obs.columns]
    if "sample_id" not in meta_cols:
        raise ValueError("adata.obs must contain 'sample_id' column for biopsy grouping.")

    samples = (
        adata.obs[meta_cols]
        .drop_duplicates()
        .rename(columns={"Patient": "subject_id", "Treatment": "timepoint", "Remission_status": "response"})
        .copy()
    )

    # Cell counts per biopsy
    cell_counts = adata.obs.groupby("sample_id").size().rename("n_myeloid_cells").reset_index()
    samples = samples.merge(cell_counts, on="sample_id", how="left")

    # Score columns that were added by score_modules
    score_cols = [c for c in adata.obs.columns if c.startswith("score_HALLMARK_")]
    if len(score_cols) == 0:
        raise ValueError("No score_HALLMARK_* columns found. Did score_modules run correctly?")

    # Aggregate to sample-level
    agg_dict = {f"mean_{c}": (c, "mean") for c in score_cols}
    agg_dict["pct_high_MITO_ROS"] = ("state_high_MITO_ROS", lambda x: float(x.mean() * 100.0))

    agg = adata.obs.groupby("sample_id").agg(**agg_dict).reset_index()

    state_summary = samples.merge(agg, on="sample_id", how="left")
    return state_summary

def main() -> None:
    args = parse_args()

    h5ad_path = Path(args.h5ad)
    gmt_path = Path(args.gmt)
    paired_path = Path(args.paired) if args.paired else None

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load
    adata = sc.read_h5ad(h5ad_path)

    # Score hallmark modules (in-place)
    gene_sets = load_gmt(gmt_path)
    used_sets = {m: gene_sets[m] for m in HALLMARKS if m in gene_sets}
    if len(used_sets) != len(HALLMARKS):
        missing = [m for m in HALLMARKS if m not in used_sets]
        raise ValueError(f"Missing hallmark sets in GMT: {missing}")

    score_modules(adata, used_sets)

    # Create state-high columns (top 20% per sample_id => q=0.8)
    score_cols = [f"score_{m}" for m in HALLMARKS]
    for sc_col in score_cols:
        if sc_col not in adata.obs.columns:
            raise ValueError(f"Expected score column not found in adata.obs: {sc_col}")
        out_col = "state_high_" + sc_col.replace("score_", "")
        assign_state_high(adata, score_col=sc_col, groupby="sample_id", q=0.8, out_col=out_col)

    # Combined state: OXPHOS_high AND ROS_high
    adata.obs["state_high_MITO_ROS"] = (
        adata.obs["state_high_HALLMARK_OXIDATIVE_PHOSPHORYLATION"]
        & adata.obs["state_high_HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]
    )

    # Write biopsy-level summary
    state_summary = build_state_summary(adata)
    state_summary_path = outdir / "state_summary.tsv"
    state_summary.to_csv(state_summary_path, sep="\t", index=False)
    print(f"Wrote: {state_summary_path}")

    # Optional: paired subject deltas table
    if paired_path and paired_path.exists():
        paired = safe_read_paired(paired_path)
        paired_state = state_summary.merge(paired, on="sample_id", how="inner").copy()

        # response group from Category
        paired_state["resp_group"] = paired_state["Category"].str.contains("Non_Remission").map(
            {True: "Non_Remission", False: "Remission"}
        )

        # Subject/timepoint averages then Post-Pre deltas
        # (Use only numeric columns to avoid categorical subtraction errors)
        subj_tp = paired_state.groupby(["subject_id", "resp_group", "timepoint"], as_index=False).agg(
            pct_high_MITO_ROS=("pct_high_MITO_ROS", "mean"),
            mean_OXPHOS=("mean_score_HALLMARK_OXIDATIVE_PHOSPHORYLATION", "mean"),
            mean_ROS=("mean_score_HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "mean"),
            mean_TNFA=("mean_score_HALLMARK_TNFA_SIGNALING_VIA_NFKB", "mean"),
            mean_INFLAM=("mean_score_HALLMARK_INFLAMMATORY_RESPONSE", "mean"),
            frac_inflamed=("Inflammation", lambda x: float((x == "Inflamed").mean()) if len(x) else 0.0),
            n_biopsies=("sample_id", "nunique"),
        )

        pre = subj_tp[subj_tp["timepoint"] == "Pre"].set_index(["subject_id", "resp_group"])
        post = subj_tp[subj_tp["timepoint"] == "Post"].set_index(["subject_id", "resp_group"])
        idx = pre.index.intersection(post.index)

        numeric_cols = [
            "pct_high_MITO_ROS",
            "mean_OXPHOS",
            "mean_ROS",
            "mean_TNFA",
            "mean_INFLAM",
            "frac_inflamed",
            "n_biopsies",
        ]

        deltas = (post.loc[idx, numeric_cols] - pre.loc[idx, numeric_cols]).reset_index()
        deltas = deltas.rename(columns={
            "pct_high_MITO_ROS": "delta_pct_high_MITO_ROS",
            "mean_OXPHOS": "delta_mean_OXPHOS",
            "mean_ROS": "delta_mean_ROS",
            "mean_TNFA": "delta_mean_TNFA",
            "mean_INFLAM": "delta_mean_INFLAM",
            "frac_inflamed": "delta_frac_inflamed",
            "n_biopsies": "delta_n_biopsies",
        })

        deltas_path = outdir / "paired_subject_deltas.tsv"
        deltas.to_csv(deltas_path, sep="\t", index=False)
        print(f"Wrote: {deltas_path}")

if __name__ == "__main__":
    main()
