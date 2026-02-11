#!/usr/bin/env python3
"""
run_discovery.py

Discovery/validation runner:
- Load myeloid AnnData (.h5ad)
- Score Hallmark modules (adds score_* columns)
- Define state_high per sample (top 20% => q=0.8)
- Build per-sample table: results/tables/.../state_summary.tsv
- OPTIONAL: if --paired is provided, compute subject-level Post-Pre deltas:
  results/tables/.../paired_subject_deltas.tsv

This script is designed to be robust across:
- TAURUS discovery (many biopsies per subject/timepoint)
- GSE298464 validation (1 biopsy per subject/timepoint)
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

# Ensure repo root is on sys.path so `import src...` works when running `python scripts/...`
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.state_scoring import load_gmt, score_modules, assign_state_high  # noqa: E402


HALLMARKS = [
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]

OXPHOS = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
ROS = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"


def _read_table_anysep(path: str) -> pd.DataFrame:
    """Read TSV/CSV without the user caring which separator it is."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Could not find: {path}")
    # Try tab first (your files are often .csv but actually tab-separated)
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path)  # comma fallback


def _norm_timepoint(x: str) -> str:
    s = str(x).strip().lower()
    if s.startswith("pre"):
        return "Pre"
    if s.startswith("post"):
        return "Post"
    # if unknown, keep original but title-case
    return str(x).strip().title()


def _ensure_score_cols(adata: sc.AnnData, hallmark_names: List[str]) -> None:
    missing = [f"score_{h}" for h in hallmark_names if f"score_{h}" not in adata.obs.columns]
    if missing:
        raise KeyError(
            "Missing expected module score columns:\n"
            + "\n".join(missing)
            + "\n\nThis usually means module scoring didn't run or used a different prefix."
        )


def _ensure_state_high_cols(adata: sc.AnnData, hallmark_names: List[str], groupby: str = "sample_id", q: float = 0.8) -> None:
    """
    Force-create state_high_<HALLMARK> boolean columns using explicit out_col names.
    Also supports alternate naming if an older helper created state_high_score_*.
    """
    for h in hallmark_names:
        expected = f"state_high_{h}"
        score_col = f"score_{h}"
        alt = f"state_high_score_{h}"

        if expected in adata.obs.columns:
            continue

        # If alt exists, copy it over
        if alt in adata.obs.columns:
            adata.obs[expected] = adata.obs[alt].astype(bool)
            continue

        # Otherwise compute it explicitly
        if score_col not in adata.obs.columns:
            raise KeyError(f"Cannot create {expected} because {score_col} is missing.")
        assign_state_high(
            adata,
            score_col=score_col,
            groupby=groupby,
            q=q,
            out_col=expected,
        )


def _build_state_summary(adata: sc.AnnData, outdir: Path) -> pd.DataFrame:
    """
    One row per sample_id. Includes:
    - metadata columns if present
    - n_myeloid_cells
    - mean_score_<HALLMARK> (per sample)
    - pct_high_<HALLMARK> (per sample; percent)
    - pct_high_MITO_ROS (combined OXPHOS_high & ROS_high; percent)
    - frac_inflamed (if Inflammation column exists; fraction)
    """
    if "sample_id" not in adata.obs.columns:
        raise KeyError("adata.obs is missing required column: sample_id")

    group = adata.obs.groupby("sample_id", observed=True)

    # ---- metadata: keep first value per sample (safe for repeated identical metadata)
    meta_candidates = [
        "sample_id",
        "subject_id",
        "Disease",
        "disease",
        "Site",
        "site",
        "timepoint",
        "response",
        "Inflammation",
        "Age",
        "Gender",
        "Ethnicity",
        "treatment",
        "sex",
        "age",
        "gsm",
        "sample_group",
    ]
    meta_cols = [c for c in meta_candidates if c in adata.obs.columns and c != "sample_id"]
    meta = group[meta_cols].first().reset_index() if meta_cols else group.size().rename("n").reset_index()[["sample_id"]]

    # ---- counts
    n_cells = group.size().rename("n_myeloid_cells").reset_index()

    # ---- means of scores
    score_cols = [f"score_{h}" for h in HALLMARKS if f"score_{h}" in adata.obs.columns]
    means = group[score_cols].mean(numeric_only=True).reset_index()
    means = means.rename(columns={c: f"mean_{c}" for c in score_cols})  # mean_score_...

    # ---- percent high per hallmark
    high_cols = [f"state_high_{h}" for h in HALLMARKS if f"state_high_{h}" in adata.obs.columns]
    pct = group[high_cols].mean(numeric_only=True).reset_index()
    pct = pct.rename(columns={c: f"pct_high_{c.replace('state_high_', '')}" for c in high_cols})
    for c in [c for c in pct.columns if c.startswith("pct_high_")]:
        pct[c] = 100.0 * pct[c]

    # ---- combined MITO_ROS-high
    mito_ros_col = "state_high_MITO_ROS"
    ox_col = f"state_high_{OXPHOS}"
    ros_col = f"state_high_{ROS}"
    if ox_col in adata.obs.columns and ros_col in adata.obs.columns:
        adata.obs[mito_ros_col] = (adata.obs[ox_col].astype(bool) & adata.obs[ros_col].astype(bool))
        mito_pct = group[mito_ros_col].mean(numeric_only=True).rename("pct_high_MITO_ROS").reset_index()
        mito_pct["pct_high_MITO_ROS"] = 100.0 * mito_pct["pct_high_MITO_ROS"]
    else:
        mito_pct = group.size().rename("pct_high_MITO_ROS").reset_index()
        mito_pct["pct_high_MITO_ROS"] = np.nan

    # ---- fraction inflamed if available
    if "Inflammation" in adata.obs.columns:
        infl = adata.obs["Inflammation"].astype(str).str.lower().str.contains("inflamed")
        infl = infl & ~adata.obs["Inflammation"].astype(str).str.lower().str.contains("non")
        adata.obs["_is_inflamed"] = infl.astype(int)
        frac_infl = group["_is_inflamed"].mean(numeric_only=True).rename("frac_inflamed").reset_index()
    else:
        frac_infl = group.size().rename("frac_inflamed").reset_index()
        frac_infl["frac_inflamed"] = np.nan

    # ---- merge
    out = meta.merge(n_cells, on="sample_id", how="left")
    out = out.merge(means, on="sample_id", how="left")
    out = out.merge(pct, on="sample_id", how="left")
    out = out.merge(mito_pct, on="sample_id", how="left")
    out = out.merge(frac_infl, on="sample_id", how="left")

    # Normalize timepoint if present
    if "timepoint" in out.columns:
        out["timepoint"] = out["timepoint"].map(_norm_timepoint)

    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "state_summary.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path}")

    return out


def _compute_paired_deltas(state_summary: pd.DataFrame, paired_df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    """
    Computes subject-level Post-Pre deltas (one row per subject_id).
    Robustly avoids subtracting categoricals by selecting numeric columns only.
    """
    required = {"sample_id", "subject_id", "timepoint"}
    if not required.issubset(set(paired_df.columns)):
        raise KeyError(
            f"--paired file must include columns: {sorted(required)}\n"
            f"Found columns: {list(paired_df.columns)}\n\n"
            "Tip: for GSE298464, you can use samples.tsv (it has subject_id/timepoint/response)."
        )

    # Join paired metadata onto metrics by sample_id
    merged = paired_df.merge(state_summary, on="sample_id", how="left", suffixes=("", "_state"))
    if merged["subject_id"].isna().any():
        raise KeyError("After merging paired file with state_summary, some rows have missing subject_id.")

    # Choose response grouping column
    if "resp_group" in merged.columns:
        merged["resp_group"] = merged["resp_group"].astype(str)
    elif "response" in merged.columns:
        merged["resp_group"] = merged["response"].astype(str)
    else:
        merged["resp_group"] = "Unknown"

    merged["timepoint"] = merged["timepoint"].map(_norm_timepoint)

    # Numeric metrics we care about (if present)
    want = {
        "pct_high_MITO_ROS": "delta_pct_high_MITO_ROS",
        f"mean_score_{OXPHOS}": "delta_mean_OXPHOS",
        f"mean_score_{ROS}": "delta_mean_ROS",
        "mean_score_HALLMARK_TNFA_SIGNALING_VIA_NFKB": "delta_mean_TNFA",
        "mean_score_HALLMARK_INFLAMMATORY_RESPONSE": "delta_mean_INFLAM",
        "frac_inflamed": "delta_frac_inflamed",
    }
    metric_cols = [c for c in want.keys() if c in merged.columns]

    # Aggregate at subject_id x resp_group x timepoint (handles multiple biopsies/timepoint)
    agg = merged.groupby(["subject_id", "resp_group", "timepoint"], observed=True).agg(
        **{c: (c, "mean") for c in metric_cols},
        n_biopsies=("sample_id", "nunique"),
    ).reset_index()

    # Pivot to get Pre and Post rows
    pre = agg.loc[agg["timepoint"] == "Pre"].set_index("subject_id")
    post = agg.loc[agg["timepoint"] == "Post"].set_index("subject_id")

    subjects = sorted(set(pre.index).intersection(set(post.index)))
    print(f"Paired subjects used: {len(subjects)}")
    if len(subjects) == 0:
        raise ValueError("No subjects have both Pre and Post after processing the paired table.")

    # Keep resp_group from Pre (should match Post)
    resp = pre.loc[subjects, "resp_group"].copy()

    # Only subtract numeric columns
    num_cols = [c for c in metric_cols + ["n_biopsies"] if c in pre.columns and c in post.columns]
    pre_num = pre.loc[subjects, num_cols].apply(pd.to_numeric, errors="coerce")
    post_num = post.loc[subjects, num_cols].apply(pd.to_numeric, errors="coerce")

    delta = post_num - pre_num
    delta = delta.reset_index()

    # Rename delta columns
    rename_map = {}
    for c in metric_cols:
        rename_map[c] = want[c]
    rename_map["n_biopsies"] = "delta_n_biopsies"
    delta = delta.rename(columns=rename_map)

    # Attach resp_group
    delta.insert(1, "resp_group", resp.values)

    # Write
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / "paired_subject_deltas.tsv"
    delta.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path}")

    return delta


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Discovery pipeline: score myeloid states and write summary TSVs."
    )
    ap.add_argument("--h5ad", required=True, help="Path to processed myeloid .h5ad")
    ap.add_argument("--paired", required=False, help="Path to paired/sample table TSV/CSV (must include sample_id, subject_id, timepoint)")
    ap.add_argument("--gmt", required=True, help="Path to hallmark_selected.gmt")
    ap.add_argument("--outdir", default="results/tables", help="Output directory for TSV tables")
    ap.add_argument("--q", type=float, default=0.8, help="Quantile cutoff for state-high (default 0.8 => top 20%)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load
    adata = sc.read_h5ad(args.h5ad)
    if "sample_id" not in adata.obs.columns:
        raise KeyError("Input h5ad must include obs['sample_id'].")

    # Score modules
    modules: Dict[str, List[str]] = load_gmt(args.gmt)

    # Keep only hallmarks we care about and that exist in gmt
    keep = {k: v for k, v in modules.items() if k in HALLMARKS}
    missing_sets = [k for k in HALLMARKS if k not in keep]
    if missing_sets:
        print("Warning: missing gene sets in GMT (will skip):", missing_sets)

    # score_modules is expected to add obs columns: score_<SETNAME>
    score_modules(adata, keep)

    # Ensure scores exist
    _ensure_score_cols(adata, [k for k in HALLMARKS if f"score_{k}" in adata.obs.columns])

    # Ensure state_high columns exist with EXACT expected names
    _ensure_state_high_cols(adata, [k for k in HALLMARKS if f"score_{k}" in adata.obs.columns], groupby="sample_id", q=args.q)

    # Build per-sample table
    state_summary = _build_state_summary(adata, outdir)

    # Optional paired deltas
    if args.paired:
        paired_df = _read_table_anysep(args.paired)

        # Normalize timepoint if present
        if "timepoint" in paired_df.columns:
            paired_df["timepoint"] = paired_df["timepoint"].map(_norm_timepoint)

        _compute_paired_deltas(state_summary, paired_df, outdir)


if __name__ == "__main__":
    main()
