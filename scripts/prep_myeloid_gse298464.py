#!/usr/bin/env python
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


def _force_plain_strings(adata: sc.AnnData) -> None:
    """
    Prevent pandas nullable StringArray from breaking .write_h5ad()
    by converting obs/var index + string columns to plain python strings (object dtype).
    """
    # Opt-in (safe in *our* environment; also used as a fallback if any StringArray remains)
    try:
        ad.settings.allow_write_nullable_strings = True
    except Exception:
        pass

    # obs index
    adata.obs.index = pd.Index(
        adata.obs.index.astype(str).to_numpy(dtype=object),
        dtype=object,
        name=adata.obs.index.name,
    )

    # var index
    adata.var.index = pd.Index(
        adata.var.index.astype(str).to_numpy(dtype=object),
        dtype=object,
        name=adata.var.index.name,
    )

    # obs columns
    for col in adata.obs.columns:
        if pd.api.types.is_string_dtype(adata.obs[col].dtype):
            adata.obs[col] = adata.obs[col].astype(str)

    # var columns (rare, but safe)
    for col in adata.var.columns:
        if pd.api.types.is_string_dtype(adata.var[col].dtype):
            adata.var[col] = adata.var[col].astype(str)


def main():
    ap = argparse.ArgumentParser(
        description="QC + marker-score + extract myeloid cells for GSE298464."
    )
    ap.add_argument("--in_h5ad", required=True, help="Input combined raw h5ad")
    ap.add_argument("--out_h5ad", required=True, help="Output myeloid-only h5ad")
    ap.add_argument("--min_genes", type=int, default=200)
    ap.add_argument("--max_genes", type=int, default=6000)
    ap.add_argument("--max_mt", type=float, default=20.0, help="Max percent mitochondrial counts")
    ap.add_argument("--min_cells", type=int, default=3, help="Min cells a gene must appear in")
    args = ap.parse_args()

    print("Reading:", args.in_h5ad)
    adata = sc.read_h5ad(args.in_h5ad)

    # Basic gene filtering
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    # Mito genes (human convention)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Cell QC
    qc_mask = (
        (adata.obs["n_genes_by_counts"] >= args.min_genes)
        & (adata.obs["n_genes_by_counts"] <= args.max_genes)
        & (adata.obs["pct_counts_mt"] <= args.max_mt)
    )
    adata = adata[qc_mask].copy()
    print("After QC shape:", adata.shape)

    # Normalize + log for marker scoring
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # A simple myeloid marker score
    # (keeps it robust across datasets without needing full clustering)
    myeloid_markers = [
        "LYZ", "LST1", "TYROBP", "FCER1G", "AIF1", "S100A8", "S100A9",
        "CTSS", "LGALS3", "MS4A7"
    ]
    present = [g for g in myeloid_markers if g in adata.var_names]
    if len(present) < 3:
        raise RuntimeError(f"Too few myeloid markers found in var_names. Found: {present}")

    sc.tl.score_genes(adata, gene_list=present, score_name="myeloid_score")

    # Pick myeloid cells by score cutoff (quantile-based = portable)
    cutoff = float(np.quantile(adata.obs["myeloid_score"], 0.80))
    adata_my = adata[adata.obs["myeloid_score"] >= cutoff].copy()

    # Basic prints
    print("Myeloid cells kept:", adata_my.n_obs)
    print("Myeloid sample_id nunique:", adata_my.obs["sample_id"].nunique())
    print("Myeloid subject_id nunique:", adata_my.obs["subject_id"].nunique())
    print("Myeloid timepoint counts:\n", adata_my.obs["timepoint"].value_counts())

    per = adata_my.obs["sample_id"].value_counts()
    print("Per-sample myeloid counts (min/median/max):", int(per.min()), int(per.median()), int(per.max()))

    # IMPORTANT: make writing stable across pandas/anndata versions
    _force_plain_strings(adata_my)

    out = args.out_h5ad
    os.makedirs(os.path.dirname(out), exist_ok=True)
    adata_my.write_h5ad(out)
    print("Wrote:", out)
    print("Final shape (cells, genes):", adata_my.shape)


if __name__ == "__main__":
    main()
