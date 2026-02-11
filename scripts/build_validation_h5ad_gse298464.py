#!/usr/bin/env python3
"""
Build a combined AnnData (.h5ad) from GEO raw 10x matrices for GSE298464.

Expected raw_dir contents (per sample):
  GSM9014962_AM5CS126_barcodes.tsv.gz
  GSM9014962_AM5CS126_features.tsv.gz
  GSM9014962_AM5CS126_matrix.mtx.gz
... etc

Expected samples_tsv columns (from your generated samples.tsv):
  sample_id, subject_id, disease, site, timepoint, response, treatment, sex, age, gsm, sample_group
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse

import anndata as ad
from anndata import AnnData


def read_lines_gz(path: Path) -> list[str]:
    with gzip.open(path, "rt", errors="ignore") as f:
        return [line.rstrip("\n\r") for line in f]


def read_features_gz(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    features.tsv.gz can be:
      gene_id \t gene_name
    or:
      gene_id \t gene_name \t feature_type
    Use gene_name as var_names and keep gene_id in var['gene_id'].
    """
    df = pd.read_csv(path, sep="\t", header=None, compression="gzip")
    if df.shape[1] < 2:
        raise ValueError(f"features file has <2 columns: {path}")
    gene_id = df.iloc[:, 0].astype(str).to_numpy()
    gene_name = df.iloc[:, 1].astype(str).to_numpy()
    return gene_id, gene_name


def read_mtx_gz(path: Path) -> sparse.spmatrix:
    with gzip.open(path, "rb") as f:
        m = mmread(f)
    if not sparse.issparse(m):
        m = sparse.csr_matrix(m)
    else:
        m = m.tocsr()
    return m


def load_one_sample(raw_dir: Path, gsm: str, sample_id: str, meta: dict) -> AnnData:
    prefix = f"{gsm}_{sample_id}"
    barcodes_p = raw_dir / f"{prefix}_barcodes.tsv.gz"
    features_p = raw_dir / f"{prefix}_features.tsv.gz"
    matrix_p = raw_dir / f"{prefix}_matrix.mtx.gz"

    for p in [barcodes_p, features_p, matrix_p]:
        if not p.exists():
            raise FileNotFoundError(f"Missing expected file: {p}")

    barcodes = read_lines_gz(barcodes_p)
    gene_id, gene_name = read_features_gz(features_p)
    m = read_mtx_gz(matrix_p)

    # 10x mtx is typically (genes x cells). We want (cells x genes).
    if m.shape[0] == len(gene_name) and m.shape[1] == len(barcodes):
        X = m.T
    elif m.shape[1] == len(gene_name) and m.shape[0] == len(barcodes):
        X = m
    else:
        raise ValueError(
            f"Matrix shape {m.shape} doesn't match features ({len(gene_name)}) and barcodes ({len(barcodes)}). "
            f"Files: {matrix_p}"
        )

    # Make cell IDs unique across samples by suffixing with sample_id
    obs_names = [f"{bc}-{sample_id}" for bc in barcodes]

    # IMPORTANT: keep index name None to avoid /obs writing issues on some versions
    obs_df = pd.DataFrame(index=pd.Index(obs_names, name=None))
    var_df = pd.DataFrame(index=pd.Index(gene_name, name="gene_symbol"))

    adata = AnnData(X=X, obs=obs_df, var=var_df)
    adata.var["gene_id"] = gene_id

    # Add metadata to every cell (cast to plain python types)
    for k, v in meta.items():
        if v is None or (isinstance(v, float) and np.isnan(v)):
            adata.obs[k] = ""
        else:
            adata.obs[k] = v

    adata.var_names_make_unique()
    return adata


def make_write_safe(adata: AnnData) -> None:
    """
    Prevent anndata write_h5ad errors with Pandas 3 string arrays / nullable strings.
    - Force obs index to object strings and remove its name
    - Force any pandas 'string' dtype columns -> object
    - Same for var index
    """
    adata.obs = adata.obs.copy()
    adata.obs.index = adata.obs.index.astype("object")
    adata.obs.index.name = None

    for c in adata.obs.columns:
        dt = adata.obs[c].dtype
        # Pandas 3 may use "string[python]" / "string[pyarrow]"
        if pd.api.types.is_string_dtype(dt) or str(dt).startswith("string"):
            adata.obs[c] = adata.obs[c].astype("object").fillna("")
        # also keep object columns clean
        elif pd.api.types.is_object_dtype(dt):
            adata.obs[c] = adata.obs[c].astype("object")

    adata.var = adata.var.copy()
    adata.var.index = adata.var.index.astype("object")


def main():
    parser = argparse.ArgumentParser(description="Build combined .h5ad for GSE298464 from raw 10x mtx files.")
    parser.add_argument("--raw_dir", required=True, help="Folder containing GSM*_barcodes/features/matrix.mtx.gz files")
    parser.add_argument("--samples_tsv", required=True, help="samples.tsv with one row per sample")
    parser.add_argument("--out", required=True, help="Output .h5ad path")
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir).resolve()
    samples_tsv = Path(args.samples_tsv).resolve()
    out = Path(args.out).resolve()

    if not raw_dir.exists():
        raise FileNotFoundError(f"--raw_dir not found: {raw_dir}")
    if not samples_tsv.exists():
        raise FileNotFoundError(f"--samples_tsv not found: {samples_tsv}")

    # Read all as strings to avoid accidental pandas nullable-string behavior
    samples = pd.read_csv(samples_tsv, sep="\t", dtype=str)

    required_cols = ["sample_id", "gsm"]
    for c in required_cols:
        if c not in samples.columns:
            raise ValueError(f"samples.tsv missing required column '{c}'. Found cols: {list(samples.columns)}")

    adatas: list[AnnData] = []

    for _, row in samples.iterrows():
        sample_id = str(row["sample_id"])
        gsm = str(row["gsm"])

        print(f"Loading {sample_id} ...")

        meta = {k: (row[k] if k in row.index else "") for k in samples.columns}
        meta["sample_id"] = sample_id
        meta["gsm"] = gsm

        adata = load_one_sample(raw_dir=raw_dir, gsm=gsm, sample_id=sample_id, meta=meta)
        make_write_safe(adata)
        adatas.append(adata)

    print("Concatenating...")

    combo = ad.concat(adatas, axis=0, join="outer", merge="first")

    combo.var_names_make_unique()
    make_write_safe(combo)

    # KEY FIX: allow writing pandas StringArray (pandas 3) into h5ad
    ad.settings.allow_write_nullable_strings = True

    out.parent.mkdir(parents=True, exist_ok=True)
    combo.write_h5ad(out)

    print(f"Wrote: {out}")
    print("Final shape (cells, genes):", combo.shape)


if __name__ == "__main__":
    main()
