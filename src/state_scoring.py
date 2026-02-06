
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd

import scanpy as sc
from anndata import AnnData


def load_gmt(gmt_path: str | Path) -> Dict[str, List[str]]:
    """
    Load a GMT file into a dict: {set_name: [gene1, gene2, ...]}.
    GMT format: set_name<TAB>description<TAB>gene1<TAB>gene2...
    """
    gmt_path = Path(gmt_path)
    gene_sets: Dict[str, List[str]] = {}
    with gmt_path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            set_name = parts[0]
            genes = [g for g in parts[2:] if g]
            gene_sets[set_name] = genes
    if not gene_sets:
        raise ValueError(f"No gene sets parsed from {gmt_path}")
    return gene_sets


def score_modules(
    adata: AnnData,
    gene_sets: Dict[str, Sequence[str]],
    score_prefix: str = "score_",
    use_raw: bool = False,
) -> None:
    """
    Adds module score columns to adata.obs named f"{score_prefix}{set_name}".
    Uses scanpy.tl.score_genes.
    """
    # Ensure gene names exist
    var_names = set(adata.raw.var_names if (use_raw and adata.raw is not None) else adata.var_names)

    for set_name, genes in gene_sets.items():
        genes_in = [g for g in genes if g in var_names]
        if len(genes_in) < 10:
            # still score, but warn by writing NaNs (keeps schema stable)
            adata.obs[f"{score_prefix}{set_name}"] = np.nan
            continue

        sc.tl.score_genes(
            adata,
            gene_list=genes_in,
            score_name=f"{score_prefix}{set_name}",
            use_raw=use_raw,
            copy=False,
        )


def add_ribosomal_score(adata: AnnData, out_col: str = "score_ribosomal") -> None:
    """RPL*/RPS* gene prefix score (simple control)."""
    genes = [g for g in adata.var_names if g.startswith("RPL") or g.startswith("RPS")]
    if len(genes) < 10:
        adata.obs[out_col] = np.nan
        return
    sc.tl.score_genes(adata, gene_list=genes, score_name=out_col, use_raw=False, copy=False)


def add_cell_cycle_scores(
    adata: AnnData,
    s_genes: Sequence[str],
    g2m_genes: Sequence[str],
    s_col: str = "score_S",
    g2m_col: str = "score_G2M",
) -> None:
    """
    Adds cell cycle scores using scanpy's canonical method.
    You will provide S and G2M gene lists later (can use scanpy defaults or a file).
    """
    # scanpy uses: sc.tl.score_genes_cell_cycle, but requires present genes
    s_in = [g for g in s_genes if g in adata.var_names]
    g2m_in = [g for g in g2m_genes if g in adata.var_names]
    if len(s_in) < 10 or len(g2m_in) < 10:
        adata.obs[s_col] = np.nan
        adata.obs[g2m_col] = np.nan
        return
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_in, g2m_genes=g2m_in)
    # scanpy writes 'S_score' and 'G2M_score' by default; copy into our names
    adata.obs[s_col] = adata.obs["S_score"]
    adata.obs[g2m_col] = adata.obs["G2M_score"]


def assign_state_high(
    adata: AnnData,
    score_col: str,
    groupby: str = "sample_id",
    q: float = 0.80,
    out_col: str | None = None,
) -> None:
    """
    Adds boolean column to adata.obs: True if score is in top q quantile within each groupby sample.
    """
    if out_col is None:
        out_col = f"statehigh_{score_col}"
    if groupby not in adata.obs.columns:
        raise KeyError(f"{groupby} not found in adata.obs")

    def _mark(group: pd.DataFrame) -> pd.Series:
        x = group[score_col]
        if x.isna().all():
            return pd.Series([False] * len(group), index=group.index)
        thr = np.nanquantile(x.values, q)
        return x >= thr

    adata.obs[out_col] = adata.obs.groupby(groupby, group_keys=False).apply(_mark)
