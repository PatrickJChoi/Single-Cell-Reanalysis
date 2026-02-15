"""Shared fixtures for unit tests."""

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import pytest


@pytest.fixture
def toy_adata():
    """AnnData with 50 cells x 600 genes, including MT- and RPL/RPS genes.

    Uses 600 genes so scanpy's calculate_qc_metrics (percent_top=[50,100,200,500])
    does not raise IndexError.
    """
    np.random.seed(42)
    n_cells, n_genes = 50, 600
    n_mt = 5
    n_ribo = 10

    mt_names = [f"MT-{i}" for i in range(n_mt)]
    ribo_names = [f"RPL{i}" for i in range(n_ribo // 2)] + [
        f"RPS{i}" for i in range(n_ribo // 2)
    ]
    other_names = [f"GENE{i:03d}" for i in range(n_genes - n_mt - n_ribo)]
    var_names = mt_names + ribo_names + other_names

    X = np.random.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32)
    # Elevate MT counts for a few "stressed" cells
    X[:10, :n_mt] += 20

    adata = AnnData(X)
    adata.var_names = pd.Index(var_names)
    adata.obs_names = pd.Index([f"cell_{i:03d}" for i in range(n_cells)])
    adata.obs["sample_id"] = (["sampleA"] * 25) + (["sampleB"] * 25)

    return adata


@pytest.fixture
def normalized_adata(toy_adata):
    """toy_adata after normalization + log1p (required for score_genes)."""
    adata = toy_adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


@pytest.fixture
def tmp_gmt(tmp_path):
    """Write a small GMT file and return its path."""
    gmt = tmp_path / "test.gmt"
    # Gene set with enough genes to pass the >=10 threshold
    genes_big = "\t".join([f"GENE{i:03d}" for i in range(20)])
    # Gene set that is too small (<10 genes found in adata)
    genes_small = "\t".join([f"RARE{i}" for i in range(5)])
    gmt.write_text(
        f"SET_BIG\tdescription\t{genes_big}\n"
        f"SET_SMALL\tdescription\t{genes_small}\n"
    )
    return gmt
