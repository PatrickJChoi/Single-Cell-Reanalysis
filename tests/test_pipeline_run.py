"""Tests for src/pipeline/run.py (QC + preprocessing functions)."""

import numpy as np
import pandas as pd
import scanpy as sc
import pytest
import sys, os
from unittest.mock import patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib
matplotlib.use("Agg")

from src.pipeline.run import load_config, qc_basic, preprocess_and_cluster


def _noop(*args, **kwargs):
    """No-op stub to replace plotting functions."""
    pass


# Mock out scanpy plotting and plt.savefig/plt.close to avoid
# seaborn/pandas version incompatibility in the test environment.
_plot_patches = [
    patch("src.pipeline.run.sc.pl.violin", _noop),
    patch("src.pipeline.run.sc.pl.scatter", _noop),
    patch("src.pipeline.run.sc.pl.umap", _noop),
    patch("src.pipeline.run.plt.savefig", _noop),
    patch("src.pipeline.run.plt.close", _noop),
]


@pytest.fixture(autouse=True)
def _mock_plots():
    """Disable all plotting for pipeline tests."""
    for p in _plot_patches:
        p.start()
    yield
    for p in _plot_patches:
        p.stop()


# ── load_config ──────────────────────────────────────────────────────

class TestLoadConfig:
    def test_loads_smoke_config(self):
        cfg = load_config("configs/smoke.yaml")
        assert "data" in cfg
        assert "qc" in cfg
        assert "preprocess" in cfg
        assert cfg["qc"]["min_genes"] == 50

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_config("configs/does_not_exist.yaml")


# ── qc_basic ─────────────────────────────────────────────────────────

class TestQcBasic:
    @pytest.fixture
    def qc_cfg(self):
        return {
            "min_genes": 5,
            "min_counts": 10,
            "max_mito_pct": 50,
            "max_genes": 600,
        }

    def test_filters_cells(self, toy_adata, qc_cfg, tmp_path):
        """QC should keep cells that pass all thresholds."""
        n_before = toy_adata.n_obs
        result = qc_basic(toy_adata, qc_cfg, str(tmp_path))
        assert 0 < result.n_obs <= n_before

    def test_qc_summary_csv_written(self, toy_adata, qc_cfg, tmp_path):
        qc_basic(toy_adata, qc_cfg, str(tmp_path))
        csv_path = tmp_path / "qc_summary.csv"
        assert csv_path.exists()
        df = pd.read_csv(csv_path)
        assert "cells_before" in df.columns
        assert "cells_after" in df.columns
        assert df["cells_before"].iloc[0] == toy_adata.n_obs

    def test_strict_filter_removes_all(self, toy_adata, tmp_path):
        """Extremely strict filter should remove most/all cells."""
        strict_cfg = {
            "min_genes": 9999,
            "min_counts": 10,
            "max_mito_pct": 50,
            "max_genes": 99999,
        }
        result = qc_basic(toy_adata, strict_cfg, str(tmp_path))
        assert result.n_obs == 0

    def test_mito_annotation(self, toy_adata, qc_cfg, tmp_path):
        """QC should annotate MT genes in var."""
        assert "mt" not in toy_adata.var.columns
        qc_basic(toy_adata, qc_cfg, str(tmp_path))
        assert "mt" in toy_adata.var.columns
        assert toy_adata.var["mt"].sum() == 5  # 5 MT- genes in toy_adata


# ── preprocess_and_cluster ───────────────────────────────────────────

class TestPreprocessAndCluster:
    @pytest.fixture
    def pp_cfg(self):
        return {
            "n_hvgs": 30,
            "n_pcs": 10,
            "leiden_resolution": 0.5,
            "random_seed": 0,
        }

    @pytest.fixture
    def qc_passed_adata(self, toy_adata, tmp_path):
        """Run QC first so adata has required annotations."""
        qc_cfg = {
            "min_genes": 5,
            "min_counts": 10,
            "max_mito_pct": 90,
            "max_genes": 600,
        }
        return qc_basic(toy_adata, qc_cfg, str(tmp_path))

    def test_adds_leiden_clusters(self, qc_passed_adata, pp_cfg, tmp_path):
        result = preprocess_and_cluster(qc_passed_adata, pp_cfg, str(tmp_path))
        assert "leiden" in result.obs.columns
        assert result.obs["leiden"].nunique() >= 1

    def test_adds_umap(self, qc_passed_adata, pp_cfg, tmp_path):
        result = preprocess_and_cluster(qc_passed_adata, pp_cfg, str(tmp_path))
        assert "X_umap" in result.obsm
        assert result.obsm["X_umap"].shape == (result.n_obs, 2)

    def test_adds_pca(self, qc_passed_adata, pp_cfg, tmp_path):
        result = preprocess_and_cluster(qc_passed_adata, pp_cfg, str(tmp_path))
        assert "X_pca" in result.obsm
        assert result.obsm["X_pca"].shape[1] == pp_cfg["n_pcs"]

    def test_deterministic_clustering(self, qc_passed_adata, pp_cfg, tmp_path):
        """Same seed should produce same clusters."""
        dir1 = tmp_path / "run1"
        dir1.mkdir()
        dir2 = tmp_path / "run2"
        dir2.mkdir()

        r1 = preprocess_and_cluster(qc_passed_adata.copy(), pp_cfg, str(dir1))
        r2 = preprocess_and_cluster(qc_passed_adata.copy(), pp_cfg, str(dir2))
        pd.testing.assert_series_equal(
            r1.obs["leiden"].reset_index(drop=True),
            r2.obs["leiden"].reset_index(drop=True),
        )
