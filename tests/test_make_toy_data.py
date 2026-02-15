"""Tests for src/pipeline/make_toy_data.py."""

import numpy as np
import scanpy as sc
import pytest
import subprocess
import sys
from pathlib import Path


class TestMakeToyData:
    def test_generates_file(self, tmp_path, monkeypatch):
        """Running make_toy_data.py creates an h5ad file with expected shape."""
        monkeypatch.chdir(tmp_path)
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve().parents[1] / "src" / "pipeline" / "make_toy_data.py")],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr
        out_path = tmp_path / "data" / "raw" / "example.h5ad"
        assert out_path.exists()

        adata = sc.read_h5ad(out_path)
        assert adata.shape == (120, 600)

    def test_expected_gene_names(self, tmp_path, monkeypatch):
        """Output should include MT- prefixed genes and GENE* genes."""
        monkeypatch.chdir(tmp_path)
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve().parents[1] / "src" / "pipeline" / "make_toy_data.py")],
            capture_output=True,
        )
        adata = sc.read_h5ad(tmp_path / "data" / "raw" / "example.h5ad")
        mt_genes = [g for g in adata.var_names if g.startswith("MT-")]
        assert len(mt_genes) == 30

    def test_deterministic_output(self, tmp_path, monkeypatch):
        """Two runs should produce identical data (seeded RNG)."""
        script = str(Path(__file__).resolve().parents[1] / "src" / "pipeline" / "make_toy_data.py")

        dir1 = tmp_path / "run1"
        dir1.mkdir()
        monkeypatch.chdir(dir1)
        subprocess.run([sys.executable, script], capture_output=True)

        dir2 = tmp_path / "run2"
        dir2.mkdir()
        monkeypatch.chdir(dir2)
        subprocess.run([sys.executable, script], capture_output=True)

        a1 = sc.read_h5ad(dir1 / "data" / "raw" / "example.h5ad")
        a2 = sc.read_h5ad(dir2 / "data" / "raw" / "example.h5ad")
        np.testing.assert_array_equal(a1.X, a2.X)

    def test_stressed_cells_have_higher_mt(self, tmp_path, monkeypatch):
        """Some cells should have elevated MT counts (stressed cells)."""
        monkeypatch.chdir(tmp_path)
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve().parents[1] / "src" / "pipeline" / "make_toy_data.py")],
            capture_output=True,
        )
        adata = sc.read_h5ad(tmp_path / "data" / "raw" / "example.h5ad")
        mt_mask = adata.var_names.str.startswith("MT-")
        mt_sums = np.array(adata.X[:, mt_mask].sum(axis=1)).ravel()
        # There should be a spread — some cells with noticeably higher MT
        assert mt_sums.max() > np.median(mt_sums) * 1.5
