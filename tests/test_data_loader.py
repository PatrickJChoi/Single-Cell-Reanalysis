"""Tests for src/data_loader.py."""

import numpy as np
import pytest
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.data_loader import load_ann_data, save_ann_data


class TestLoadAnnData:
    def test_load_h5ad(self, toy_adata, tmp_path):
        path = tmp_path / "test.h5ad"
        toy_adata.write_h5ad(path)
        loaded = load_ann_data(path)
        assert loaded.shape == toy_adata.shape
        assert list(loaded.var_names) == list(toy_adata.var_names)

    def test_load_csv(self, toy_adata, tmp_path):
        path = tmp_path / "test.csv"
        import pandas as pd

        df = pd.DataFrame(
            toy_adata.X, index=toy_adata.obs_names, columns=toy_adata.var_names
        )
        df.to_csv(path)
        loaded = load_ann_data(path)
        assert loaded.shape == toy_adata.shape

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_ann_data(tmp_path / "nonexistent.h5ad")


class TestSaveAnnData:
    def test_round_trip(self, toy_adata, tmp_path):
        path = tmp_path / "out.h5ad"
        save_ann_data(toy_adata, path)
        loaded = load_ann_data(path)
        assert loaded.shape == toy_adata.shape
        np.testing.assert_array_almost_equal(loaded.X, toy_adata.X)
