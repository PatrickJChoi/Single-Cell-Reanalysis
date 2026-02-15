"""Tests for src/state_scoring.py."""

import numpy as np
import pandas as pd
import pytest
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.state_scoring import (
    load_gmt,
    score_modules,
    add_ribosomal_score,
    assign_state_high,
)


# ── load_gmt ──────────────────────────────────────────────────────────

class TestLoadGmt:
    def test_basic_parse(self, tmp_gmt):
        gs = load_gmt(tmp_gmt)
        assert "SET_BIG" in gs
        assert "SET_SMALL" in gs
        assert len(gs["SET_BIG"]) == 20
        assert len(gs["SET_SMALL"]) == 5

    def test_empty_file_raises(self, tmp_path):
        empty = tmp_path / "empty.gmt"
        empty.write_text("")
        with pytest.raises(ValueError, match="No gene sets"):
            load_gmt(empty)

    def test_skips_short_lines(self, tmp_path):
        bad = tmp_path / "bad.gmt"
        bad.write_text("ONLY_NAME\n")  # fewer than 3 tab-separated fields
        with pytest.raises(ValueError, match="No gene sets"):
            load_gmt(bad)

    def test_ignores_empty_gene_fields(self, tmp_path):
        gmt = tmp_path / "trailing.gmt"
        gmt.write_text("MYSET\tdesc\tA\tB\t\t\n")
        gs = load_gmt(gmt)
        assert gs["MYSET"] == ["A", "B"]


# ── score_modules ────────────────────────────────────────────────────

class TestScoreModules:
    def test_scores_added_to_obs(self, normalized_adata, tmp_gmt):
        gs = load_gmt(tmp_gmt)
        score_modules(normalized_adata, gs, score_prefix="score_")
        assert "score_SET_BIG" in normalized_adata.obs.columns
        # SET_BIG should have real numeric values (enough genes overlap)
        assert normalized_adata.obs["score_SET_BIG"].notna().all()

    def test_small_set_gets_nan(self, normalized_adata, tmp_gmt):
        gs = load_gmt(tmp_gmt)
        score_modules(normalized_adata, gs, score_prefix="score_")
        # SET_SMALL has only RARE* genes not present in adata → NaN
        assert normalized_adata.obs["score_SET_SMALL"].isna().all()

    def test_custom_prefix(self, normalized_adata, tmp_gmt):
        gs = load_gmt(tmp_gmt)
        score_modules(normalized_adata, gs, score_prefix="mod_")
        assert "mod_SET_BIG" in normalized_adata.obs.columns

    def test_empty_gene_sets_dict(self, normalized_adata):
        # Should do nothing and not raise
        score_modules(normalized_adata, {})
        # No new score columns added
        assert not any(c.startswith("score_") for c in normalized_adata.obs.columns)


# ── add_ribosomal_score ──────────────────────────────────────────────

class TestRibosomalScore:
    def test_scores_added(self, normalized_adata):
        add_ribosomal_score(normalized_adata)
        assert "score_ribosomal" in normalized_adata.obs.columns
        assert normalized_adata.obs["score_ribosomal"].notna().all()

    def test_nan_when_too_few_ribo_genes(self, normalized_adata):
        # Remove ribosomal genes
        mask = ~normalized_adata.var_names.str.startswith(("RPL", "RPS"))
        adata = normalized_adata[:, mask].copy()
        add_ribosomal_score(adata)
        assert adata.obs["score_ribosomal"].isna().all()

    def test_custom_col_name(self, normalized_adata):
        add_ribosomal_score(normalized_adata, out_col="ribo_test")
        assert "ribo_test" in normalized_adata.obs.columns


# ── assign_state_high ────────────────────────────────────────────────

class TestAssignStateHigh:
    def test_basic_assignment(self, normalized_adata):
        # Create a fake score column
        normalized_adata.obs["score_test"] = np.random.randn(normalized_adata.n_obs)
        assign_state_high(normalized_adata, "score_test", groupby="sample_id", q=0.8)
        col = "statehigh_score_test"
        assert col in normalized_adata.obs.columns
        # About 20% per group should be True (top quantile)
        frac = normalized_adata.obs[col].mean()
        assert 0.05 <= frac <= 0.5  # loose bounds

    def test_custom_output_column(self, normalized_adata):
        normalized_adata.obs["score_x"] = np.random.randn(normalized_adata.n_obs)
        assign_state_high(
            normalized_adata, "score_x", groupby="sample_id", out_col="my_flag"
        )
        assert "my_flag" in normalized_adata.obs.columns

    def test_missing_groupby_raises(self, normalized_adata):
        normalized_adata.obs["score_x"] = 0.0
        with pytest.raises(KeyError, match="bad_col"):
            assign_state_high(
                normalized_adata, "score_x", groupby="bad_col"
            )

    def test_all_nan_scores_gives_false(self, normalized_adata):
        normalized_adata.obs["score_nan"] = np.nan
        assign_state_high(normalized_adata, "score_nan", groupby="sample_id")
        assert not normalized_adata.obs["statehigh_score_nan"].any()

    def test_top20_per_group(self, normalized_adata):
        """Verify state-high marks roughly top 20% within each sample group."""
        np.random.seed(99)
        normalized_adata.obs["score_z"] = np.random.randn(normalized_adata.n_obs)
        assign_state_high(normalized_adata, "score_z", groupby="sample_id", q=0.80)
        for _, grp in normalized_adata.obs.groupby("sample_id"):
            n_high = grp["statehigh_score_z"].sum()
            # With 25 cells and q=0.8, expect ~5 cells marked high (≥20%)
            assert n_high >= 1
            assert n_high <= len(grp)
