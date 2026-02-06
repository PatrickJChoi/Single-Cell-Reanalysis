import sys
from pathlib import Path

# --- make repo root importable so "from src..." works ---
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from anndata import AnnData

from src.state_scoring import load_gmt, score_modules, assign_state_high, add_ribosomal_score


def main() -> None:
    # 1) Dummy AnnData (does NOT need real gene names for this smoketest)
    X = np.random.poisson(1.0, size=(200, 300))
    adata = AnnData(X)
    adata.var_names = [f"gene{i}" for i in range(300)]

    # Minimal metadata
    adata.obs["sample_id"] = ["s1"] * 100 + ["s2"] * 100
    adata.obs["subject_id"] = ["p1"] * 100 + ["p2"] * 100
    adata.obs["timepoint"] = ["pre"] * 50 + ["post"] * 50 + ["pre"] * 50 + ["post"] * 50
    adata.obs["response"] = ["responder"] * 100 + ["nonresponder"] * 100
    adata.obs["celltype"] = ["myeloid"] * 200

    # 2) Load hallmark sets and score
    gmt_path = ROOT / "src" / "gene_sets" / "hallmark_selected.gmt"
    gene_sets = load_gmt(gmt_path)

    score_modules(adata, gene_sets, score_prefix="score_")
    add_ribosomal_score(adata, out_col="score_ribosomal")

    # 3) Assign state-high for each hallmark score (top 20% within sample)
    for set_name in gene_sets.keys():
        assign_state_high(
            adata,
            score_col=f"score_{set_name}",
            groupby="sample_id",
            q=0.80,
            out_col=f"statehigh_{set_name}",
        )

    # 4) Build a minimal summary table (smoketest version)
    rows = []
    for (sample_id, celltype), g in adata.obs.groupby(["sample_id", "celltype"]):
        row = {
            "sample_id": sample_id,
            "celltype": celltype,
            "subject_id": g["subject_id"].iloc[0],
            "timepoint": g["timepoint"].mode().iloc[0],
            "response": g["response"].mode().iloc[0],
            "n_cells": int(len(g)),
            "mean_score_ribosomal": float(pd.to_numeric(g["score_ribosomal"], errors="coerce").mean()),
        }

        for set_name in gene_sets.keys():
            score_col = f"score_{set_name}"
            row[f"mean_{set_name}"] = float(pd.to_numeric(g[score_col], errors="coerce").mean())
            row[f"frac_statehigh_{set_name}"] = float(g[f"statehigh_{set_name}"].mean())

        rows.append(row)

    out = pd.DataFrame(rows)

    # 5) Write output
    outdir = ROOT / "results" / "tables"
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "state_summary.tsv"
    out.to_csv(outpath, sep="\t", index=False)

    print(f"Wrote {outpath}")
    print(out.head())


if __name__ == "__main__":
    main()
