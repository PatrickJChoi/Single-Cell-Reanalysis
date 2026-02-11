import os
import yaml
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad

# Needed because pandas 3 uses nullable string arrays; opt-in for writing
ad.settings.allow_write_nullable_strings = True

def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)

def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def qc_basic(adata, qc_cfg: dict, out_dir: str):
    # Annotate mito genes by name (MT- prefix)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    before = adata.n_obs

    adata = adata[adata.obs.n_genes_by_counts >= qc_cfg["min_genes"]].copy()
    adata = adata[adata.obs.total_counts >= qc_cfg["min_counts"]].copy()
    adata = adata[adata.obs.pct_counts_mt <= qc_cfg["max_mito_pct"]].copy()
    adata = adata[adata.obs.n_genes_by_counts <= qc_cfg["max_genes"]].copy()

    after = adata.n_obs

    summary = {
        "cells_before": before,
        "cells_after": after,
        "min_genes": qc_cfg["min_genes"],
        "min_counts": qc_cfg["min_counts"],
        "max_mito_pct": qc_cfg["max_mito_pct"],
        "max_genes": qc_cfg["max_genes"],
    }
    pd.DataFrame([summary]).to_csv(os.path.join(out_dir, "qc_summary.csv"), index=False)

    # QC plots
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    plt.savefig(os.path.join(out_dir, "qc_violin.png"), dpi=200, bbox_inches="tight")
    plt.close()

    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False)
    plt.savefig(os.path.join(out_dir, "qc_counts_vs_mito.png"), dpi=200, bbox_inches="tight")
    plt.close()

    return adata

def preprocess_and_cluster(adata, pp_cfg: dict, out_dir: str):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=pp_cfg["n_hvgs"], subset=True)
    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(
        adata,
        n_comps=pp_cfg["n_pcs"],
        svd_solver="arpack",
        random_state=pp_cfg["random_seed"],
    )
    sc.pp.neighbors(adata, n_pcs=pp_cfg["n_pcs"], random_state=pp_cfg["random_seed"])
    sc.tl.umap(adata, random_state=pp_cfg["random_seed"])
    sc.tl.leiden(
        adata,
        resolution=pp_cfg["leiden_resolution"],
        random_state=pp_cfg["random_seed"],
    )

    sc.pl.umap(adata, color=["leiden"], show=False)
    plt.savefig(os.path.join(out_dir, "umap_leiden.png"), dpi=200, bbox_inches="tight")
    plt.close()

    return adata

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)

    ensure_dir("outputs")
    ensure_dir(cfg["reports"]["out_dir"])

    adata = sc.read(cfg["data"]["input_path"])
    adata = qc_basic(adata, cfg["qc"], cfg["reports"]["out_dir"])
    adata = preprocess_and_cluster(adata, cfg["preprocess"], cfg["reports"]["out_dir"])
    adata.write(cfg["data"]["output_h5ad"])

    print("DONE")
    print("Wrote:", cfg["data"]["output_h5ad"])
    print("Reports in:", cfg["reports"]["out_dir"])

if __name__ == "__main__":
    main()
