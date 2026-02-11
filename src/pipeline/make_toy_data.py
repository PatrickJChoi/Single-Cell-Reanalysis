import numpy as np
import scanpy as sc
import anndata as ad
ad.settings.allow_write_nullable_strings = True
from pathlib import Path

def main():
    np.random.seed(0)

    n_cells = 120
    n_genes = 600
    n_mt = 30  # number of mitochondrial genes

    # Create gene names (some MT-*)
    gene_names = [f"GENE{i:04d}" for i in range(n_genes - n_mt)]
    mt_names = [f"MT-{i:03d}" for i in range(n_mt)]
    var_names = np.array(mt_names + gene_names, dtype=str)

    # Poisson counts with slightly higher MT for some cells
    X = np.random.poisson(lam=1.2, size=(n_cells, n_genes)).astype(np.float32)
    stressed = np.random.choice(n_cells, size=20, replace=False)
    X[stressed, :n_mt] += np.random.poisson(lam=3.0, size=(20, n_mt)).astype(np.float32)

    adata = sc.AnnData(X)
    adata.var_names = var_names
    adata.obs_names = [f"cell_{i:03d}" for i in range(n_cells)]

    out_path = Path("data/raw/example.h5ad")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_path.as_posix())
    print(f"Wrote toy dataset to {out_path} with shape {adata.shape}")

if __name__ == "__main__":
    main()
