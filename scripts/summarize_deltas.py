from pathlib import Path
import pandas as pd
import numpy as np

def summarize(path: str, label: str) -> pd.DataFrame:
    p = Path(path)
    df = pd.read_csv(p, sep="\t")
    delta_cols = [c for c in df.columns if c.lower().startswith("delta_")]

    rows = []
    for c in delta_cols:
        x = pd.to_numeric(df[c], errors="coerce").dropna().astype(float)
        if len(x) == 0:
            continue
        rows.append({
            "metric": c,
            "n": len(x),
            "mean": x.mean(),
            "median": x.median(),
            "std": x.std(ddof=1) if len(x) > 1 else 0.0,
            "frac_pos": (x > 0).mean(),
            "frac_neg": (x < 0).mean(),
            "min": x.min(),
            "max": x.max(),
        })

    out = pd.DataFrame(rows)
    out = out.sort_values("mean", key=lambda s: s.abs(), ascending=False)

    print("\n" + "=" * 90)
    print(f"{label}: {p.as_posix()}")
    print("=" * 90)
    print(f"Shape: {df.shape}")
    print(f"Delta columns: {len(delta_cols)}")
    if out.empty:
        print("No delta columns found.")
    else:
        print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    return out.set_index("metric")

def main():
    disc = summarize("results/tables/paired_subject_deltas.tsv", "DISCOVERY")
    val  = summarize("results/tables/validation_gse298464/paired_subject_deltas.tsv", "VALIDATION")

    common = disc.index.intersection(val.index)
    if len(common) == 0:
        print("\nNo common delta metrics between discovery and validation.")
        return

    same_sign = (np.sign(disc.loc[common, "mean"]) == np.sign(val.loc[common, "mean"]))
    sign_tbl = pd.DataFrame({
        "disc_mean": disc.loc[common, "mean"],
        "val_mean":  val.loc[common, "mean"],
        "same_sign": same_sign.values
    }).sort_values("disc_mean", key=lambda s: s.abs(), ascending=False)

    print("\n" + "=" * 90)
    print("DISCOVERY vs VALIDATION: sign consistency of MEAN deltas")
    print("=" * 90)
    print(sign_tbl.to_string(float_format=lambda v: f"{v:.4f}"))
    print(f"\nSign-consistent metrics: {int(same_sign.sum())}/{len(common)}")

if __name__ == "__main__":
    main()
