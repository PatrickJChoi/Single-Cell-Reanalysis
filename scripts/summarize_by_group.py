import pandas as pd
import numpy as np
from pathlib import Path

PATH = Path("results/tables/paired_subject_deltas_with_metadata.tsv")

def main():
    df = pd.read_csv(PATH, sep="\t")

    if "response" not in df.columns:
        raise ValueError("No 'response' column found. Did you create paired_subject_deltas_with_metadata.tsv?")

    delta_cols = [c for c in df.columns if c.lower().startswith("delta_")]
    groups = [g for g in df["response"].dropna().unique()]

    print("File:", PATH.as_posix())
    print("Rows:", len(df))
    print("Response groups:", groups)
    print("Delta columns:", delta_cols)

    # Group means/medians
    means = df.groupby("response")[delta_cols].mean(numeric_only=True).round(4)
    medians = df.groupby("response")[delta_cols].median(numeric_only=True).round(4)

    print("\n" + "="*90)
    print("GROUP MEANS (by response)")
    print("="*90)
    print(means.to_string())

    print("\n" + "="*90)
    print("GROUP MEDIANS (by response)")
    print("="*90)
    print(medians.to_string())

    # If exactly 2 groups, compute effect sizes (Cohen's d)
    if len(groups) == 2:
        g1, g2 = groups[0], groups[1]
        rows = []
        for c in delta_cols:
            x = pd.to_numeric(df.loc[df["response"] == g1, c], errors="coerce").dropna().astype(float)
            y = pd.to_numeric(df.loc[df["response"] == g2, c], errors="coerce").dropna().astype(float)
            if len(x) < 2 or len(y) < 2:
                continue
            pooled = np.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2)
            d = np.nan if pooled == 0 else (x.mean() - y.mean()) / pooled
            rows.append([c, x.mean(), y.mean(), x.mean() - y.mean(), d, len(x), len(y)])

        if rows:
            out = pd.DataFrame(rows, columns=["metric", f"mean_{g1}", f"mean_{g2}", "diff", "cohens_d", f"n_{g1}", f"n_{g2}"])
            out = out.reindex(out["diff"].abs().sort_values(ascending=False).index)
            print("\n" + "="*90)
            print(f"TOP DIFFERENCES (two-group: {g1} vs {g2})")
            print("="*90)
            print(out.head(10).round(4).to_string(index=False))
        else:
            print("\nNot enough data for effect sizes (need >=2 samples per group).")
    else:
        print("\nNote: effect sizes computed only when there are exactly 2 response groups.")

if __name__ == "__main__":
    main()
