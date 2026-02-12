from pathlib import Path
import pandas as pd
import numpy as np

def read_tsv(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t")

def summarize_paired(df: pd.DataFrame, label: str):
    print("\n" + "="*80)
    print(label)
    print("="*80)
    print(f"Shape: {df.shape}")

    # Identify ID-ish columns
    id_cols = [c for c in ["subject_id", "sample_id", "timepoint", "response"] if c in df.columns]
    if id_cols:
        print("ID columns present:", id_cols)

    # Response counts
    if "response" in df.columns:
        print("\nResponse counts:")
        print(df["response"].value_counts(dropna=False).to_string())
    else:
        print("\nNo 'response' column found in this table.")

    # Delta columns (numeric + contains 'delta' somewhere)
    delta_cols = []
    for c in df.columns:
        if "delta" in c.lower():
            try:
                if pd.api.types.is_numeric_dtype(df[c]):
                    delta_cols.append(c)
            except Exception:
                pass

    print("\nDelta columns found:", delta_cols if delta_cols else "(none found)")

    # Group means/medians by response
    if delta_cols and "response" in df.columns:
        means = df.groupby("response")[delta_cols].mean(numeric_only=True).round(4)
        medians = df.groupby("response")[delta_cols].median(numeric_only=True).round(4)
        print("\nGroup MEANS (by response):")
        print(means.to_string())
        print("\nGroup MEDIANS (by response):")
        print(medians.to_string())

        # If exactly 2 groups, compute simple effect sizes (Cohen's d) per metric
        groups = [g for g in df["response"].dropna().unique()]
        if len(groups) == 2:
            g1, g2 = groups[0], groups[1]
            rows = []
            for c in delta_cols:
                x = df.loc[df["response"] == g1, c].dropna().astype(float)
                y = df.loc[df["response"] == g2, c].dropna().astype(float)
                if len(x) >= 2 and len(y) >= 2:
                    pooled = np.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2)
                    d = np.nan if pooled == 0 else (x.mean() - y.mean()) / pooled
                    rows.append([c, float(x.mean()), float(y.mean()), float(x.mean() - y.mean()), float(d), len(x), len(y)])
            if rows:
                out = pd.DataFrame(rows, columns=[f"metric", f"mean_{g1}", f"mean_{g2}", "diff", "cohens_d", f"n_{g1}", f"n_{g2}"])
                out = out.reindex(out["diff"].abs().sort_values(ascending=False).index)
                print(f"\nTop absolute differences (two-group, {g1} vs {g2}):")
                print(out.head(10).round(4).to_string(index=False))
    else:
        print("\nSkipping group summary (need delta columns + 'response').")

def main():
    # Discovery paired deltas
    disc_paired = Path("results/tables/paired_subject_deltas.tsv")
    if disc_paired.exists():
        df = read_tsv(disc_paired)
        summarize_paired(df, "DISCOVERY: results/tables/paired_subject_deltas.tsv")
    else:
        print("Missing:", disc_paired)

    # Discovery state summary (prints top rows + columns)
    disc_state = Path("results/tables/state_summary.tsv")
    if disc_state.exists():
        df2 = read_tsv(disc_state)
        print("\n" + "="*80)
        print("DISCOVERY: results/tables/state_summary.tsv")
        print("="*80)
        print(f"Shape: {df2.shape}")
        print("Columns:", list(df2.columns))
        print("\nHead (first 10 rows):")
        print(df2.head(10).to_string(index=False))
    else:
        print("Missing:", disc_state)

    # Validation (try to find paired/delta-like tables)
    val_dir = Path("results/tables/validation_gse298464")
    if val_dir.exists():
        print("\n" + "="*80)
        print("VALIDATION: results/tables/validation_gse298464/")
        print("="*80)
        tsvs = sorted(val_dir.glob("*.tsv"))
        if not tsvs:
            print("No TSV files found in validation folder.")
            return
        print("TSV files found:")
        for p in tsvs:
            print(" -", p.name)

        # Prefer a paired/deltas table if present
        candidates = [p for p in tsvs if "paired" in p.name.lower() or "delta" in p.name.lower()]
        if candidates:
            p = candidates[0]
            try:
                dfv = read_tsv(p)
                summarize_paired(dfv, f"VALIDATION: {p.as_posix()}")
            except Exception as e:
                print("Could not read candidate validation table:", p.name, "error:", e)
        else:
            print("No obvious paired/delta TSV found for validation summary.")
    else:
        print("\nValidation folder not found (this is OK):", val_dir)

if __name__ == "__main__":
    main()
