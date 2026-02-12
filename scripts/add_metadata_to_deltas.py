from pathlib import Path
import pandas as pd

DISC_PAIRED = Path("results/tables/paired_subject_deltas.tsv")
DISC_SAMPLES = Path("data/discovery/taurus_v3/samples.tsv")
DISC_OUT = Path("results/tables/paired_subject_deltas_with_metadata.tsv")

VAL_PAIRED = Path("results/tables/validation_gse298464/paired_subject_deltas.tsv")
VAL_SAMPLES = Path("results/tables/validation_gse298464/samples.tsv")  # if exists
VAL_OUT = Path("results/tables/validation_gse298464/paired_subject_deltas_with_metadata.tsv")

def read_table(path: Path) -> pd.DataFrame:
    # Try tab first, then whitespace
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path, sep=r"\s+")

def normalize_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]
    return df

def build_subject_metadata(samples: pd.DataFrame) -> pd.DataFrame:
    """
    Reduce sample-level metadata to subject-level.
    For columns with a single unique value per subject, keep it.
    For response-like columns, keep the most frequent non-null value.
    """
    samples = samples.copy()
    if "subject_id" not in samples.columns:
        raise ValueError("samples table must contain a 'subject_id' column")

    # Prefer these columns if present
    preferred = ["response", "Disease", "Inflammation", "Site", "Gender", "Age", "Ethnicity"]
    keep_cols = ["subject_id"] + [c for c in preferred if c in samples.columns]

    # If response isn't present, still try to keep any column that looks like it
    if "response" not in keep_cols:
        for c in samples.columns:
            if c.lower() in ("response", "responder", "group", "outcome"):
                keep_cols.append(c)
                break

    meta = samples[keep_cols].copy()

    def pick_value(s: pd.Series):
        s = s.dropna()
        if len(s) == 0:
            return None
        # most frequent value
        return s.value_counts().idxmax()

    out = meta.groupby("subject_id").agg({c: pick_value for c in meta.columns if c != "subject_id"}).reset_index()
    return out

def merge_and_write(paired_path: Path, samples_path: Path, out_path: Path):
    if not paired_path.exists():
        print("Missing:", paired_path)
        return

    paired = normalize_cols(read_table(paired_path))
    if "subject_id" not in paired.columns:
        raise ValueError(f"{paired_path} must have a subject_id column")

    if not samples_path.exists():
        print(f"NOTE: samples file not found ({samples_path}). Writing unchanged copy with no metadata.")
        paired.to_csv(out_path, sep="\t", index=False)
        print("Wrote:", out_path.as_posix())
        return

    samples = normalize_cols(read_table(samples_path))

    # If samples.tsv is missing subject_id but has sample_id, we can’t map reliably here.
    if "subject_id" not in samples.columns:
        raise ValueError(f"{samples_path} does not contain subject_id, so we cannot merge metadata.")

    subj_meta = build_subject_metadata(samples)

    merged = paired.merge(subj_meta, on="subject_id", how="left")

    # Put metadata columns right after subject_id
    cols = list(merged.columns)
    meta_cols = [c for c in subj_meta.columns if c != "subject_id"]
    ordered = ["subject_id"] + meta_cols + [c for c in cols if c not in (["subject_id"] + meta_cols)]
    merged = merged[ordered]

    merged.to_csv(out_path, sep="\t", index=False)
    print("Wrote:", out_path.as_posix())
    # quick coverage report
    for c in meta_cols:
        frac = merged[c].notna().mean()
        print(f"  {c}: {frac:.2%} non-null")

def main():
    print("== Discovery ==")
    merge_and_write(DISC_PAIRED, DISC_SAMPLES, DISC_OUT)

    print("\n== Validation ==")
    # validation samples file may not exist; script will handle
    merge_and_write(VAL_PAIRED, VAL_SAMPLES, VAL_OUT)

if __name__ == "__main__":
    main()
