import argparse
import gzip
from pathlib import Path

import pandas as pd


def parse_family_soft(soft_gz: str) -> pd.DataFrame:
    samples = []
    cur = None

    def flush():
        nonlocal cur
        if cur is None:
            return
        row = {"gsm": cur["gsm"]}
        for k, v in cur.get("kv", {}).items():
            row[k] = v
        if "title" in cur:
            row["title"] = cur["title"]
        if "source_name" in cur:
            row["source_name"] = cur["source_name"]
        if cur.get("characteristics_free"):
            row["characteristics_free"] = " | ".join(cur["characteristics_free"])
        samples.append(row)
        cur = None

    with gzip.open(soft_gz, "rt", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")

            if line.startswith("^SAMPLE = "):
                flush()
                gsm = line.split("=", 1)[1].strip()
                cur = {"gsm": gsm, "kv": {}, "characteristics_free": []}
                continue

            if cur is None:
                continue

            if line.startswith("!Sample_title = "):
                cur["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                cur["source_name"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                val = line.split("=", 1)[1].strip()
                if ":" in val:
                    k, v = val.split(":", 1)
                    k = k.strip().lower().replace(" ", "_")
                    v = v.strip()
                    cur["kv"][k] = v
                else:
                    cur["characteristics_free"].append(val)

    flush()
    return pd.DataFrame(samples)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--soft", required=True, help="Path to GSE*_family.soft.gz")
    ap.add_argument("--manifest", required=True, help="raw_manifest.tsv (gsm, sample_code)")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    df = parse_family_soft(args.soft)
    man = pd.read_csv(args.manifest, sep="\t")
    df = df.merge(man, on="gsm", how="left")

    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(outpath, sep="\t", index=False)

    print(f"Wrote: {outpath}  rows={df.shape[0]} cols={df.shape[1]}")
    print("First columns:", list(df.columns)[:30])
    print(df.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
