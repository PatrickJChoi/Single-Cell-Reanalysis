import pandas as pd
from pathlib import Path

INPATH = Path("data/validation/gse298464/samples_geo.tsv")
OUTPATH = Path("data/validation/gse298464/samples.tsv")

df = pd.read_csv(INPATH, sep="\t")

# Parse timepoint + response from sample_group (e.g., Pre-R1, Post-NR3)
def parse_timepoint(x: str) -> str:
    x = str(x)
    if x.startswith("Pre-"):
        return "Pre"
    if x.startswith("Post-"):
        return "Post"
    return "Unknown"

def parse_resp(x: str) -> str:
    x = str(x)
    # NR = no-response = Non_Remission
    if "-NR" in x:
        return "Non_Remission"
    # R = response = Remission
    if "-R" in x:
        return "Remission"
    return "Unknown"

def parse_subject(x: str) -> str:
    # Pre-R1 -> R1 ; Post-NR4 -> NR4
    x = str(x)
    x = x.replace("Pre-", "").replace("Post-", "")
    # prefix with UC_ so it won't collide with other diseases later
    return f"UC_{x}"

df["timepoint"] = df["sample_group"].map(parse_timepoint)
df["response"] = df["sample_group"].map(parse_resp)
df["subject_id"] = df["sample_group"].map(parse_subject)

# Standardized fields for your pipeline
out = pd.DataFrame({
    "sample_id": df["sample_code"],     # we use AM5CS### as sample_id
    "subject_id": df["subject_id"],
    "disease": "UC",
    "site": df["tissue"],               # "Colon"
    "timepoint": df["timepoint"],       # Pre/Post
    "response": df["response"],         # Remission/Non_Remission
    "treatment": "IFX",                 # infliximab (anti-TNF)
    "sex": df["sex"],
    "age": df["age"],
    "gsm": df["gsm"],
    "sample_group": df["sample_group"],
})

OUTPATH.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUTPATH, sep="\t", index=False)

print(f"Wrote: {OUTPATH}  rows={out.shape[0]} cols={out.shape[1]}")
print("Unique sample_id:", out["sample_id"].nunique())
print("Unique subjects :", out["subject_id"].nunique())
print("\nSubjects x timepoint table:")
print(out.pivot_table(index="subject_id", columns="timepoint", values="sample_id", aggfunc="count", fill_value=0))
print("\nHead:")
print(out.head(8).to_string(index=False))
