from pathlib import Path
from gseapy import Msigdb

TARGET = [
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]

DBVER = "2023.1.Hs"   # human MSigDB release
CATEGORY = "h.all"    # hallmark collection

outdir = Path("src/gene_sets")
outdir.mkdir(parents=True, exist_ok=True)
outfile = outdir / "hallmark_selected.gmt"

msig = Msigdb()
gmt = msig.get_gmt(category=CATEGORY, dbver=DBVER)

missing = [k for k in TARGET if k not in gmt]
if missing:
    # Show a few keys to help debug if versions/categories differ
    preview = list(gmt.keys())[:20]
    raise ValueError(
        f"Missing sets in MSigDB download: {missing}\n"
        f"First 20 available keys: {preview}"
    )

with outfile.open("w", encoding="utf-8") as f:
    for set_name in TARGET:
        genes = gmt[set_name]
        f.write("\t".join([set_name, f"MSigDB {DBVER} {CATEGORY}"] + genes) + "\n")

print(f"Wrote {len(TARGET)} sets to {outfile}")
