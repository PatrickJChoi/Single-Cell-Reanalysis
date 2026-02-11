import os
import pandas as pd
import matplotlib.pyplot as plt

OUTDIR = r"results\figures"
TABLE  = r"results\tables\validation_gse298464\paired_subject_deltas.tsv"

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    d = pd.read_csv(TABLE, sep="\t")
    if "resp_group" not in d.columns or "delta_pct_high_MITO_ROS" not in d.columns:
        raise KeyError(f"Missing expected columns in {TABLE}. Found: {list(d.columns)}")

    nr = d.loc[d["resp_group"] == "Non_Remission", "delta_pct_high_MITO_ROS"].dropna()
    r  = d.loc[d["resp_group"] == "Remission",     "delta_pct_high_MITO_ROS"].dropna()

    print("Validation deltas:",
          "Non_Remission n =", len(nr), "median =", float(nr.median()),
          "| Remission n =", len(r), "median =", float(r.median()))

    plt.figure()
    plt.boxplot([nr, r], tick_labels=["Non_Remission", "Remission"], showfliers=True)
    plt.axhline(0, linewidth=1)
    plt.ylabel("Δ(Post − Pre) % MITO_ROS-high (myeloid)")
    plt.title("GSE298464 validation: MITO_ROS state change by response")

    outpath = os.path.join(OUTDIR, "figV2_gse298464_mito_ros_delta.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()

    print("Saved:", outpath)

if __name__ == "__main__":
    main()
