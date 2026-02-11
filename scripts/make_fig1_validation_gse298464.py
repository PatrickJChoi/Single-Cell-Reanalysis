import os
import pandas as pd
import matplotlib.pyplot as plt

OUTDIR = r"results\figures"
TABLE  = r"results\tables\validation_gse298464\state_summary.tsv"

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    s = pd.read_csv(TABLE, sep="\t")

    # Key column
    col = "pct_high_MITO_ROS"
    if col not in s.columns:
        raise KeyError(f"Missing {col} in {TABLE}. Found: {list(s.columns)}")

    # Use response + timepoint
    if "response" not in s.columns or "timepoint" not in s.columns:
        raise KeyError("state_summary.tsv must have 'response' and 'timepoint' columns.")

    # Split into groups/timepoints
    pre_nr  = s.loc[(s["response"] == "Non_Remission") & (s["timepoint"] == "Pre"),  col].dropna()
    post_nr = s.loc[(s["response"] == "Non_Remission") & (s["timepoint"] == "Post"), col].dropna()
    pre_r   = s.loc[(s["response"] == "Remission") & (s["timepoint"] == "Pre"),  col].dropna()
    post_r  = s.loc[(s["response"] == "Remission") & (s["timepoint"] == "Post"), col].dropna()

    print("Counts:",
          "NR Pre", len(pre_nr), "NR Post", len(post_nr),
          "| R Pre", len(pre_r), "R Post", len(post_r))

    # Figure: 4-group boxplot
    plt.figure()
    plt.boxplot(
        [pre_nr, post_nr, pre_r, post_r],
        tick_labels=["NR Pre", "NR Post", "R Pre", "R Post"],
        showfliers=True
    )
    plt.ylabel("% MITO_ROS-high (myeloid)")
    plt.title("GSE298464 validation: MITO_ROS-high by response and timepoint")
    outpath = os.path.join(OUTDIR, "figV1_gse298464_mito_ros_by_group.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()

    print("Saved:", outpath)

if __name__ == "__main__":
    main()
