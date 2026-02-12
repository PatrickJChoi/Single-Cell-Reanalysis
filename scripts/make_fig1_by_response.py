from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

IN_TSV = Path("results/tables/paired_subject_deltas_with_metadata.tsv")
OUT_DIR = Path("results/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Overwrite the same filename so docs stays consistent
OUT_PNG = OUT_DIR / "fig1_response_deltas.png"

# Metrics to show (nice, interpretable order)
METRICS = [
    "delta_pct_high_MITO_ROS",
    "delta_mean_ROS",
    "delta_mean_OXPHOS",
    "delta_mean_TNFA",
    "delta_mean_INFLAM",
    "delta_frac_inflamed",
]

# Pretty labels for axis
LABELS = {
    "delta_pct_high_MITO_ROS": "MITO/ROS-high (%)",
    "delta_mean_ROS": "ROS score",
    "delta_mean_OXPHOS": "OXPHOS score",
    "delta_mean_TNFA": "TNFA/NFkB score",
    "delta_mean_INFLAM": "Inflammatory response",
    "delta_frac_inflamed": "Frac. inflamed",
}

GROUPS = ["Non_Remission", "Remission"]

def main():
    df = pd.read_csv(IN_TSV, sep="\t")
    if "response" not in df.columns:
        raise ValueError("Missing 'response' column in paired_subject_deltas_with_metadata.tsv")

    # Keep only available metrics
    metrics = [m for m in METRICS if m in df.columns]
    if not metrics:
        raise ValueError("None of the expected delta metrics were found in the input TSV.")

    # Ensure group order
    df["response"] = pd.Categorical(df["response"], categories=GROUPS, ordered=True)

    # Build per-metric per-group arrays
    data_non = []
    data_rem = []
    for m in metrics:
        x_non = pd.to_numeric(df.loc[df["response"] == "Non_Remission", m], errors="coerce").dropna().astype(float).values
        x_rem = pd.to_numeric(df.loc[df["response"] == "Remission", m], errors="coerce").dropna().astype(float).values
        data_non.append(x_non)
        data_rem.append(x_rem)

    n_non = len(df.loc[df["response"] == "Non_Remission"])
    n_rem = len(df.loc[df["response"] == "Remission"])

    # Plot
    plt.figure(figsize=(10, 4.8))
    ax = plt.gca()

    # Positions: 2 boxes per metric
    gap = 3.0
    pos_non = [i * gap + 1.0 for i in range(len(metrics))]
    pos_rem = [i * gap + 2.0 for i in range(len(metrics))]
    centers = [i * gap + 1.5 for i in range(len(metrics))]

    ax.boxplot(data_non, positions=pos_non, widths=0.65, showfliers=False)
    ax.boxplot(data_rem, positions=pos_rem, widths=0.65, showfliers=False)

        # Jitter points with different markers (no forced colors)
    rng = np.random.default_rng(0)

    xs_non, ys_non = [], []
    xs_rem, ys_rem = [], []

    for i in range(len(metrics)):
        if len(data_non[i]) > 0:
            xs_non.append(np.full(len(data_non[i]), pos_non[i]) + rng.normal(0, 0.06, size=len(data_non[i])))
            ys_non.append(data_non[i])
        if len(data_rem[i]) > 0:
            xs_rem.append(np.full(len(data_rem[i]), pos_rem[i]) + rng.normal(0, 0.06, size=len(data_rem[i])))
            ys_rem.append(data_rem[i])

    if xs_non:
        ax.scatter(np.concatenate(xs_non), np.concatenate(ys_non), s=14, marker="o", alpha=0.9, label="Non_Remission")
    if xs_rem:
        ax.scatter(np.concatenate(xs_rem), np.concatenate(ys_rem), s=14, marker="s", alpha=0.9, label="Remission")

    ax.axhline(0, linewidth=1)
    ax.set_ylabel("Post − baseline (paired delta)")
    ax.set_title(f"Discovery: paired deltas by response group (n={n_non} Non_Remission, n={n_rem} Remission)")

    ax.set_xticks(centers)
    ax.set_xticklabels([LABELS.get(m, m.replace("delta_", "")) for m in metrics], rotation=25, ha="right")

    ax.legend(loc="best")

    ax.axhline(0, linewidth=1)
    ax.set_ylabel("Post − baseline (paired delta)")
    ax.set_title(f"Discovery: paired deltas by response group (n={n_non} Non_Remission, n={n_rem} Remission)")

    ax.set_xticks(centers)
    ax.set_xticklabels([LABELS.get(m, m.replace("delta_", "")) for m in metrics], rotation=25, ha="right")

    ax.legend(loc="best")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=250)
    print("Wrote:", OUT_PNG.as_posix())

if __name__ == "__main__":
    main()
