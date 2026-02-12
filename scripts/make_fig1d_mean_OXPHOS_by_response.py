from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TABLE_META = Path("results/tables/paired_subject_deltas_with_metadata.tsv")
TABLE_BASE = Path("results/tables/paired_subject_deltas.tsv")

OUTFIG = Path("results/figures/fig1d_delta_mean_OXPHOS.png")


def load_df():
    if TABLE_META.exists():
        return pd.read_csv(TABLE_META, sep="\t")
    return pd.read_csv(TABLE_BASE, sep="\t")


def main():
    df = load_df()

    if "response" not in df.columns:
        raise ValueError(
            "No 'response' column found. Run:\n"
            "  python scripts\\add_metadata_to_deltas.py"
        )

    col = "delta_mean_OXPHOS"
    if col not in df.columns:
        raise ValueError(f"Missing column: {col}")

    groups = ["Non_Remission", "Remission"]
    df["response"] = pd.Categorical(df["response"], categories=groups, ordered=True)

    vals = {
        g: pd.to_numeric(df.loc[df["response"] == g, col], errors="coerce")
            .dropna().astype(float).values
        for g in groups
    }

    data = [vals[g] for g in groups]
    labels = [f"{g.replace('_',' ')} (n={len(vals[g])})" for g in groups]

    fig, ax = plt.subplots(figsize=(7.2, 5.0))

    # Remove hollow outlier markers; we plot all points ourselves
    ax.boxplot(data, tick_labels=labels, showfliers=False)

    rng = np.random.default_rng(0)
    markers = {"Non_Remission": "o", "Remission": "s"}

    for i, g in enumerate(groups, start=1):
        y = vals[g]
        if len(y) == 0:
            continue
        x = np.full(len(y), i, dtype=float) + rng.normal(0, 0.06, size=len(y))
        ax.scatter(x, y, s=40, alpha=0.9, marker=markers[g])

    ax.axhline(0, linewidth=1)
    ax.set_title("TAURUS discovery: Mean OXPHOS change by response")
    ax.set_ylabel("Delta(Post - Pre) mean OXPHOS module score")

    ax.text(
        0.5, 0.94, "Paired delta = post − baseline",
        transform=ax.transAxes, ha="center", va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=2),
    )

    OUTFIG.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(OUTFIG, dpi=250)
    print("Wrote:", OUTFIG.as_posix())


if __name__ == "__main__":
    main()
