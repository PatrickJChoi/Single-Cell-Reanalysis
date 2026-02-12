import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


TABLE = "results/tables/validation_gse298464/state_summary.tsv"
OUTDIR = "results/figures"


def _find_col(df: pd.DataFrame, target: str) -> str:
    """Case-insensitive column finder."""
    for c in df.columns:
        if c.strip().lower() == target.strip().lower():
            return c
    raise KeyError(f"Required column '{target}' not found. Available: {list(df.columns)}")


def _find_pct_mito_ros_col(df: pd.DataFrame) -> str:
    # Prefer exact match used elsewhere in the repo
    for cand in ["pct_high_MITO_ROS", "pct_high_mito_ros"]:
        for c in df.columns:
            if c.strip().lower() == cand.lower():
                return c

    # Fallback: any pct_high*mito*ros column
    for c in df.columns:
        lc = c.lower()
        if "pct" in lc and "high" in lc and "mito" in lc and "ros" in lc:
            return c

    raise KeyError(
        "Could not find a MITO/ROS percent-high column. "
        "Expected something like 'pct_high_MITO_ROS'. "
        f"Available: {list(df.columns)}"
    )


def _norm_response(x: str) -> str:
    s = str(x).strip().lower().replace(" ", "_")
    if "non" in s and "remission" in s:
        return "Non_Remission"
    if s in {"nr", "nonremission", "non_remission"}:
        return "Non_Remission"
    if s in {"r", "remission"}:
        return "Remission"
    # keep original-ish if unexpected
    return str(x).strip()


def _norm_timepoint(x: str) -> str:
    s = str(x).strip().lower().replace(" ", "_")
    if s in {"pre", "baseline", "t0"}:
        return "Pre"
    if s in {"post", "after", "t1"}:
        return "Post"
    return str(x).strip()


def main() -> None:
    if not os.path.exists(TABLE):
        raise FileNotFoundError(f"Missing validation table: {TABLE}")

    df = pd.read_csv(TABLE, sep="\t")

    # required columns (case-insensitive)
    response_col = _find_col(df, "response")
    timepoint_col = _find_col(df, "timepoint")
    y_col = _find_pct_mito_ros_col(df)

    # normalize values so ordering is stable even if case differs
    df = df.copy()
    df["_response"] = df[response_col].map(_norm_response)
    df["_timepoint"] = df[timepoint_col].map(_norm_timepoint)

    order = [
        ("Non_Remission", "Pre", "NR Pre"),
        ("Non_Remission", "Post", "NR Post"),
        ("Remission", "Pre", "R Pre"),
        ("Remission", "Post", "R Post"),
    ]

    data = []
    counts = []
    for r, t, _ in order:
        vals = pd.to_numeric(
            df.loc[(df["_response"] == r) & (df["_timepoint"] == t), y_col],
            errors="coerce",
        ).dropna()
        data.append(vals.values)
        counts.append(int(vals.shape[0]))

    # Plot
    os.makedirs(OUTDIR, exist_ok=True)
    outpath = os.path.join(OUTDIR, "figV1_gse298464_mito_ros_by_group.png")

    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)

    positions = np.arange(1, len(order) + 1)

    # Boxplot WITHOUT hollow outlier markers
    ax.boxplot(
        data,
        tick_labels=[f"{lab} (n={n})" for (_, _, lab), n in zip(order, counts)],
        positions=positions,
        showfliers=False,     # <-- this removes hollow circles
        widths=0.55,
    )

    # Overlay filled scatter points (so you still see all samples)
    rng = np.random.default_rng(0)
    for i, ((resp, tp, _), vals) in enumerate(zip(order, data), start=1):
        if len(vals) == 0:
            continue
        jitter = rng.normal(0, 0.05, size=len(vals))

        # style by response
        if resp == "Non_Remission":
            marker = "o"
            face = "C0"
        else:
            marker = "s"
            face = "C1"

        ax.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            s=45,
            marker=marker,
            facecolors=face,
            edgecolors=face,
            alpha=0.9,
            zorder=3,
        )

    ax.axhline(0, linewidth=1)

    ax.set_ylabel("% MITO_ROS-high (myeloid)")
    ax.set_title(
        "GSE298464 validation: % MITO_ROS-high by response and timepoint\n"
        "Paired delta = post − baseline",
        pad=14,
    )

    # Legend (response only)
    h1 = ax.scatter([], [], marker="o", s=60, facecolors="C0", edgecolors="C0", label="Non_Remission")
    h2 = ax.scatter([], [], marker="s", s=60, facecolors="C1", edgecolors="C1", label="Remission")
    ax.legend(handles=[h1, h2], loc="upper right", frameon=True)

    fig.savefig(outpath, dpi=200)
    plt.close(fig)

    print(f"Wrote: {outpath}")


if __name__ == "__main__":
    main()
