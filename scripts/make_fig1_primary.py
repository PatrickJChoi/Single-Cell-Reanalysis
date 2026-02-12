from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

TABLE_WITH_META = Path("results/tables/paired_subject_deltas_with_metadata.tsv")
TABLE_BASE = Path("results/tables/paired_subject_deltas.tsv")

FIG_DIR = Path("results/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_FIG1 = FIG_DIR / "fig1_mito_ros_delta.png"
OUT_FIG1B = FIG_DIR / "fig1b_mito_ros_vs_inflammation.png"

def load_deltas():
    if TABLE_WITH_META.exists():
        df = pd.read_csv(TABLE_WITH_META, sep="\t")
        src = TABLE_WITH_META.as_posix()
    else:
        df = pd.read_csv(TABLE_BASE, sep="\t")
        src = TABLE_BASE.as_posix()
    return df, src

def main():
    df, src = load_deltas()

    # Required columns
    ycol = "delta_pct_high_MITO_ROS"
    xcol = "delta_frac_inflamed"

    if ycol not in df.columns:
        raise ValueError(f"Missing required column: {ycol} (source: {src})")

    # ----------------------------
    # Figure 1 (S1): MITO/ROS-high delta by response (if available) else pooled
    # ----------------------------
    plt.figure(figsize=(7.2, 5.0))
    ax = plt.gca()

    y_all = pd.to_numeric(df[ycol], errors="coerce").dropna()

    if "response" in df.columns:
        groups = ["Non_Remission", "Remission"]
        df_plot = df.copy()
        df_plot["response"] = pd.Categorical(df_plot["response"], categories=groups, ordered=True)

        data = []
        labels = []
        for g in groups:
            vals = pd.to_numeric(
                df_plot.loc[df_plot["response"] == g, ycol],
                errors="coerce"
            ).dropna().astype(float).values
            data.append(vals)
            labels.append(f"{g.replace('_',' ')} (n={len(vals)})")

        # Use tick_labels to avoid matplotlib deprecation warning
        ax.boxplot(data, tick_labels=labels, showfliers=False)

        # Add jitter points (two default colors by plotting each group once)
        rng = np.random.default_rng(0)
        for i, vals in enumerate(data, start=1):
            if len(vals) == 0:
                continue
            x = np.full(len(vals), i) + rng.normal(0, 0.05, size=len(vals))
            ax.scatter(x, vals, s=24, alpha=0.9)

        ax.set_title("TAURUS discovery: MITO/ROS state change by response")
    else:
        # Pooled: single box + points
        ax.boxplot([y_all.values], tick_labels=["All subjects"], showfliers=False)
        rng = np.random.default_rng(0)
        x = np.full(len(y_all.values), 1.0) + rng.normal(0, 0.05, size=len(y_all.values))
        ax.scatter(x, y_all.values, s=24, alpha=0.9)
        ax.set_title("TAURUS discovery: MITO/ROS state change (pooled)")

    ax.axhline(0, linewidth=1)
    ax.set_ylabel("Paired delta (post − baseline): % MITO/ROS-high (myeloid)")
    # Move subtitle slightly down so it never feels cramped
    ax.text(0.5, 0.94, "Paired delta = post − baseline", transform=ax.transAxes, ha="center", va="top")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(OUT_FIG1, dpi=250)
    print("Wrote:", OUT_FIG1.as_posix())

    # ----------------------------
    # Figure 1b: MITO/ROS delta vs inflammation delta (color by response if available)
    # ----------------------------
    if xcol not in df.columns:
        print(f"NOTE: {xcol} not found in table; skipping Figure 1b.")
        return

    x = pd.to_numeric(df[xcol], errors="coerce")
    y = pd.to_numeric(df[ycol], errors="coerce")
    mask = x.notna() & y.notna()
    df_sc = df.loc[mask].copy()
    df_sc[xcol] = x.loc[mask].astype(float).values
    df_sc[ycol] = y.loc[mask].astype(float).values

    plt.figure(figsize=(7.2, 5.0))
    ax = plt.gca()

    if "response" in df_sc.columns:
        groups = ["Non_Remission", "Remission"]
        df_sc["response"] = pd.Categorical(df_sc["response"], categories=groups, ordered=True)

        # Plot each group once (Matplotlib assigns default colors automatically)
        for g, marker in [("Non_Remission", "o"), ("Remission", "s")]:
            sub = df_sc[df_sc["response"] == g]
            if len(sub) == 0:
                continue
            ax.scatter(
                sub[xcol],
                sub[ycol],
                s=55,
                marker=marker,
                alpha=0.9,
                label=g.replace("_", " ")
            )
        ax.legend(loc="best")
    else:
        ax.scatter(df_sc[xcol], df_sc[ycol], s=55, alpha=0.9)

    ax.axhline(0, linewidth=1)
    ax.axvline(0, linewidth=1)

    ax.set_title("TAURUS discovery: MITO/ROS change vs inflammation change")
    ax.set_xlabel("Paired delta (post − baseline): fraction inflamed biopsies")
    ax.set_ylabel("Paired delta (post − baseline): % MITO/ROS-high (myeloid)")
    ax.text(0.5, 0.94, "Paired delta = post − baseline", transform=ax.transAxes, ha="center", va="top")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(OUT_FIG1B, dpi=250)
    print("Wrote:", OUT_FIG1B.as_posix())

if __name__ == "__main__":
    main()
