from __future__ import annotations

from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(".")
VALDIR = ROOT / "results" / "tables" / "validation_gse298464"
OUTDIR = ROOT / "results" / "figures"
DOCS_FIGDIR = ROOT / "docs" / "assets" / "figures"

OUTDIR.mkdir(parents=True, exist_ok=True)


def _copy_to_docs(fig_path: Path) -> None:
    if DOCS_FIGDIR.exists():
        DOCS_FIGDIR.mkdir(parents=True, exist_ok=True)
        shutil.copy2(fig_path, DOCS_FIGDIR / fig_path.name)


def _jitter(n: int, width: float = 0.06) -> np.ndarray:
    rng = np.random.default_rng(0)
    return rng.uniform(-width, width, size=n)


def _boxplot(ax, data, labels, **kwargs):
    """Matplotlib 3.9 renamed labels -> tick_labels; keep compatible."""
    try:
        return ax.boxplot(data, tick_labels=labels, **kwargs)
    except TypeError:
        return ax.boxplot(data, labels=labels, **kwargs)


def _pick_col(df: pd.DataFrame, prefer: list[str], contains: list[str] | None = None) -> str | None:
    cols = list(df.columns)
    low = {c.lower(): c for c in cols}

    for p in prefer:
        if p.lower() in low:
            return low[p.lower()]

    if contains:
        for c in cols:
            cl = c.lower()
            if all(k in cl for k in contains):
                return c

    return None


def _mode(series: pd.Series) -> str | None:
    s = series.dropna().astype(str)
    if s.empty:
        return None
    vc = s.value_counts()
    return vc.index[0]


def _norm_resp(x: str) -> str:
    x = str(x).strip()
    xl = x.lower()
    if "non" in xl and "rem" in xl:
        return "Non_Remission"
    if xl.startswith("non"):
        return "Non_Remission"
    if "rem" in xl:
        return "Remission"
    return x


def make_v1() -> Path:
    """
    V1: % MITO_ROS-high by response and timepoint (NR Pre, NR Post, R Pre, R Post)
    Expects: state_summary.tsv with columns:
      - response
      - timepoint
      - pct_high_MITO_ROS
    """
    path = VALDIR / "state_summary.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")

    df = pd.read_csv(path, sep="\t")

    resp_col = _pick_col(df, prefer=["response"], contains=["response"])
    tp_col = _pick_col(df, prefer=["timepoint"], contains=["timepoint"])
    if resp_col is None or tp_col is None:
        raise ValueError("state_summary.tsv must contain 'response' and 'timepoint' columns.")

    if "pct_high_MITO_ROS" not in df.columns:
        raise ValueError("state_summary.tsv must contain 'pct_high_MITO_ROS' column.")

    df = df.copy()
    df[resp_col] = df[resp_col].astype(str).map(_norm_resp)
    df[tp_col] = df[tp_col].astype(str)

    # shorten timepoint display
    def tp_short(x: str) -> str:
        xl = x.lower()
        if "pre" in xl or "base" in xl:
            return "Pre"
        if "post" in xl or "wk" in xl or "week" in xl:
            return "Post"
        return x

    df["tp_short"] = df[tp_col].map(tp_short)

    order = [("Non_Remission", "Pre"), ("Non_Remission", "Post"), ("Remission", "Pre"), ("Remission", "Post")]
    labels = []
    data = []
    grp_for_pts = []

    for r, t in order:
        sub = df[(df[resp_col] == r) & (df["tp_short"] == t)]
        vals = pd.to_numeric(sub["pct_high_MITO_ROS"], errors="coerce").dropna().to_numpy()
        data.append(vals)
        short = "NR" if r == "Non_Remission" else "R"
        labels.append(f"{short} {t} (n={len(vals)})")
        grp_for_pts.append(r)

    fig, ax = plt.subplots(figsize=(10.5, 6.5))

    _boxplot(ax, data, labels, showfliers=False, widths=0.6)

    # overlay points
    for i, (vals, r) in enumerate(zip(data, grp_for_pts), start=1):
        if len(vals) == 0:
            continue
        x = np.full(len(vals), i, dtype=float) + _jitter(len(vals))
        marker = "o" if r == "Non_Remission" else "s"
        ax.scatter(x, vals, marker=marker, s=60, alpha=0.9, label=r if i == 1 or i == 3 else None)

    ax.axhline(0, linewidth=1)
    ax.set_ylabel("% MITO_ROS-high (myeloid)")
    ax.set_title(
        "GSE298464 validation: % MITO_ROS-high by response and timepoint\nPaired delta = post − baseline (see V2)",
        fontsize=16,
    )
    plt.tight_layout()

    outpath = OUTDIR / "figV1_gse298464_mito_ros_by_group.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

    _copy_to_docs(outpath)
    return outpath


def make_v2() -> Path:
    """
    V2: Paired delta (post - baseline) of %MITO_ROS-high by response.
    Expects: paired_subject_deltas.tsv with columns:
      - subject_id
      - delta_pct_high_MITO_ROS
    If 'response' is missing, we infer it from state_summary.tsv and merge by subject_id.
    """
    deltas_path = VALDIR / "paired_subject_deltas.tsv"
    if not deltas_path.exists():
        raise FileNotFoundError(f"Missing {deltas_path}")

    df = pd.read_csv(deltas_path, sep="\t")

    if "subject_id" not in df.columns:
        raise ValueError("paired_subject_deltas.tsv must contain 'subject_id' for V2.")

    if "delta_pct_high_MITO_ROS" not in df.columns:
        raise ValueError("paired_subject_deltas.tsv must contain 'delta_pct_high_MITO_ROS' for V2.")

    df = df.copy()

    # If response is missing, infer from state_summary.tsv
    if "response" not in df.columns:
        state_path = VALDIR / "state_summary.tsv"
        if not state_path.exists():
            raise ValueError("V2 needs 'response'. Add it to paired_subject_deltas.tsv or provide state_summary.tsv.")

        st = pd.read_csv(state_path, sep="\t")

        resp_col = _pick_col(st, prefer=["response"], contains=["response"])
        subj_col = _pick_col(st, prefer=["subject_id"], contains=["subject", "id"])
        if resp_col is None:
            raise ValueError("state_summary.tsv must contain a response column to infer V2 groups.")

        st = st.copy()
        st[resp_col] = st[resp_col].astype(str).map(_norm_resp)

        # Best case: state_summary has subject_id
        if subj_col is not None:
            st[subj_col] = st[subj_col].astype(str)
            resp_map = (
                st.groupby(subj_col)[resp_col]
                .apply(_mode)
                .dropna()
                .rename("response")
                .reset_index()
                .rename(columns={subj_col: "subject_id"})
            )
            df["subject_id"] = df["subject_id"].astype(str)
            df = df.merge(resp_map, on="subject_id", how="left")
        else:
            # Fallback: try to match subject_id as substring of sample_id
            sample_col = _pick_col(st, prefer=["sample_id"], contains=["sample"])
            if sample_col is None:
                raise ValueError("Couldn't infer response: state_summary.tsv has no subject_id or sample_id-like column.")

            st[sample_col] = st[sample_col].astype(str)
            df["subject_id"] = df["subject_id"].astype(str)

            inferred = []
            for sid in df["subject_id"].unique():
                sub = st[st[sample_col].str.contains(sid, regex=False)]
                inferred.append((sid, _mode(sub[resp_col])))

            resp_map = pd.DataFrame(inferred, columns=["subject_id", "response"]).dropna()
            df = df.merge(resp_map, on="subject_id", how="left")

    # normalize response labels
    df["response_norm"] = df["response"].map(_norm_resp) if "response" in df.columns else df["response_norm"].map(_norm_resp)

    # drop subjects we still couldn't label
    df = df.dropna(subset=["response_norm"])

    groups = ["Non_Remission", "Remission"]
    data = []
    labels = []

    for g in groups:
        vals = pd.to_numeric(df.loc[df["response_norm"] == g, "delta_pct_high_MITO_ROS"], errors="coerce").dropna().to_numpy()
        data.append(vals)
        labels.append(f"{g.replace('_',' ')} (n={len(vals)})")

    fig, ax = plt.subplots(figsize=(10.5, 6.5))

    _boxplot(ax, data, labels, showfliers=False, widths=0.6)

    for i, (vals, g) in enumerate(zip(data, groups), start=1):
        if len(vals) == 0:
            continue
        x = np.full(len(vals), i, dtype=float) + _jitter(len(vals))
        marker = "o" if g == "Non_Remission" else "s"
        ax.scatter(x, vals, marker=marker, s=70, alpha=0.9)

    ax.axhline(0, linewidth=1)
    ax.set_ylabel("Δ(Post − Pre) % MITO_ROS-high (myeloid)")
    ax.set_title(
        "GSE298464 validation: MITO_ROS state change by response\nPaired delta = post − baseline",
        fontsize=16,
    )
    plt.tight_layout()

    outpath = OUTDIR / "figV2_gse298464_mito_ros_delta.png"
    fig.savefig(outpath, dpi=150)
    plt.close(fig)

    _copy_to_docs(outpath)
    return outpath


if __name__ == "__main__":
    p1 = make_v1()
    p2 = make_v2()
    print(f"Wrote: {p1}")
    print(f"Wrote: {p2}")
