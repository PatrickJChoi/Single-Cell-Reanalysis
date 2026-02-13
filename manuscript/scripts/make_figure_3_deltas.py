import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import os

# -----------------------------
# Setup & Load Data
# -----------------------------
# Create output directory if it doesn't exist
os.makedirs("manuscript/figures", exist_ok=True)

df = pd.read_csv("results/tables/paired_subject_deltas_with_metadata.tsv", sep="\t")

# -----------------------------
# Global Style (Nature-like)
# -----------------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# Force the specific colors: Blue for Remission, Red for Non-Remission
palette = {"Remission": "#2c7bb6", "Non_Remission": "#d7191c"}
order = ["Remission", "Non_Remission"]

fig, axes = plt.subplots(1, 2, figsize=(7, 3.5))

# -----------------------------
# Helper Function for "Money Plots"
# -----------------------------
def plot_with_stat(ax, y_col, title, y_label):
    # Boxplot with specific colors
    sns.boxplot(
        data=df, x="resp_group", y=y_col, hue="resp_group",
        order=order, palette=palette, width=0.5, ax=ax, fliersize=0, 
        linewidth=1, legend=False
    )
    # Strip plot (dots) on top
    sns.stripplot(
        data=df, x="resp_group", y=y_col, hue="resp_group",
        order=order, palette=palette, ax=ax, size=4, jitter=True,
        edgecolor="black", linewidth=0.5, alpha=0.9, legend=False
    )
    
    # Calculate Mann-Whitney U p-value
    g1 = df[df["resp_group"]=="Remission"][y_col]
    g2 = df[df["resp_group"]=="Non_Remission"][y_col]
    stat, p = mannwhitneyu(g1, g2)
    
    # Draw statistical bracket
    y_max = df[y_col].max()
    y_min = df[y_col].min()
    y_range = y_max - y_min
    y_line = y_max + (y_range * 0.1) # Position above max point
    
    # Draw the bracket lines
    ax.plot([0, 0, 1, 1], [y_line, y_line, y_line, y_line], c="black", lw=1)
    ax.text(0.5, y_line, f"p={p:.3f}", ha="center", va="bottom", fontsize=8)
    
    # Aesthetics
    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.set_xlabel("")
    ax.axhline(0, linestyle="--", color="gray", linewidth=1)
    sns.despine(ax=ax)

# -----------------------------
# Generate Subplots
# -----------------------------

# Panel A: The Biological Signal (MITO-High)
plot_with_stat(axes[0], "delta_pct_high_MITO_ROS", "Myeloid State Change", "Δ % MITO-High Cells")

# Panel B: The Clinical Signal (Inflammation or ROS)
# Using 'delta_mean_ROS' matches your specific caption "ROS scores"
plot_with_stat(axes[1], "delta_mean_ROS", "Pathway Change", "Δ ROS Hallmark Score")

plt.tight_layout()

# -----------------------------
# Save
# -----------------------------
output_path = "manuscript/figures/Figure3_DeltaDistributions.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")

print(f"Figure 3 saved to {output_path} with Remission/Non-Remission split!")