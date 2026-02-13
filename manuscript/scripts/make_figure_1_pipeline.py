import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

# Create directory if missing
os.makedirs("manuscript/figures", exist_ok=True)

# ---------------------------
# Global Nature-style settings
# ---------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})

# INCREASED WIDTH to 16 to fit the last box
fig, ax = plt.subplots(figsize=(16, 3))
ax.set_xlim(0, 16) 
ax.set_ylim(0, 3)
ax.axis("off")

def draw_box(x, y, text, facecolor):
    box = FancyBboxPatch(
        (x, y), 2.2, 1.0,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1,
        edgecolor="black",
        facecolor=facecolor
    )
    ax.add_patch(box)
    ax.text(x + 1.1, y + 0.5, text,
            ha="center", va="center", fontsize=9)

def draw_arrow(x1, x2):
    # Arrow connects right side of box 1 to left side of box 2
    arrow = FancyArrowPatch(
        (x1, 1.5), (x2, 1.5),
        arrowstyle="->",
        linewidth=1,
        color="black",
        mutation_scale=15  # Makes the arrow head slightly bigger/clearer
    )
    ax.add_patch(arrow)

# Node positions (Spread out slightly more if needed, but these work with width 16)
xs = [0.5, 3, 5.5, 8, 10.5, 13]

# Draw nodes
draw_box(xs[0], 1, "Raw scRNA-seq Data\n(Taurus v3 Cohort)", "white")
draw_box(xs[1], 1, "QC & Filtering\n(>200 genes, <20% Mito)", "#f0f0f0")
draw_box(xs[2], 1, "Hallmark Scoring\n(ROS & OXPHOS Sets)", "#f0f0f0")
draw_box(xs[3], 1, "State Definition\n(Top 10% MITO-High)", "#f0f0f0")
draw_box(xs[4], 1, "Paired Delta Calc\n(Post - Baseline)", "#f0f0f0")
draw_box(xs[5], 1, "Statistical Stratification\n(Remission vs. Non-Remission)", "#cfe8f3")

# Draw arrows
for i in range(len(xs) - 1):
    draw_arrow(xs[i] + 2.2, xs[i+1])

plt.tight_layout()
plt.savefig("manuscript/figures/figure_1_pipeline.png", dpi=300, bbox_inches="tight")
plt.close()

print("Figure 1 saved with corrected width!")