import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

# -----------------------------
# Setup
# -----------------------------
os.makedirs("manuscript/figures", exist_ok=True)

# -----------------------------
# Load Data
# -----------------------------
df = pd.read_csv("results/tables/paired_subject_deltas_with_metadata.tsv", sep="\t")

# CRITICAL STEP: Sort by Response Group so Remission patients are grouped together
df = df.sort_values("resp_group", ascending=False) # Remission first, then Non-Remission

# -----------------------------
# Select Features
# -----------------------------
features = [
    "delta_pct_high_MITO_ROS",
    "delta_mean_ROS",
    "delta_mean_OXPHOS",
    "delta_mean_TNFA",
    "delta_mean_INFLAM",
]

heatmap_data = df[features]
response_group = df["resp_group"]

# -----------------------------
# Styling (Nature-Style)
# -----------------------------
plt.close("all")
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9
})

# Create Figure with extra space on the left for the color bar
fig, ax = plt.subplots(figsize=(5, 6))

# Draw Heatmap
# RdBu_r is better than coolwarm for "Delta" (Blue=Negative/Drop, Red=Positive/Increase)
sns.heatmap(
    heatmap_data,
    cmap="RdBu_r",
    center=0,
    linewidths=0.0, 
    cbar_kws={"label": "Δ Score", "shrink": 0.8},
    yticklabels=False, # Hide messy subject IDs
    ax=ax
)

# -----------------------------
# Add The "Response" Color Bar
# -----------------------------
# This draws a colored strip on the left to show which rows are Remission vs Non-Remission
for i, (idx, row_resp) in enumerate(response_group.items()):
    color = "#2c7bb6" if row_resp == "Remission" else "#d7191c" # Blue vs Red
    # Add a rectangle patch to the left of the heatmap
    rect = patches.Rectangle((-0.05, i), 0.05, 1, transform=ax.get_yaxis_transform(), 
                             clip_on=False, facecolor=color)
    ax.add_patch(rect)

# Add label for the color bar
ax.text(-0.06, 0.5, "Response Group", rotation=90, va="center", ha="right", 
        transform=ax.transAxes, fontsize=9, fontweight="bold")

# -----------------------------
# Labels & Save
# -----------------------------
ax.set_title("Pathway-Level Paired Deltas\n(Sorted by Response)", pad=20)
ax.set_xlabel("")
ax.set_ylabel("Subjects")

plt.tight_layout()

plt.savefig(
    "manuscript/figures/Figure2_Heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

print("Figure 2 saved! (Sorted by Remission/Non-Remission)")