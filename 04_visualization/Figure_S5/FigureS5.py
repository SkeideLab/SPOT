#!/usr/bin/env python3
"""
Load roi_diff_results.csv and generate RainCloud plot.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ptitprince as pt
import os

# Get the full path of the script
file_path = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------
# Load data
# ----------------------------------------------------
df = pd.read_csv(f"{file_path}/roi_diff_results.csv")

# ----------------------------------------------------
# Color and label mapping
# ----------------------------------------------------
colors = {
    'tri-2 \n prenatal': 'limegreen',
    'tri-3 \n prenatal': 'mediumaquamarine',
    'preterm \n neonatal': 'gold',
    'fullterm \n neonatal': 'yellowgreen',
    'adolescent': 'cadetblue',
    'adult': 'steelblue'
}

group_color_map = {
    '2nd': 'tri-2 \n prenatal',
    '3rd': 'tri-3 \n prenatal',
    'preterm': 'preterm \n neonatal',
    'fullterm': 'fullterm \n neonatal',
    'adolescent': 'adolescent',
    'adult': 'adult'
}

df["Group"] = df["Group"].map(group_color_map)

# ----------------------------------------------------
# Compute % within threshold
# ----------------------------------------------------
percent_pass = {}
for g in df["Group"].unique():
    subset = df[df["Group"] == g]
    percent_pass[g] = np.mean(subset["MeanROI_diff"] <= 0.1) * 100

# ----------------------------------------------------
# Plot RainCloud
# ----------------------------------------------------
plt.figure(figsize=(12, 6))
ax = plt.gca()

pt.RainCloud(
    x="Group",
    y="MeanROI_diff",
    data=df,
    palette=colors,
    bw=.2,
    width_viol=.6,
    ax=ax,
    orient="h",
    alpha=.65,
    dodge=True,
    box_showfliers=False
)

# Add % text
for i, g in enumerate(df["Group"].unique()):
    ax.text(1.02, i, f"{percent_pass[g]:.1f}%", va="center", fontsize=12)

ax.set_xlabel("Averaged $R^2$ difference (first - second)", fontsize=15)
ax.set_ylabel(None)
ax.axvline(x=0.1, color='k', linestyle='--', linewidth=1)
ax.set_xlim(-0.5, 2)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.savefig(f"{file_path}/Fig_S5.png", dpi=300, bbox_inches='tight')
plt.show()
