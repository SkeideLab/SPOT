"""
Statistical comparison of ccf results for six groups
"""
import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
from multiprocessing import Pool, cpu_count
import math
import seaborn as sns
import matplotlib.pyplot as plt
import ptitprince as pt  # For raincloud plots
from statannotations.Annotator import Annotator
import os

# Get the full path of the script
file_path = os.path.dirname(os.path.abspath(__file__))



param = ["eccentricity", "polarangle"]
hemispheres = ["L", "R"]

colors = {
    'tri-2 \n prenatal': 'limegreen',
    'tri-3 \n prenatal': 'mediumaquamarine',
    'preterm \n neonatal': 'gold',
    'full-term \n neonatal': 'yellowgreen',
    'adolescent': 'cadetblue',
    'adult': 'steelblue'
}
group_name = ['tri-2 \n prenatal',
    'tri-3 \n prenatal',
    'preterm \n neonatal',
    'full-term \n neonatal',
    'adolescent',
    'adult',
    ]
for param_index, param in enumerate(param):
    fig, axs = plt.subplots(2, 2, figsize=(12, 12), sharex=True)  # Adjust size as needed
    for hemi_index, hemi in enumerate(hemispheres):
        if hemi == "L":
            hemi_lower="Left"
        elif hemi == "R":
            hemi_lower="Right"
        # Initialize lists to hold combined data for both "raw" and "combat"
        combined_group_names = []
        combined_correlation_values = []
        combined_cov_type = []  # To distinguish between raw and combat data
        df = pd.DataFrame()
        for area_index, area in enumerate(["V2", "V3"]):
            for model_idx, model in enumerate(["neonates_young","neonates_old", "fetal_young", "fetal_old", "adolescent", "adult"]):        
                msd_value=pd.read_csv(f'{file_path}/MSD_{model}_{area}.csv')
                Benson = msd_value.loc[:,f"{hemi}_{param}_real_benson"]
                Simulated = msd_value.loc[:,f"{hemi}_{param}_real_simulated"]
                differences = Benson - Simulated
                print(differences)
                data_differences = pd.DataFrame({'MSD differences': differences, 'Group': group_name[model_idx]})
                df = pd.concat([df,data_differences], ignore_index=True)

            # Add vertical line at x=0
            axs[area_index, hemi_index].axvline(x=0, color='k', linestyle='--', linewidth=1)

            # Generate raincloud plot
            pt.RainCloud(x="Group", y='MSD differences', data=df, 
                            palette=colors, bw=.2, width_viol=.6, 
                            ax=axs[area_index, hemi_index], orient="h", alpha=.65, dodge=True, box_showfliers=False)
            if area_index ==0:
                axs[area_index, hemi_index].set_title(hemi_lower, fontsize=18)

            axs[area_index, hemi_index].set_xlabel("$\Delta$ MSD", fontsize=18)  # X-axis label font size
            axs[area_index, hemi_index].set_ylabel(None)

            # Set tick parameters
            if hemi_index == 0:
                axs[area_index, hemi_index].tick_params(axis='both', which='major', labelsize=18)  # Major tick labels
                axs[area_index, hemi_index].tick_params(axis='both', which='minor', labelsize=16)  # Minor tick labels
            elif hemi_index == 1:
                axs[area_index, hemi_index].set_yticks([])
                axs[area_index, hemi_index].tick_params(axis='x', which='major', labelsize=18)  # Minor tick labels

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.text(0.56, 0.95, "V2", ha='center', va='center', fontsize=24, weight='bold')
    fig.text(0.56, 0.486, "V3", ha='center', va='center', fontsize=24, weight='bold')
    plt.subplots_adjust(hspace=0.16)
    #plt.show()
    if param =="eccentricity":
        fig_name = "S8"
    elif param == "polarangle":
        fig_name = "S11"
    fig.savefig(f"{file_path}/Figure{fig_name}.png")