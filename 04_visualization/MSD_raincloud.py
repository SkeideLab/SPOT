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


param = ["sigma"]
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
    fig, axs = plt.subplots(1, 2, figsize=(12, 8), sharex=True)  # Adjust size as needed
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
        for model_idx, model in enumerate(["fetal<29", "fetal>29", "neonates<37", "neonates>37", "12-16y", "18-21y"]):        
            if model =="neonates<37":
                msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_sigma_neonates_young.csv")
            elif model == "neonates>37":
                msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_sigma_neonates_old.csv")
            elif model == "fetal>29":
                msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_sigma_fetal_old.csv")
            elif model == "fetal<29":
                msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_sigma_fetal_young.csv")    
            elif model == "12-16y":
                msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_sigma_HCP_old.csv")
            elif model == "18-21y":
                msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_sigma_HCP_young.csv")   
            Benson = msd_value.loc[:,f"{hemi}_{param}_real_benson"]
            Simulated = msd_value.loc[:,f"{hemi}_{param}_real_simulated"]
            differences = Benson - Simulated
            data_differences = pd.DataFrame({'MSD differences': differences, 'Group': group_name[model_idx]})
            df = pd.concat([df,data_differences], ignore_index=True)


        # Generate raincloud plot
        pt.RainCloud(x="Group", y='MSD differences', data=df, 
                        palette=colors, bw=.2, width_viol=.6, 
                        ax=axs[hemi_index], orient="h", alpha=.65, dodge=True, box_showfliers=False)
        axs[hemi_index].set_title(hemi_lower, fontsize=20)
        axs[hemi_index].set_ylabel(None)
        axs[hemi_index].set_xlabel("$\Delta$ MSD", fontsize=20)  # X-axis label font size

        # Set tick parameters
        if hemi_index == 0:
            axs[hemi_index].tick_params(axis='both', which='major', labelsize=20)  # Major tick labels
            axs[hemi_index].tick_params(axis='both', which='minor', labelsize=16)  # Minor tick labels
        elif hemi_index == 1:
            axs[hemi_index].set_yticks([])
            axs[hemi_index].tick_params(axis='x', which='major', labelsize=20)  # Minor tick labels


    plt.tight_layout()
    #plt.show()
    fig.savefig(f"/data/p_02915/dhcp_derivatives_SPOT/Figures_v3/MSD_{param}.png")