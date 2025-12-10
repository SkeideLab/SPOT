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

def boferroni_correction(p):
    corrected_p_values = p * 15  # Bonferroni adjustment
    corrected_p_values = np.minimum(p, 1)  # Cap at 1.0
    return corrected_p_values

def flatten(arr):
    return arr.flatten()


# Create a figure with subplots
areas =  ["V2d","V2v",]
hemispheres = ["L", "R"]
axes = ["S"]
colors = ["gold", "yellowgreen", "limegreen", "mediumaquamarine", "cadetblue", "steelblue"]
# Initialize lists to hold combined data for both "raw" and "combat"
combined_group_names = []
combined_correlation_values = []
combined_cov_type = []  # To distinguish between raw and combat data
param = "polarangle"
# Load the adjusted/combat harmonized data
for axis in axes: 
    fig, axs = plt.subplots(2, 2, figsize=(15,12))  # Adjust size as needed
    for hemi_index, hemi in enumerate(hemispheres):
        if hemi == "L":
            hemi_lower = "Left"
        elif hemi == "R":
            hemi_lower = "Right"
        for area_index, area in enumerate(areas):  
            combat_harmonized = pd.read_csv(f"{file_path}/raw_hemi-{hemi}_area-{area}_{param}_{axis}.csv", index_col=None).to_numpy()
            if area =="V2v":
                title = "ventral"
                area1="v"
                raw_harmonized = pd.read_csv(f"{file_path}/raw_hemi-{hemi}_area-V3v_{param}_{axis}.csv", index_col=None).to_numpy()
            elif area =="V2d":
                area1="d"
                raw_harmonized = pd.read_csv(f"{file_path}/raw_hemi-{hemi}_area-V3d_{param}_{axis}.csv", index_col=None).to_numpy()
                title = "dorsal"
            covars = pd.read_csv(f"{file_path}/covars_hemi-L.csv")
            # Load original model statistics:
            model_stats = pd.read_csv(f"/data/p_02915/SPOT/01_results_gradient_hemi-{hemi}_area-{area1}_{param}_S.csv")
            original_p_values = model_stats.iloc[:,4].values   # first column

            # Determine significance flag
            significant = original_p_values < 0.001

            # Initialize dictionaries to store averages for each group
            groups_mean = {
                'tri-2 \n prenatal': np.nanmean(combat_harmonized[:, covars[covars["group"] == "2nd"].index], axis=0),
                'tri-3 \n prenatal': np.nanmean(combat_harmonized[:, covars[covars["group"] == "3rd"].index], axis=0),
                'preterm \n neonatal': np.nanmean(combat_harmonized[:, covars[covars["group"] == "preterm"].index], axis=0),
                'full-term \n neonatal': np.nanmean(combat_harmonized[:, covars[covars["group"] == "fullterm"].index], axis=0),
                'adolescent': np.nanmean(combat_harmonized[:, covars[covars["group"] == "adolescent"].index], axis=0),
                'adult': np.nanmean(combat_harmonized[:, covars[covars["group"] == "adult"].index], axis=0)
            }

            # Initialize lists to store combined data
            combined_group_names = []
            combined_correlation_values = []
            combined_data_type = []

            # Flatten the group dictionary into lists for combat_harmonized
            for group_name, group_mean_values in groups_mean.items():
                combined_group_names.extend([group_name] * len(group_mean_values))  # Repeat group name for each value
                combined_correlation_values.extend(group_mean_values)  # Add all the mean values
                combined_data_type.extend(['V2'] * len(group_mean_values))  # Indicate this is combat data

            # Calculate raw_harmonized means for the same groups
            raw_groups_mean = {
                'tri-2 \n prenatal': np.nanmean(raw_harmonized[:, covars[covars["group"] == "2nd"].index], axis=0),
                'tri-3 \n prenatal': np.nanmean(raw_harmonized[:, covars[covars["group"] == "3rd"].index], axis=0),
                'preterm \n neonatal': np.nanmean(raw_harmonized[:, covars[covars["group"] == "preterm"].index], axis=0),
                'full-term \n neonatal': np.nanmean(raw_harmonized[:, covars[covars["group"] == "fullterm"].index], axis=0),
                'adolescent': np.nanmean(raw_harmonized[:, covars[covars["group"] == "adolescent"].index], axis=0),
                'adult': np.nanmean(raw_harmonized[:, covars[covars["group"] == "adult"].index], axis=0)
            }

            # Flatten the raw group dictionary into lists
            for group_name, group_mean_values in raw_groups_mean.items():
                combined_group_names.extend([group_name] * len(group_mean_values))  # Repeat group name for each value
                combined_correlation_values.extend(group_mean_values)  # Add all the mean values
                combined_data_type.extend(['V3'] * len(group_mean_values))  # Indicate this is raw data

            # Create a DataFrame for the adjusted/combat and raw data
            df_combined = pd.DataFrame({
                'Group': combined_group_names,
                'Direction': combined_correlation_values,
                'Data Type': combined_data_type  # New column to differentiate data types
            })
            df_combined = df_combined.fillna(0)

            print(df_combined)

            ax = axs[area_index, hemi_index]
            ax.set_title(f"{hemi_lower} {title}", fontsize=22)

            

            # Use hue to differentiate between Combat and Raw data
            pt.RainCloud(x="Group", y="Direction", data=df_combined, 
                        hue='Data Type', bw=.2, width_viol=.6, ax=ax, orient="h", 
                        alpha=.65, dodge=True, 
                        box_showfliers=False)

            if area_index == 0:
                ax.set_xlabel("Gradient magnitude", fontsize=20) 
            elif area_index == 1:
                ax.set_xlabel("Gradient magnitude", fontsize=20)  # X-axis label font size
            ax.set_ylabel(None)            
            if hemi_index == 0:                
                ax.tick_params(axis='x', which='major', labelsize=18)  # Major tick labels
                ax.tick_params(axis='y', which='major', labelsize=20)  # Minor tick labels
            elif hemi_index == 1:
                ax.set_yticks([])
                ax.tick_params(axis='x', which='major', labelsize=18)  # Major tick labels
            ax.set_xlim(-5, 5)  
            ax.legend_.remove()
            hue_order=['V2', "V3"]
            pairs=[
                (("tri-2 \n prenatal", "V2"), ("tri-2 \n prenatal", "V3")),
                (("tri-3 \n prenatal", "V2"), ("tri-3 \n prenatal", "V3")),
                (("preterm \n neonatal", "V2"), ("preterm \n neonatal", "V3")),
                (("full-term \n neonatal", "V2"), ("full-term \n neonatal", "V3")),
                (("adolescent", "V2"), ("adolescent", "V3")),
                (("adult", "V2"), ("adult", "V3")),
                ]
                # Loop through all groups
            unique_groups = df_combined["Group"].unique()

            for i, group in enumerate(unique_groups):
                if significant[i]:
                    ax.text(
                        0,                       # near right edge
                        i-0.5,                         # y-position at group index
                        "*",                       # mark
                        fontsize=15,
                        color="black",
                        va="center"
                    )


    
    plt.tight_layout()
    handles, labels = plt.gca().get_legend_handles_labels()
    # Create a legend and position it at the center of the plot
    legend = fig.legend(handles[0:2], labels[0:2], title="", loc='center', bbox_to_anchor=(0.05, 0.5), fontsize=18, frameon=False)
    # Set the alpha for the legend box
    legend.get_frame().set_alpha(0.5)  # Set the opacity to 50%
    fig.text(0.14, 0.52, "ventral", ha='center', va='center', fontsize=18)
    fig.text(0.51, 0.52, "dorsal", ha='center', va='center', fontsize=18)
    fig.text(0.59, 0.52, "ventral", ha='center', va='center', fontsize=18)
    fig.text(0.96, 0.52, "dorsal", ha='center', va='center', fontsize=18)
    fig.text(0.14, 0.025, "ventral", ha='center', va='center', fontsize=18)
    fig.text(0.51, 0.025, "dorsal", ha='center', va='center', fontsize=18)
    fig.text(0.59, 0.025, "ventral", ha='center', va='center', fontsize=18)
    fig.text(0.96, 0.025, "dorsal", ha='center', va='center', fontsize=18)



    # Set tick parameters
    
    #plt.tight_layout(rect=[0, 0, 1.08, 1]) 
    plt.tight_layout()
    #plt.show()
    fig.savefig(f"{file_path}/Figure3.png", dpi = 300)