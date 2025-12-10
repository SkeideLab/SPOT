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

def boferroni_correction(p_values):
    """Bonferroni correction."""
    corrected_p_values = p_values * 15  # Bonferroni adjustment
    corrected_p_values = np.minimum(corrected_p_values, 1)  # Cap at 1.0
    return corrected_p_values

def flatten(arr):
    return arr.flatten()


# Create figure with subplots (one for tSNR and one for sSNR)
fig, axs = plt.subplots(2, 1, figsize=(8, 16))  # Adjust size as needed
colors = {
    'tri-2 \n prenatal': 'limegreen',
    'tri-3 \n prenatal': 'mediumaquamarine',
    'preterm \n neonatal': 'gold',
    'fullterm \n neonatal': 'yellowgreen',
    'adolescent': 'cadetblue',
    'adult': 'steelblue'
}
# Initialize lists to hold combined data for both "raw" and "combat"
combined_group_names = []
combined_correlation_values = []
combined_cov_type = []  # To distinguish between raw and combat data

# Load the adjusted/combat harmonized data (for both tSNR and sSNR)
method = "raw"
for param in ["t", "s"]:  # Loop for both tSNR and sSNR
    combat_harmonized = pd.read_csv(f"{file_path}/{method}_{param}snr.csv", index_col=None).to_numpy()
    covars = pd.read_csv(f"{file_path}/covars_hemi-L.csv")
    groups_flatt={}
    # Group data into categories (flattened)
    groups_flatt = {
        'tri-2 \n prenatal': flatten(combat_harmonized[covars[covars["group"] == "2nd"].index]),
        'tri-3 \n prenatal': flatten(combat_harmonized[covars[covars["group"] == "3rd"].index]),
        'preterm \n neonatal': flatten(combat_harmonized[covars[covars["group"] == "preterm"].index]),
        'fullterm \n neonatal': flatten(combat_harmonized[covars[covars["group"] == "fullterm"].index]),
        'adolescent': flatten(combat_harmonized[covars[covars["group"] == "adolescent"].index]),
        'adult': flatten(combat_harmonized[covars[covars["group"] == "adult"].index])
    }
    combined_group_names = []
    combined_correlation_values =[]
    # Flatten the group dictionary into lists
    for group_name, values in groups_flatt.items():
        combined_group_names.extend([group_name] * len(values))  # Repeat group name for each value
        combined_correlation_values.extend(values)  # Add all the correlation values
    
    # Create a DataFrame for the adjusted/combat data
    df_combined = pd.DataFrame({
        'Group': combined_group_names,
        f'{param}SNR': combined_correlation_values,
    })
    #df_combined['Group'] = df_combined['Group'].astype('category')
    #df_combined[f'{param}SNR'] = pd.to_numeric(df_combined[f'{param}SNR'], errors='coerce')

    pt.RainCloud(x="Group", y=f'{param}SNR', data=df_combined,
                palette=colors, bw=.2, width_viol=.6, ax=axs[0 if param == 't' else 1], orient="h",
                alpha=.65, dodge=True, box_showfliers=False)
    # Adjusting plot limits and titles
    if param == 't':
        axs[0].set_xlim(0, 125)
    else:
        axs[1].set_xlim(0, 6)
    axs[0 if param == 't' else 1].set_xlabel(f'{param}SNR', fontsize=25)
    axs[0 if param == 't' else 1].tick_params(axis='both', which='major', labelsize=25)
    axs[0 if param == 't' else 1].tick_params(axis='both', which='minor', labelsize=18)

    # Load bootstrapped p-values for statistical comparison
    bootstrap_p_value = pd.read_csv(f"{file_path}/01_results_{param}snr_raw_v2.csv")
    bootstrap_p_value['Original p-value'] = pd.to_numeric(bootstrap_p_value['Original p-value'], errors='coerce')

    print(bootstrap_p_value.shape)
    # Bootstrapped p-values for group comparisons
    bootstrapped_p_values = {
        ("tri-2 \n prenatal", "tri-3 \n prenatal"): boferroni_correction(bootstrap_p_value.at[0, "Original p-value"]),
        ("tri-2 \n prenatal", "preterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[1, "Original p-value"]),
        ("tri-2 \n prenatal", "fullterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[2, "Original p-value"]),
        ("tri-2 \n prenatal", "adolescent"): boferroni_correction(bootstrap_p_value.at[3, "Original p-value"]),
        ("tri-2 \n prenatal", "adult"): boferroni_correction(bootstrap_p_value.at[4, "Original p-value"]),
        ("tri-3 \n prenatal", "preterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[5, "Original p-value"]),
        ("tri-3 \n prenatal", "fullterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[6, "Original p-value"]),
        ("tri-3 \n prenatal", "adolescent"): boferroni_correction(bootstrap_p_value.at[7, "Original p-value"]),
        ("tri-3 \n prenatal", "adult"): boferroni_correction(bootstrap_p_value.at[8, "Original p-value"]),
        ("preterm \n neonatal", "fullterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[9, "Original p-value"]),
        ("preterm \n neonatal", "adolescent"): boferroni_correction(bootstrap_p_value.at[10, "Original p-value"]),
        ("preterm \n neonatal", "adult"): boferroni_correction(bootstrap_p_value.at[11, "Original p-value"]),
        ("fullterm \n neonatal", "adolescent"): boferroni_correction(bootstrap_p_value.at[12, "Original p-value"]),
        ("fullterm \n neonatal", "adult"): boferroni_correction(bootstrap_p_value.at[13, "Original p-value"]),
        ("adolescent", "adult"): boferroni_correction(bootstrap_p_value.at[14, "Original p-value"]),
    }

    # Set significance threshold
    significance_threshold = 0.01
    significant_comparisons = [comp for comp, p in bootstrapped_p_values.items() if p < significance_threshold]
    significant_p_values = [p for comp, p in bootstrapped_p_values.items() if p < significance_threshold]

    # Annotator for statistical testing only for significant comparisons
    if significant_comparisons:
        annotator = Annotator(axs[0 if param == 't' else 1], significant_comparisons, 
                      data=df_combined, x=f'{param}SNR', y="Group", orient='h')
        annotator.configure(text_format='star', pvalue_thresholds=[(1, ""), (0.001, "*")], loc="outside")
        annotator.set_pvalues_and_annotate(significant_p_values)
    axs[0].set_ylabel(None)
    axs[1].set_ylabel(None)

 
# Adjust plot layout
#plt.tight_layout()  # Keep margins
fig.savefig(f"{file_path}/FigureS3.png", dpi=300, bbox_inches='tight')