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

def boferroni_correction(p):
    corrected_p_values = p * 15  # Bonferroni adjustment
    corrected_p_values = np.minimum(p, 1)  # Cap at 1.0
    return corrected_p_values

def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        np.logical_or(
            visparc.agg_data() == labels_area[0], visparc.agg_data(
            ) == labels_area[1]
        )
    )[0]
    return indices_area


def flatten(arr):
    return arr.flatten()


VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-maxprob_seg.label.gii")
LABELS_V2 = (2, 3)

# Create a figure with subplots

areas = ["V2", "V3", "V2_V3"]
hemispheres = ["L", "R"]

colors = ["gold", "yellowgreen", "limegreen", "mediumaquamarine", "cadetblue", "steelblue"]

fig, axs = plt.subplots(figsize=(16, 8))  # Adjust size as needed
# Initialize lists to hold combined data for both "raw" and "combat"
combined_group_names = []
combined_correlation_values = []
combined_cov_type = []  # To distinguish between raw and combat data

# Load the adjusted/combat harmonized data
method = "raw"
combat_harmonized = pd.read_csv(f"/data/p_02915/SPOT/{method}_tnsr.csv", index_col=None).to_numpy()
covars = pd.read_csv(f"/data/p_02915/SPOT/covars_hemi-L.csv")

# Initialize lists to store combined data for adjusted data
# Group data into categories
groups_flatt = {
    'tri-2 \n prenatal': combat_harmonized[:, covars[covars["group"] == "2nd"].index].flatten(),
    'tri-3 \n prenatal': combat_harmonized[:, covars[covars["group"] == "3rd"].index].flatten(),
    'preterm \n neonatal': combat_harmonized[:, covars[covars["group"] == "preterm"].index].flatten(),
    'fullterm \n neonatal': combat_harmonized[:, covars[covars["group"] == "fullterm"].index].flatten(),
    'adolescent': combat_harmonized[:, covars[covars["group"] == "adolescent"].index].flatten(),
    'adult': combat_harmonized[:, covars[covars["group"] == "adult"].index].flatten()
}

# Flatten the group dictionary into lists
for group_name, values in groups_flatt.items():
    combined_group_names.extend([group_name] * len(values))  # Repeat group name for each value
    combined_correlation_values.extend(values)  # Add all the correlation values
    combined_cov_type.extend(['adjusted'] * len(values))  # Only adjusted data

# Create a DataFrame for the adjusted/combat data
df_combined = pd.DataFrame({
    'Group': combined_group_names,
    'tSNR': combined_correlation_values,
    'Covariate Type': combined_cov_type  # Only "adjusted"
})

# Generate horizontal raincloud plot for adjusted data
pt.RainCloud(x="Group", y='tSNR', data=df_combined,
             palette=colors, bw=.2, width_viol=.6, ax=axs, orient="h", 
             alpha=.65, dodge=True,
            box_showfliers=False)

axs.set_xlim(0,125)
# Load the bootstrapped p-values
bootstrap_p_value = pd.read_csv("/data/p_02915/SPOT/Result/01_results_tsnr_raw.csv")

# Bootstrapped p-values for group comparisons (adjusting this part to work for adjusted data)
bootstrapped_p_values = {
    ("tri-2 \n prenatal", "tri-3 \n prenatal"): boferroni_correction(bootstrap_p_value.at[1, "Original p-value"]),
    ("tri-2 \n prenatal", "preterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[2, "Original p-value"]),
    ("tri-2 \n prenatal", "fullterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[3, "Original p-value"]),
    ("tri-2 \n prenatal", "adolescent"): boferroni_correction(bootstrap_p_value.at[4, "Original p-value"]),
    ("tri-2 \n prenatal", "adult"): boferroni_correction(bootstrap_p_value.at[5, "Original p-value"]),
    ("tri-3 \n prenatal", "preterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[6, "Original p-value"]),
    ("tri-3 \n prenatal", "fullterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[7, "Original p-value"]),
    ("tri-3 \n prenatal", "adolescent"): boferroni_correction(bootstrap_p_value.at[8, "Original p-value"]),
    ("tri-3 \n prenatal", "adult"): boferroni_correction(bootstrap_p_value.at[9, "Original p-value"]),
    ("preterm \n neonatal", "fullterm \n neonatal"): boferroni_correction(bootstrap_p_value.at[10, "Original p-value"]),
    ("preterm \n neonatal", "adolescent"): boferroni_correction(bootstrap_p_value.at[11, "Original p-value"]),
    ("preterm \n neonatal", "adult"): boferroni_correction(bootstrap_p_value.at[12, "Original p-value"]),
    ("fullterm \n neonatal", "adolescent"): boferroni_correction(bootstrap_p_value.at[13, "Original p-value"]),
    ("fullterm \n neonatal", "adult"): boferroni_correction(bootstrap_p_value.at[14, "Original p-value"]),
    ("adolescent", "adult"): boferroni_correction(bootstrap_p_value.at[15, "Original p-value"]),
}

# Set the significance threshold
significance_threshold = 0.01

# Filter comparisons that are significant
significant_comparisons = [comp for comp, p in bootstrapped_p_values.items() if p < significance_threshold]
significant_p_values = [p for comp, p in bootstrapped_p_values.items() if p < significance_threshold]

# Initialize the annotator for statistical testing only for significant comparisons
if significant_comparisons:
    # Use x="tSNR" and y="Group" since it's horizontal
    annotator = Annotator(axs, significant_comparisons, data=df_combined, x="tSNR", y="Group", orient='h')

    # Configure the annotator (no need for perform_stat_test or custom_p_values here)
    annotator.configure(text_format='star', pvalue_thresholds=[(1, ""), (0.001, "*")], loc="outside")

    # Use the bootstrapped p-values when applying annotations
    annotator.set_pvalues_and_annotate(significant_p_values)

handles, labels = plt.gca().get_legend_handles_labels()
#plt.legend(handles[0:2], labels[0:2], title="")
plt.gca().set_ylabel(None)
axs.set_xlabel("tSNR", fontsize=25)  # X-axis label font size

# Set tick parameters
axs.tick_params(axis='both', which='major', labelsize=25)  # Major tick labels
axs.tick_params(axis='both', which='minor', labelsize=18)  # Minor tick labels

plt.tight_layout(rect=[0, 0, 1.08, 1]) 
#plt.tight_layout()
#plt.show()
fig.savefig("/data/p_02915/dhcp_derivatives_SPOT/Figures_v3/tSNR_v2.png")