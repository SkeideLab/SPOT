import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from neuromaps.datasets import fetch_fsaverage
import matplotlib.pyplot as plt
from nilearn import plotting
import seaborn as sns
from matplotlib import cm
import matplotlib.gridspec as gridspec
import scipy.stats as stats
import os
 
# Get the full path of the script
file_path = os.path.dirname(os.path.abspath(__file__))

# Compute R^2 values
def r2_score(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)  # Residual sum of squares
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)  # Total sum of squares
    return 1 - (ss_res / ss_tot)

def r_squared_significance(R2, n, k):
    """
    Test the significance of R-squared using an F-test.
    
    Parameters:
    R2 : float - R-squared value from the regression
    n  : int   - Number of observations (sample size)
    k  : int   - Number of predictors (independent variables)
    
    Returns:
    F-statistic and p-value
    """
    if R2 == 1:  # To avoid division by zero error
        return float('inf'), 0.0  
    
    # Compute the F-statistic
    F_stat = (R2 / k) / ((1 - R2) / (n - k - 1))
    
    # Compute p-value
    p_value = 1 - stats.f.cdf(F_stat, dfn=k, dfd=n - k - 1)
    
    return F_stat, p_value

def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
     # Ensure labels_area is a list
    if not isinstance(labels_area, list):
        labels_area = [labels_area]
    
    # Collect indices for all labels in labels_area
    indices_area = np.concatenate([
        np.nonzero(visparc.agg_data() == label)[0]
        for label in labels_area
    ])

    return indices_area


VISPARC_PATH = (
    "{file_path}/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")

left_temp = nib.load(VISPARC_PATH.format(file_path=file_path,hemi="L"))
indices_left_v2v = get_indices_roi([3], left_temp)
indices_left_v2d = get_indices_roi([4], left_temp)
indices_left_v3v = get_indices_roi([5], left_temp)
indices_left_v3d = get_indices_roi([6], left_temp)
right_temp = nib.load(VISPARC_PATH.format(file_path=file_path, hemi="R"))
indices_right_v2v = get_indices_roi([3], right_temp)
indices_right_v2d = get_indices_roi([4], right_temp)
indices_right_v3v = get_indices_roi([5], right_temp)
indices_right_v3d = get_indices_roi([6], right_temp)

df = pd.DataFrame()
indices_left_v1 = get_indices_roi([1, 2, 3], left_temp)
indices_right_v1 = get_indices_roi([1, 2, 3], right_temp)

# Assuming `bin` is defined and `df` contains your data
groups = ["2nd", "3rd", "preterm",
          "fullterm", 'adolescent', 'adult']
text_label = {'preterm': 'preterm neonatal', 'fullterm': 'full-term neonatal',
              '2nd': 'tri-2nd prenatal', '3rd': 'tri-3rd prenatal', 'adolescent': 'adolescent', 'adult': 'adult'}
left_label = ["a", "b", "c", "d", "e", "f", "g"]
colors = {'preterm': 'limegreen', 'fullterm': 'mediumaquamarine', '2nd': 'gold',
          '3rd': 'yellowgreen', 'adolescent': 'cadetblue', 'adult': 'steelblue'}
vmax = 10

# Extract the label data
fig, axs = plt.subplots(4, 2, figsize=(12, 15), sharex=True)

for i, group in enumerate(groups):    
    if group == "preterm":
        average_left = f"{file_path}/Averaged_younger_n_L_label-sigma_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_younger_n_R_label-sigma_desc-real_roi-v2th00_metric.gii"
        eccentricity_left = f"{file_path}/Averaged_younger_n_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        eccentricity_right = f"{file_path}/Averaged_younger_n_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "fullterm":
        average_left = f"{file_path}/Averaged_older_n_L_label-sigma_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_older_n_R_label-sigma_desc-real_roi-v2th00_metric.gii"
        eccentricity_left = f"{file_path}/Averaged_older_n_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        eccentricity_right = f"{file_path}/Averaged_older_n_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "2nd":
        average_left = f"{file_path}/Averaged_younger_fetal_L_label-sigma_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_younger_fetal_R_label-sigma_desc-real_roi-v2th00_metric.gii"
        eccentricity_left = f"{file_path}/Averaged_younger_fetal_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        eccentricity_right = f"{file_path}/Averaged_younger_fetal_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "3rd":
        average_left = f"{file_path}/Averaged_older_fetal_L_label-sigma_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_older_fetal_R_label-sigma_desc-real_roi-v2th00_metric.gii"
        eccentricity_left = f"{file_path}/Averaged_older_fetal_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        eccentricity_right = f"{file_path}/Averaged_older_fetal_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "adolescent":
        average_left = f"{file_path}/Averaged_young_L_label-sigma_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_young_R_label-sigma_desc-real_roi-v2th00_metric.gii"
        eccentricity_left = f"{file_path}/Averaged_young_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        eccentricity_right = f"{file_path}/Averaged_young_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "adult":
        average_left = f"{file_path}/Averaged_old_L_label-sigma_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_old_R_label-sigma_desc-real_roi-v2th00_metric.gii"
        eccentricity_left = f"{file_path}/Averaged_old_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        eccentricity_right = f"{file_path}/Averaged_old_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"


    # Load sigma values
    left_sigma = surface.load_surf_data(average_left)
    right_sigma = surface.load_surf_data(average_right)
    left_eccentricity = surface.load_surf_data(eccentricity_left)
    right_eccentricity = surface.load_surf_data(eccentricity_right)
  

    # Define processing function
    def process_data(ax, indices, x_labels, y_ticks, ylim, hemi, area, group, average_x, average_y, color):
        ecc = average_x[indices]
        sig = average_y[indices]
        
        df = pd.DataFrame({"eccentricity": ecc, "sigma": sig}).sort_values(by="eccentricity").reset_index(drop=True)
        # Remove rows with NaN values in 'distance' or 'sigma' columns
        df = df.dropna(subset=["eccentricity", "sigma"])

        # Ensure that there are no invalid values (like inf) in the data
        df = df[np.isfinite(df["eccentricity"]) & np.isfinite(df["sigma"])]
        # Bin and compute stats
        num_bins = 21
        bins = np.linspace(0, 20, num_bins + 1)
        df["eccentricity_bin"] = pd.cut(df["eccentricity"], bins=bins, labels=0.5 * (bins[:-1] + bins[1:]))
        stats = df.groupby("eccentricity_bin")["sigma"].agg(["mean", "std"]).reset_index()
        stats["eccentricity_bin"] = stats["eccentricity_bin"].astype(float)        
        # Plot
        #ax.scatter(df["eccentricity"], df["sigma"], color="k", s=0.5, alpha=0.1, label="Raw Data")  # Scatter plot
        #ax.plot(stats["eccentricity_bin"], stats["mean"], label="Averaged sigma", color="k")
        ax.set_xlim(0, 20)
        ax.set_ylim(*ylim)
        if group == "adult":            
            ax.set_yticks(y_ticks) 
            ax.set_yticklabels(y_ticks, fontsize=15)
            ax.set_ylabel("sigma", labelpad=5, fontsize=15)
            #if area =="V3v":
            ax.set_xticks(x_labels)
            ax.set_xticklabels(x_labels, fontsize=15)
            ax.set_xlabel("eccentricity", labelpad=2, fontsize=15)
            #else:
            #    ax.set_xticks([])
            
        if area == "V2d":
            if hemi =="L":
                ax.set_title("Left", fontsize=20)
            elif hemi =="R":
                ax.set_title("Right", fontsize=20)

        x = np.linspace(0, 20, 21)

        linear_fit = np.polyfit(df["eccentricity"], df["sigma"], 1)
        linear_values = np.polyval(linear_fit, df["eccentricity"])
        r2_linear = r2_score(df["sigma"], linear_values)

        #print(linear_values, stats["mean"])
        # Calculate the linear fit values using the mean slope and y-shift
        recon_values = np.polyval(linear_fit, x)
        F_stat, p_val = r_squared_significance(r2_linear, len(linear_values), 1)
        print(f"F-statistic: {F_stat:.4f}")
        print(f"P-value: {p_val:.4f}")
        
        # Plot the linear fit on the corresponding subplot
        #ax.plot(stats["eccentricity_bin"], linear_values, color="green", linestyle="--")
        # Plot a black dashed line on top
        ax.plot(x, recon_values, color=color, linewidth=2, linestyle="-", label=group, zorder=1)
    
    
    process_data(axs[0, 0], indices_left_v2d, [0, 10, 20], [0, 0.05, 0.1, 0.15, 0.2], (0, 0.2), "L", "V2d", group, left_eccentricity, left_sigma, colors[group])
    process_data(axs[1, 0], indices_left_v2v, [0, 10, 20], [0, 0.05, 0.1, 0.15, 0.2], (0, 0.2), "L", "V2v", group, left_eccentricity, left_sigma, colors[group])
    process_data(axs[2, 0], indices_left_v3d, [0, 10, 20], [0, 0.05, 0.1, 0.15, 0.2], (0, 0.2),"L", "V3d", group, left_eccentricity, left_sigma, colors[group])
    process_data(axs[3, 0], indices_left_v3v, [0, 10, 20], [0, 0.05, 0.1, 0.15, 0.2], (0, 0.2), "L", "V3v", group, left_eccentricity, left_sigma, colors[group])
    process_data(axs[0, 1], indices_right_v2d, [0, 10, 20], [0, 0.05, 0.1, 0.15, 0.2], (0, 0.2), "R", "V2d", group, right_eccentricity, right_sigma, colors[group])
    process_data(axs[1, 1], indices_right_v2v, [0, 10, 20], [0, 0.05, 0.1, 0.15, 0.2], (0, 0.2), "R", "V2v", group, right_eccentricity, right_sigma, colors[group])
    process_data(axs[2, 1], indices_right_v3d, [0, 10, 20], [0, 0.05, 0.1, 0.15, 0.2], (0, 0.2), "R", "V3d", group, right_eccentricity, right_sigma, colors[group])
    process_data(axs[3, 1], indices_right_v3v, [0, 10, 20], [0, 0.05, 0.1, 0.15, 0.2], (0, 0.2), "R", "V3v", group, right_eccentricity, right_sigma, colors[group])

plt.tight_layout(rect=[0.1, 0.1, 1, 0.95])
fig.text(0.56, 0.945, "V2", ha='center', va='center', fontsize=24, weight='bold')
fig.text(0.56, 0.52, "V3", ha='center', va='center', fontsize=24, weight='bold')
fig.text(0.05, 0.85, "dorsal", ha='center', va='center', fontsize=20, weight='bold')
fig.text(0.05, 0.64, "ventral", ha='center', va='center', fontsize=20, weight='bold')
fig.text(0.05, 0.43, "dorsal", ha='center', va='center', fontsize=20, weight='bold')
fig.text(0.05, 0.215, "ventral", ha='center', va='center', fontsize=20, weight='bold')
# Add legend at the bottom
handles = [plt.Line2D([0], [0], color=colors[group], linewidth=2, linestyle="-", label=text_label[group]) for group in groups]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.56, 0.02), ncol=3, fontsize=18, frameon=False)

plt.subplots_adjust(hspace=0.4, wspace=0.3)
plt.show()
fig.savefig(f"{file_path}/FigureS4.png", dpi=300, bbox_inches='tight')
 