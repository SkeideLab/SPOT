import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from neuromaps.datasets import fetch_fsaverage
import matplotlib.pyplot as plt
from scipy import stats


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
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii")
LABELS_V2 = [2, 3]


left_temp = nib.load(VISPARC_PATH.format(hemi="L"))
indices_left = get_indices_roi([2, 3], left_temp)
indices_left_v2 = get_indices_roi([2], left_temp)
indices_left_v3 = get_indices_roi([3], left_temp)
right_temp = nib.load(VISPARC_PATH.format(hemi="R"))
indices_right = get_indices_roi([2, 3], right_temp)
indices_right_v2 = get_indices_roi([2], right_temp)
indices_right_v3 = get_indices_roi([3], right_temp)


surfaces = fetch_fsaverage(density="164k")
lh, rh = surfaces["inflated"]
# Load sulcal depth data
# Load the left hemisphere sulcal depth data
lh_sulc = surface.load_surf_data(
    "/data/p_02915/templates/template_fsaverage/fsaverage/surf/lh.curv")
# Load the right hemisphere sulcal depth data
rh_sulc = surface.load_surf_data(
    "/data/p_02915/templates/template_fsaverage/fsaverage/surf/rh.curv")
lh_sulc_map_binary = np.where(lh_sulc < 0.0, 0.25, 0.6)
rh_sulc_map_binary = np.where(rh_sulc < 0.0, 0.25, 0.6)# Load the surface files for left and right hemispheres
indices_left_v1 = get_indices_roi([1, 2, 3], left_temp)
indices_right_v1 = get_indices_roi([1, 2, 3], right_temp)

groups = ["2nd", "3rd", "preterm",
          "fullterm", 'adolescent', 'adult', "benson"]
text_label = {'preterm': 'preterm neonatal', 'fullterm': 'full-term neonatal',
              '2nd': 'tri-2nd prenatal', '3rd': 'tri-3rd prenatal', 'adolescent': 'adolescent', 'adult': 'adult', "benson": "template"}
left_label = ["a", "b", "c", "d", "e", "f", "g"]
colors = {'preterm': 'limegreen', 'fullterm': 'mediumaquamarine', '2nd': 'gold',
          '3rd': 'yellowgreen', 'adolescent': 'cadetblue', 'adult': 'steelblue', "benson":"#656364"}
vmax = 180
# Load the GIFTI file (atlas file)
atlas_lh = nib.load('/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-L_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii')
atlas_rh = nib.load('/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-R_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii')
template_lh = nib.load(lh)  # Load the left hemisphere surface
template_rh = nib.load(rh)

# Extract the label data
atlas_data_l = atlas_lh.darrays[0].data
atlas_data_r = atlas_rh.darrays[0].data  # Extract the data array
# Specify the labels of interest
left_coords, left_faces = nib.freesurfer.io.read_geometry(f"/data/p_02915/templates/template_fsaverage/fsaverage/surf/lh.sphere")
right_coords, right_faces = nib.freesurfer.io.read_geometry(f"/data/p_02915/templates/template_fsaverage/fsaverage/surf/rh.sphere")
def divide_indices_by_z(coords, indices):
    z_values = coords[indices][:, 2]
    middle_z = np.median(z_values)
    below_indices = indices[z_values < middle_z]
    above_indices = indices[z_values >= middle_z]
    return below_indices, above_indices

# For left hemisphere V2
below_middle_indices_left_v2, above_middle_indices_left_v2 = divide_indices_by_z(left_coords, indices_left_v2)

# For left hemisphere V3
below_middle_indices_left_v3, above_middle_indices_left_v3 = divide_indices_by_z(left_coords, indices_left_v3)

# For right hemisphere V2
below_middle_indices_right_v2, above_middle_indices_right_v2 = divide_indices_by_z(right_coords, indices_right_v2)

# For right hemisphere V3
below_middle_indices_right_v3, above_middle_indices_right_v3 = divide_indices_by_z(right_coords, indices_right_v3)


distance_left= pd.read_csv(f"/data/p_02915/dhcp_derivatives_SPOT/L_polarangle_distance.csv")   
distance_right= pd.read_csv(f"/data/p_02915/dhcp_derivatives_SPOT/R_polarangle_distance.csv") 

# Compute R^2 values
def r2_score(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)  # Residual sum of squares
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)  # Total sum of squares
    return 1 - (ss_res / ss_tot)
fit = pd.DataFrame()

# Initialize empty lists
sub_id_list, group_list, hemi_list, area_list, r2_linear_list, slope_list, y_shift_list= [], [], [], [], [], [], []

for group in ["preterm", "fullterm", "2nd", "3rd", "adolescent", "adult", "benson"]:
    if group == "preterm":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_v2.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fullterm":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_v2.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "2nd":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "3rd":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "adolescent":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "adult":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "benson":
        PREFIX_MODEL=("/data/p_02915/SPOT/01_dataprep/retinotopy/"
                      "templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-angleretinotbenson2014_seg.shape.gii")

    for hemi in ["L", "R"]:
        visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
        indices_v2 = get_indices_roi(LABELS_V2, visparc)
        group_name = f"{group}_{hemi}"
        if hemi == "L":
            above_indices = [above_middle_indices_left_v2, above_middle_indices_left_v3]
            below_indices = [below_middle_indices_left_v2, below_middle_indices_left_v3]
            distance_data = distance_left
        else:
            above_indices = [above_middle_indices_right_v2, above_middle_indices_right_v3]
            below_indices = [below_middle_indices_right_v2, below_middle_indices_right_v3]
            distance_data = distance_right
        
        for j, area, indices_set in zip(range(2), ["dorsal", "ventral"], [above_indices, below_indices]):
            distance = np.concatenate([distance_data.loc[distance_data["Vertex"].isin(indices), "Distance_to_Border"] for indices in indices_set])
            indices = np.concatenate(indices_set)
            
            subject_iterator = subject_info.iterrows() if group != "benson" else [(0, None)]
            
            for index, row in subject_iterator:
                sub_id = subject_info.at[index, "sub_id"] if group != "benson" else 0
                prefix_model = PREFIX_MODEL.format(
                    sub=sub_id.replace('sub-', '') if group not in ["benson", "adolescent", "adult"] else sub_id,
                    ses=subject_info.at[index, "sess_id"].replace('ses-', '') if group not in ["benson", "adolescent", "adult"] else "",
                    hemi=hemi,
                )
                                
                polarangle = surface.load_surf_data(prefix_model)[indices]

                df = pd.DataFrame({"distance": distance, "polarangle": polarangle}).sort_values(by="distance").reset_index(drop=True)
                # Filter data to include only distances between -15 and 15
                df_filtered = df[(df["distance"] >= -15) & (df["distance"] <= 15)]
                bins = np.linspace(df["distance"].min(), df["distance"].max(), 31)
                df["distance_bin"] = pd.cut(df["distance"], bins=bins, labels=0.5 * (bins[:-1] + bins[1:]))
                # Step 3: Compute mean and std for each bin
                #stats = df.groupby("distance_bin")["polarangle"].agg(["mean", "std"]).reset_index()
                #stats["distance_bin"] = stats["distance_bin"].astype(float)  # Convert bin labels to float

                #linear_fit = np.polyfit(df["distance"], df["polarangle"], 1)
                #quadratic_fit = np.polyfit(df["distance"]**2, df["polarangle"], 1)
                linear_fit = np.polyfit(abs(df_filtered["distance"]), df_filtered["polarangle"], 1)
                linear_values = np.polyval(linear_fit, abs(df_filtered["distance"]))
                #linear_fit = np.polyfit(abs(stats["distance_bin"]), stats["mean"], 1)
                #linear_values = np.polyval(linear_fit, abs(stats["distance_bin"]))
                
                #r2_linear = r2_score(df["polarangle"], np.polyval(linear_fit, df["distance"]))
                #r2_quadratic = r2_score(df["polarangle"], quadratic_fit[0] * df["distance"]**2 + quadratic_fit[1])
                
                r2_linear = r2_score(df_filtered["polarangle"], linear_values)
                
                sub_id_list.append(sub_id)
                group_list.append(group)
                hemi_list.append(hemi)
                area_list.append(area)
                r2_linear_list.append(r2_linear)
                slope_list.append(linear_fit[0])
                y_shift_list.append(linear_fit[1])
                
fit = pd.DataFrame({
    "sub_id": sub_id_list,
    "Group": group_list,
    "hemi": hemi_list,
    "area": area_list,
    "linear_fit": r2_linear_list,
    "slope": slope_list,
    "y_shift": y_shift_list,
})
fit.to_csv('/data/p_02915/SPOT/fitting_linear_quad.csv')

# Step 4: Perform statistical tests for each group
alpha = 0.01

print("\n--- Statistical Tests Per Group, Hemisphere, and Area ---")
for (group, hemi, area), subset in fit.groupby(["Group", "hemi", "area"]):
    group_slopes = subset["slope"].values

    if len(group_slopes) < 3:
        print(f"\nGroup: {group}, Hemisphere: {hemi}, Area: {area} - Not enough data for statistical testing.")
        continue

    # Step 4.1: Normality test
    shapiro_test = stats.shapiro(group_slopes)
    print(f"\nGroup: {group}, Hemisphere: {hemi}, Area: {area}")
    print(f"Shapiro-Wilk test statistic: {shapiro_test.statistic}, p-value: {shapiro_test.pvalue}")

    # Step 4.2: Choose the alternative hypothesis based on area
    alternative_hypothesis = "greater" if area.lower() == "dorsal" else "less"

    # Step 4.3: Choose the appropriate test
    if shapiro_test.pvalue > alpha:  # Data is normal
        print(f"Data is normally distributed. Running one-sample t-test (alternative='{alternative_hypothesis}')...")
        t_stat, p_value = stats.ttest_1samp(group_slopes, 0, alternative=alternative_hypothesis)
        print(f"T-test statistic: {t_stat}, p-value: {p_value}")
    else:  # Data is not normal
        print(f"Data is not normally distributed. Running Wilcoxon signed-rank test (alternative='{alternative_hypothesis}')...")
        w_stat, p_value = stats.wilcoxon(group_slopes, alternative=alternative_hypothesis)
        print(f"Wilcoxon test statistic: {w_stat}, p-value: {p_value}")

    # Step 4.4: Interpretation
    if p_value < alpha:
        print(f"Conclusion for {group}, {hemi}, {area}: Reject H0. The slopes are significantly {'positive' if alternative_hypothesis == 'greater' else 'negative'}.")
    else:
        print(f"Conclusion for {group}, {hemi}, {area}: Fail to reject H0. No significant evidence that slopes are {'positive' if alternative_hypothesis == 'greater' else 'negative'}.")

# Prepare the plot
fig, ax = plt.subplots(2, 2)
x = np.linspace(-15, 15, 31)

for i, hemi in enumerate(["L", "R"]):  # Left and Right hemispheres
    for j, area in enumerate(["dorsal", "ventral"]):  # Dorsal and Ventral areas
        ax[j, i].set_title(f'{hemi} {area}')
        for group in ["2nd", "3rd", "preterm", "fullterm", "adolescent", "adult", "benson"]:
            # Filter the data for the current group, hemi, and area
            subset = fit.loc[(fit["Group"] == group) & (fit["hemi"] == hemi) & (fit["area"] == area)].reset_index(drop=True)
            
            if group == "benson":                
                slope = subset["slope"].values[0]
                y = subset["y_shift"].values[0]
            else:            
                # Perform bootstrapping for the slope and y-shift
                bootstrapped_slope = []
                bootstrapped_y = []
                
                for _ in range(10000):  # Bootstrapping iterations
                    slope_sample = np.random.choice(subset["slope"], size=len(subset["slope"]), replace=True)
                    bootstrapped_slope.append(np.mean(slope_sample))
                    y_sample = np.random.choice(subset["y_shift"], size=len(subset["y_shift"]), replace=True)
                    bootstrapped_y.append(np.mean(y_sample))
                
                # Calculate the mean of the bootstrapped samples
                slope = np.mean(bootstrapped_slope)
                y = np.mean(bootstrapped_y)
            
            # Calculate the linear fit values using the mean slope and y-shift
            linear_values = np.polyval([slope, y], abs(x))
            
            # Plot the linear fit on the corresponding subplot
            ax[j, i].plot(x, linear_values, label=f"{text_label[group]}, a = {slope:.3f}", color=colors[group])

plt.subplots_adjust(right=0.8, wspace = 0.6) 
# Display the legend and show the plot
fig.suptitle("y = a|x| + b")
for i in range(2):
    for j in range(2):
        ax[j, i].legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        ax[j, i].set_ylabel("polar angle", rotation=90, labelpad=5, fontsize=10, ha="right")
        ax[j, i].set_xlabel("distance (mm)", labelpad=5, fontsize=10)
plt.show()

