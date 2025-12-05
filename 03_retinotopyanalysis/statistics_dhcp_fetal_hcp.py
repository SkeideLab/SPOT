import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import chisquare
from scipy.stats import bootstrap
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
# Get the screen dimensions using tkinter
import tkinter as tk

root = tk.Tk()
screen_width_px = root.winfo_screenwidth()
screen_height_px = root.winfo_screenheight()
root.destroy()

# Convert pixels to inches (assuming 100 dpi for conversion)
dpi = 100
screen_width_in = screen_width_px / dpi
screen_height_in = screen_height_px / dpi

# Set desired figure size to half of screen width and full screen height
fig_width_in = screen_width_in / 2
fig_height_in = screen_height_in


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


def chi_squared_test(dist1, dist2, bin):
    # Compute histogram for dist1
    counts1, bin_edges = np.histogram(dist1, bins=bin)
    counts2, _ = np.histogram(dist2, bins=bin)

    # Manually normalize the histograms
    # / len(dist1) + 0.01  # Adding a small value to avoid zero counts
    normalized_counts1 = counts1 + 5
    # / len(dist2) + 0.01  # Adding a small value to avoid zero counts
    normalized_counts2 = counts2 + 5

    # Combine the histograms into a contingency table
    contingency_table = np.array([normalized_counts1, normalized_counts2])

    # Calculate chi-squared statistic and p-value
    chi2, p, _, _ = chi2_contingency(contingency_table)
    return chi2, p


def boot_strap_chi(group1_data, group2_data, num_bootstraps, bin, group, group_1):
    num_bootstraps = num_bootstraps
    chi_squared_values = []
    # print(group1_data.shape)
    # print(group2_data.shape)

    for _ in range(num_bootstraps):
        # Resample subjects with replacement
        resampled_indices = np.random.choice(
            range(group1_data.shape[0]), size=group1_data.shape[0], replace=True)
        resampled_indices1 = np.random.choice(
            range(group2_data.shape[0]), size=group2_data.shape[0], replace=True)
        resampled_group1 = group1_data[resampled_indices]
        resampled_group2 = group2_data[resampled_indices1]
        # Flatten the resampled data
        resampled_group1_flat = resampled_group1[resampled_group1 != 0]
        resampled_group2_flat = resampled_group2[resampled_group2 != 0]

        # Calculate chi-squared statistic for the resampled data
        chi_squared, _ = chi_squared_test(
            resampled_group1_flat, resampled_group2_flat, bin)
        chi_squared_values.append(chi_squared)

    # Calculate confidence interval
    lower_percentile = (1 - 0.95) / 2 * 100
    upper_percentile = (1 + 0.95) / 2 * 100
    lower_bound = np.percentile(chi_squared_values, lower_percentile)
    upper_bound = np.percentile(chi_squared_values, upper_percentile)

    # Calculate p-value
    group1_data = group1_data[group1_data != 0]
    group2_data = group2_data[group2_data != 0]
    # print(group1_data.shape)
    # print(group2_data.shape)
    observed_chi_squared, p = chi_squared_test(group1_data, group2_data, bin)
    print(f"{group} vs. {group_1}")
    # print(format(p, '.50f'))
    p_value = np.mean(np.array(chi_squared_values) >= observed_chi_squared)

    print(f"Observed Chi-Square Statistic: {observed_chi_squared:.3f}")
    print(
        f"Bootstrap Confidence Interval: ({lower_bound:.3f}, {upper_bound:.3f})")
    # print("Bootstrapped p-value:", p_value)


VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")
LABELS_V2 = (3, 4)

bin = np.linspace(-0.7, 0.7, 29)
# bin = np.linspace(3, 25, 10)
print(bin)
df = pd.DataFrame()
for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]:
    if group == "neonates<37":
        print(group)
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-real_rss.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "neonates>37":
        print(group)
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-real_rss.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal<29":
        print(group)
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-real_rss.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal>29":
        print(group)
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-real_rss.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "12-16y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-real_rss.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "18-21y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-real_rss.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
        sub_num = len(subject_info["sub_id"])

    for hemi in ["L", "R"]:
        visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
        indices_v2 = get_indices_roi(LABELS_V2, visparc)
        group_name = f"{group}_{hemi}"
        print(group_name)
        parameters = []
    # FORMAT PATHS FOR INPUT AND OUTPUT
        for index, row in subject_info.iterrows():
            if group == "12-16y" or group == "18-21y":
                sub_id = subject_info.at[index, "sub_id"]
                prefix_model = PREFIX_MODEL.format(
                    sub=sub_id,
                    hemi=hemi,
                )
                input_path = prefix_model

                ccf = surface.load_surf_data(input_path)
                ccf_v0 = ccf[indices_v2].astype(np.float64)
                parameters.append(ccf_v0)
            else:
                sub_id = subject_info.at[index, "sub_id"]
                sess_id = subject_info.at[index, "sess_id"]
                sub = sub_id.replace('sub-', '')
                ses = sess_id.replace('ses-', '')
                prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                input_path = prefix_model

                ccf = surface.load_surf_data(input_path)
                ccf_v0 = ccf[indices_v2].astype(np.float64)
                parameters.append(ccf_v0)

        if group == "neonates<37" and hemi == "L":
            neonates_y_L = np.array(parameters)
        elif group == "neonates>37" and hemi == "L":
            neonates_o_L = np.array(parameters)
        elif group == "fetal<29" and hemi == "L":
            fetal_y_L = np.array(parameters)
        elif group == "fetal>29" and hemi == "L":
            fetal_o_L = np.array(parameters)
        elif group == "12-16y" and hemi == "L":
            hcp_y_L = np.array(parameters)
        elif group == "18-21y" and hemi == "L":
            hcp_o_L = np.array(parameters)
        elif group == "neonates<37" and hemi == "R":
            neonates_y_R = np.array(parameters)
        elif group == "neonates>37" and hemi == "R":
            neonates_o_R = np.array(parameters)
        elif group == "fetal<29" and hemi == "R":
            fetal_y_R = np.array(parameters)
        elif group == "fetal>29" and hemi == "R":
            fetal_o_R = np.array(parameters)
        elif group == "12-16y" and hemi == "R":
            hcp_y_R = np.array(parameters)
        elif group == "18-21y" and hemi == "R":
            hcp_o_R = np.array(parameters)

        data = np.array(parameters).reshape(-1)
        data_nozero = data[data != 0]
        data_nozero_flatten = data.flatten()
        data_r = pd.DataFrame({'r': data_nozero_flatten, 'Group': group_name})
        df = pd.concat([df, data_r], ignore_index=True)

counts1, _ = np.histogram(fetal_y_L, bins=bin)
counts2, _ = np.histogram(fetal_o_L, bins=bin)
counts3, _ = np.histogram(neonates_y_L, bins=bin)
counts4, _ = np.histogram(neonates_o_L, bins=bin)
counts5, _ = np.histogram(hcp_y_L, bins=bin)
counts6, _ = np.histogram(hcp_o_L, bins=bin)

# Manually normalize the histograms
normalized_counts1 = counts1 + 5
normalized_counts2 = counts2 + 5
normalized_counts3 = counts3 + 5
normalized_counts4 = counts4 + 5
normalized_counts5 = counts5 + 5
normalized_counts6 = counts6 + 5
table = [normalized_counts1, normalized_counts2, normalized_counts3,
         normalized_counts4, normalized_counts5, normalized_counts6]
contigency_table = np.array(table)
chi2, p, dof, expected = chi2_contingency(contigency_table)

print(np.array_equal(counts6, counts5))
print(np.array_equal(hcp_y_L, hcp_o_L))
print(np.array_equal(hcp_y_R, hcp_o_R))
# print(f"Chi2: {chi2}")
# print(f"P-value: {p}")
# print(f"Degrees of freedom: {dof}")
# print("Expected frequencies:")
# print(expected)

for hemi in ["left", "right"]:
    if hemi == "left":
        print(hemi)
        for group in ["fetal<29", "fetal>29", "neonates<37", "neonates>37",  "12-16y", "18-21y"]:
            for group_1 in ["fetal<29", "fetal>29", "neonates<37", "neonates>37",  "12-16y", "18-21y"]:
                if group == "fetal<29" and group_1 == "fetal>29":
                    group1_data = fetal_y_L
                    group2_data = fetal_o_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal<29" and group_1 == "neonates<37":
                    group1_data = fetal_y_L
                    group2_data = neonates_y_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal<29" and group_1 == "neonates>37":
                    group1_data = fetal_y_L
                    group2_data = neonates_o_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal<29" and group_1 == "12-16y":
                    group1_data = fetal_y_L
                    group2_data = hcp_y_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal<29" and group_1 == "18-21y":
                    group1_data = fetal_y_L
                    group2_data = hcp_o_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal>29" and group_1 == "neonates<37":
                    group1_data = fetal_o_L
                    group2_data = neonates_y_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal>29" and group_1 == "neonates>37":
                    group1_data = fetal_o_L
                    group2_data = neonates_o_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal>29" and group_1 == "12-16y":
                    group1_data = fetal_o_L
                    group2_data = hcp_y_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal>29" and group_1 == "18-21y":
                    group1_data = fetal_o_L
                    group2_data = hcp_o_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates<37" and group_1 == "neonates>37":
                    group1_data = neonates_y_L
                    group2_data = neonates_o_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates<37" and group_1 == "12-16y":
                    group1_data = neonates_y_L
                    group2_data = hcp_y_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates<37" and group_1 == "18-21y":
                    group1_data = neonates_y_L
                    group2_data = hcp_o_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates>37" and group_1 == "12-16y":
                    group1_data = neonates_o_L
                    group2_data = hcp_y_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates>37" and group_1 == "18-21y":
                    group1_data = neonates_o_L
                    group2_data = hcp_o_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "12-16y" and group_1 == "18-21y":
                    group1_data = hcp_y_L
                    group2_data = hcp_o_L
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
    elif hemi == "right":
        print(hemi)
        for group in ["fetal<29", "fetal>29", "neonates<37", "neonates>37",  "12-16y", "18-21y"]:
            for group_1 in ["fetal<29", "fetal>29", "neonates<37", "neonates>37",  "12-16y", "18-21y"]:
                if group == "fetal<29" and group_1 == "fetal>29":
                    group1_data = fetal_y_R
                    group2_data = fetal_o_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal<29" and group_1 == "neonates<37":
                    group1_data = fetal_y_R
                    group2_data = neonates_y_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal<29" and group_1 == "neonates>37":
                    group1_data = fetal_y_R
                    group2_data = neonates_o_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal<29" and group_1 == "12-16y":
                    group1_data = fetal_y_R
                    group2_data = hcp_y_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal<29" and group_1 == "18-21y":
                    group1_data = fetal_y_R
                    group2_data = hcp_o_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal>29" and group_1 == "neonates<37":
                    group1_data = fetal_o_R
                    group2_data = neonates_y_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal>29" and group_1 == "neonates>37":
                    group1_data = fetal_o_R
                    group2_data = neonates_o_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal>29" and group_1 == "12-16y":
                    group1_data = fetal_o_R
                    group2_data = hcp_y_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "fetal>29" and group_1 == "18-21y":
                    group1_data = fetal_o_R
                    group2_data = hcp_o_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates<37" and group_1 == "neonates>37":
                    group1_data = neonates_y_R
                    group2_data = neonates_o_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates<37" and group_1 == "12-16y":
                    group1_data = neonates_y_R
                    group2_data = hcp_y_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates<37" and group_1 == "18-21y":
                    group1_data = neonates_y_R
                    group2_data = hcp_o_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates>37" and group_1 == "12-16y":
                    group1_data = neonates_o_R
                    group2_data = hcp_y_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "neonates>37" and group_1 == "18-21y":
                    group1_data = neonates_o_R
                    group2_data = hcp_o_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
                elif group == "12-16y" and group_1 == "18-21y":
                    group1_data = hcp_y_R
                    group2_data = hcp_o_R
                    boot_strap_chi(group1_data, group2_data,
                                   10000, bin, group, group_1)
