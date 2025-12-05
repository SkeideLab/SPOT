'''
Statistical comparison of gradients among six groups
'''
import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
from multiprocessing import Pool, cpu_count
import math


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


def get_indices_roi_v1(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        visparc.agg_data() == labels_area
    )[0]
    return indices_area


def flatten(arr):
    return arr.flatten()


VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")
LABELS_V2 = (3, 4, 5, 6)
original_mean = []
bootstrap_mean_all = []
bootstrap_p_value = []
bootstrap_ci_lows = []
bootstrap_ci_highs = []
index_label = []
df = pd.DataFrame()
for param in ["eccentricity", "polarangle"]:
    results_list = []
    print(param)
    for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]:
        if group == "neonates<37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )

            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_v2.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "neonates>37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )

            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_v2.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "fetal<29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )

            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "fetal>29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )

            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "12-16y":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{hemi}_{param}_gradient.gii"
            )

            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "18-21y":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{hemi}_{param}_gradient.gii"
            )

            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
            sub_num = len(subject_info["sub_id"])

        for hemi in ["L", "R"]:
            visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
            indices_v2 = get_indices_roi(LABELS_V2, visparc)
            indices_v2v = get_indices_roi_v1(LABELS_V2[0], visparc)
            indices_v2d = get_indices_roi_v1(LABELS_V2[1], visparc)
            indices_v3v = get_indices_roi_v1(LABELS_V2[2], visparc)
            indices_v3d = get_indices_roi_v1(LABELS_V2[3], visparc)
            group_name = f"{group}_{hemi}"
            parameters = []
            parameters_v2_ventral = []
            parameters_v2_dorsal = []
            parameters_v3_ventral = []
            parameters_v3_dorsal = []
        # FORMAT PATHS FOR INPUT AND OUTPUT
            for index, row in subject_info.iterrows():
                if group == "12-16y" or group == "18-21y":
                    sub_id = subject_info.at[index, "sub_id"]
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub_id,
                        hemi=hemi,
                        param=param,
                    )
                    input_path = prefix_model

                    ccf = surface.load_surf_data(input_path)
                    ccf_v0 = ccf[:, indices_v2].astype(np.float64)
                    ccf_v1 = ccf[:, indices_v2v].astype(np.float64)
                    ccf_v2 = ccf[:, indices_v2d].astype(np.float64)
                    ccf_v3 = ccf[:, indices_v3v].astype(np.float64)
                    ccf_v4 = ccf[:, indices_v3d].astype(np.float64)
                    parameters.append(ccf_v0)
                    parameters_v2_ventral.append(ccf_v1)
                    parameters_v2_dorsal.append(ccf_v2)
                    parameters_v3_ventral.append(ccf_v3)
                    parameters_v3_dorsal.append(ccf_v4)
                else:
                    sub_id = subject_info.at[index, "sub_id"]
                    sess_id = subject_info.at[index, "sess_id"]
                    sub = sub_id.replace('sub-', '')
                    ses = sess_id.replace('ses-', '')
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                        param=param,
                    )
                    input_path = prefix_model

                    ccf = surface.load_surf_data(input_path)
                    ccf_v0 = ccf[:, indices_v2].astype(np.float64)
                    ccf_v1 = ccf[:, indices_v2v].astype(np.float64)
                    ccf_v2 = ccf[:, indices_v2d].astype(np.float64)
                    ccf_v3 = ccf[:, indices_v3v].astype(np.float64)
                    ccf_v4 = ccf[:, indices_v3d].astype(np.float64)
                    parameters.append(ccf_v0)
                    parameters_v2_ventral.append(ccf_v1)
                    parameters_v2_dorsal.append(ccf_v2)
                    parameters_v3_ventral.append(ccf_v3)
                    parameters_v3_dorsal.append(ccf_v4)
            n_bootstraps = 10000

            R_v2v = []
            A_v2v = []
            S_v2v = []
            R_v2d = []
            A_v2d = []
            S_v2d = []

            R_v3v = []
            A_v3v = []
            S_v3v = []
            R_v3d = []
            A_v3d = []
            S_v3d = []
            for j in range(len(parameters_v2_ventral)):
                R_v2v.append(parameters_v2_ventral[j][0, :])
                A_v2v.append(parameters_v2_ventral[j][1, :])
                S_v2v.append(parameters_v2_ventral[j][2, :])
                R_v2d.append(parameters_v2_dorsal[j][0, :])
                A_v2d.append(parameters_v2_dorsal[j][1, :])
                S_v2d.append(parameters_v2_dorsal[j][2, :])

                R_v3v.append(parameters_v3_ventral[j][0, :])
                A_v3v.append(parameters_v3_ventral[j][1, :])
                S_v3v.append(parameters_v3_ventral[j][2, :])
                R_v3d.append(parameters_v3_dorsal[j][0, :])
                A_v3d.append(parameters_v3_dorsal[j][1, :])
                S_v3d.append(parameters_v3_dorsal[j][2, :])
                

            for testdata in [R_v2v, A_v2v, S_v2v, R_v2d, A_v2d, S_v2d, R_v3v, A_v3v, S_v3v, R_v3d, A_v3d, S_v3d]:
                if testdata is R_v2v:
                    test_name = "R_v2v"
                elif testdata is A_v2v:
                    test_name = "A_v2v"
                elif testdata is S_v2v:
                    test_name = "S_v2v"
                elif testdata is R_v2d:
                    test_name = "R_v2d"
                elif testdata is A_v2d:
                    test_name = "A_v2d"
                elif testdata is S_v2d:
                    test_name = "S_v2d"
                elif testdata is R_v3v:
                    test_name = "R_v3v"
                elif testdata is A_v3v:
                    test_name = "A_v3v"
                elif testdata is S_v3v:
                    test_name = "S_v3v"
                elif testdata is R_v3d:
                    test_name = "R_v3d"
                elif testdata is A_v3d:
                    test_name = "A_v3d"
                elif testdata is S_v3d:
                    test_name = "S_v3d"
                bootstrap_mean_list = []
                for _ in range(n_bootstraps):
                    # Resample the combined data
                    resampled_indices = np.random.choice(
                        len(testdata), size=len(testdata), replace=True)
                    # Convert testdata to a NumPy array if it isn't already one
                    if not isinstance(testdata, np.ndarray):
                        testdata = np.array(testdata)
                    resampled_data = testdata[resampled_indices]
                    # Compute the Kruskal-Wallis statistic for the resampled groups
                    bootstrap_mean = np.nanmean(resampled_data)
                    bootstrap_mean_list.append(bootstrap_mean)

                # Calculate the bootstrap p-value for Kruskal-Wallis test
                bootstrap_statistics = np.array(bootstrap_mean_list)
                # bootstrap_mean_p_value = np.sum(bootstrap_statistics-np.nanmean(testdata) >= np.nanmean(testdata)) / n_bootstraps
                test_statistic = np.abs(
                    np.mean(bootstrap_statistics) - np.nanmean(testdata))
                bootstrap_mean_p_value = 2 * min(np.mean(bootstrap_statistics >= np.nanmean(
                    testdata)), np.mean(bootstrap_statistics <= np.nanmean(testdata)))

                print("\nOriginal mean:", np.nanmean(testdata))
                print("Bootstrapped p-value:", bootstrap_mean_p_value)
                bootstrap_ci_low = np.percentile(bootstrap_statistics, 2.5)
                bootstrap_ci_high = np.percentile(bootstrap_statistics, 97.5)

                print(
                    f"Bootstrap CI for mean value - Low: {bootstrap_ci_low}, High: {bootstrap_ci_high}")
                original_mean.append(np.nanmean(testdata))
                bootstrap_mean_all.append(np.nanmean(bootstrap_statistics))
                bootstrap_p_value.append(bootstrap_mean_p_value)
                bootstrap_ci_lows.append(bootstrap_ci_low)
                bootstrap_ci_highs.append(bootstrap_ci_high)
                index_label.append(f"{group}_{hemi}_{param}_{test_name}")


results_dict = {
    "Original mean": original_mean,
    "bootstrap mean": bootstrap_mean_all,
    'Bootstrapped p-value': bootstrap_p_value,
    'Bootstrap CI low': bootstrap_ci_lows,
    'Bootstrap CI high': bootstrap_ci_highs
}

# Create a DataFrame from the dictionary
results_df = pd.DataFrame(results_dict, index=index_label)

# Save the DataFrame to a CSV file
results_df.to_csv(f'/data/p_02915/SPOT/01_results_gradient.csv')
