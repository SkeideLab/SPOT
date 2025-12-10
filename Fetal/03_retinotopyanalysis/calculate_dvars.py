"""
Cacluate dvars
"""
import pandas as pd
import numpy as np
import nibabel as nib
from nilearn import surface

# Function to convert elements to float, setting non-convertible elements to np.nan


def to_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.nan


def flatten(arr):
    return arr.flatten()


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


LABELS_V1 = (1, 2)
LABELS_V2 = (3, 4)

for group in ["12-16y", "18-21y", "neonates<37", "neonates>37", "fetal<29", "fetal>29"]:
    results_list = []
    if group == "neonates<37":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/dhcp_surface/sub-{sub}/ses-{ses}/func/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "neonates>37":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/dhcp_surface/sub-{sub}/ses-{ses}/func/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal<29":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/func/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        )
        VISPARC_PATH = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-visualtopographywang2015_label-maxprob_dparc.label.gii")
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal>29":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/func/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        )
        VISPARC_PATH = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-visualtopographywang2015_label-maxprob_dparc.label.gii")
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "12-16y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/"
            "{sub}_hemi-{hemi}-Atlas-MSMAll_hp0_clean_bold.func.gii"
        )
        VISPARC_PATH = (
            "/data/p_02915/SPOT/00_HCP/template/hemi-{hemi}_mesh-native_dens-native_desc-visualtopographywang2015_label-maxprob_dparc.label.gii")
        subject_info = pd.read_csv(
            '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "18-21y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/"
            "{sub}_hemi-{hemi}-Atlas-MSMAll_hp0_clean_bold.func.gii"
        )
        VISPARC_PATH = (
            "/data/p_02915/SPOT/00_HCP/template/hemi-{hemi}_mesh-native_dens-native_desc-visualtopographywang2015_label-maxprob_dparc.label.gii")
        subject_info = pd.read_csv(
            '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
        sub_num = len(subject_info["sub_id"])

    group_name = f"{group}"
    parameters_L = []
    parameters_R = []
    # FORMAT PATHS FOR INPUT AND OUTPUT
    for index, row in subject_info.iterrows():
        if group == "12-16y" or group == "18-21y":
            for hemi in ["L", "R"]:
                visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
                indices_v2 = get_indices_roi(LABELS_V2, visparc)
                indices_v1 = get_indices_roi(LABELS_V1, visparc)
                sub_id = subject_info.at[index, "sub_id"]
                prefix_model = PREFIX_MODEL.format(
                    sub=sub_id,
                    hemi=hemi,
                )
                input_path = prefix_model

                if hemi == "L":
                    ccf = surface.load_surf_data(input_path).astype(np.float64)
                    ccf_v1 = ccf[indices_v1]
                    ccf_v2 = ccf[indices_v2]
                    left = np.vstack((ccf_v1, ccf_v2))

                    diff_data = np.diff(left, axis=1)

                    # Square the differences
                    squared_diff_data = diff_data ** 2
                    # Compute the voxel-wise variance of the time series
                    voxel_variance = np.var(diff_data, axis=1, ddof=1)

                    # To avoid division by zero, set zero variance to a small number (epsilon)
                    epsilon = 1e-8
                    voxel_variance[voxel_variance == 0] = epsilon

                    # Standardize the squared differences
                    standardized_squared_diff = squared_diff_data / \
                        voxel_variance[..., np.newaxis]

                    # Mean of standardized squared differences across voxels
                    mean_standardized_squared_diff = np.mean(
                        standardized_squared_diff, axis=0)

                    # Square root of mean standardized squared differences
                    sdvars = np.sqrt(mean_standardized_squared_diff)
                    # Save the DVARS values to a CSV file
                    output_path = f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub_id}_{hemi}_dvars.csv'
                    pd.DataFrame(sdvars, columns=['DVARS']).to_csv(
                        output_path, index=False)
                    parameters_L.append(np.mean(sdvars))
                elif hemi == "R":
                    ccf = surface.load_surf_data(input_path).astype(np.float64)
                    ccf_v1 = ccf[indices_v1]
                    ccf_v2 = ccf[indices_v2]
                    right = np.vstack((ccf_v1, ccf_v2))

                    diff_data = np.diff(right, axis=1)

                    # Square the differences
                    squared_diff_data = diff_data ** 2
                    # Compute the voxel-wise variance of the time series
                    voxel_variance = np.var(diff_data, axis=1, ddof=1)

                    # To avoid division by zero, set zero variance to a small number (epsilon)
                    epsilon = 1e-8
                    voxel_variance[voxel_variance == 0] = epsilon

                    # Standardize the squared differences
                    standardized_squared_diff = squared_diff_data / \
                        voxel_variance[..., np.newaxis]

                    # Mean of standardized squared differences across voxels
                    mean_standardized_squared_diff = np.mean(
                        standardized_squared_diff, axis=0)

                    # Square root of mean standardized squared differences
                    sdvars = np.sqrt(mean_standardized_squared_diff)
                    # Save the DVARS values to a CSV file
                    output_path = f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub_id}_{hemi}_dvars.csv'
                    pd.DataFrame(sdvars, columns=['DVARS']).to_csv(
                        output_path, index=False)
                    parameters_R.append(np.mean(sdvars))
        else:
            for hemi in ["L", "R"]:
                sub_id = subject_info.at[index, "sub_id"]
                sess_id = subject_info.at[index, "sess_id"]
                sub = sub_id.replace('sub-', '')
                ses = sess_id.replace('ses-', '')
                visparc = nib.load(VISPARC_PATH.format(
                    hemi=hemi, sub=sub, ses=ses))
                indices_v2 = get_indices_roi(LABELS_V2, visparc)
                indices_v1 = get_indices_roi(LABELS_V1, visparc)
                prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                input_path = prefix_model

                if hemi == "L":
                    ccf = surface.load_surf_data(input_path).astype(np.float64)
                    ccf_v1 = ccf[indices_v1]
                    ccf_v2 = ccf[indices_v2]
                    left = np.vstack((ccf_v1, ccf_v2))

                    diff_data = np.diff(left, axis=1)

                    # Square the differences
                    squared_diff_data = diff_data ** 2
                    # Compute the voxel-wise variance of the time series
                    voxel_variance = np.var(diff_data, axis=1, ddof=1)

                    # To avoid division by zero, set zero variance to a small number (epsilon)
                    epsilon = 1e-8
                    voxel_variance[voxel_variance == 0] = epsilon

                    # Standardize the squared differences
                    standardized_squared_diff = squared_diff_data / \
                        voxel_variance[..., np.newaxis]

                    # Mean of standardized squared differences across voxels
                    mean_standardized_squared_diff = np.mean(
                        standardized_squared_diff, axis=0)

                    # Square root of mean standardized squared differences
                    sdvars = np.sqrt(mean_standardized_squared_diff)
                    # Save the DVARS values to a CSV file
                    output_path = f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_dvars.csv'
                    pd.DataFrame(sdvars, columns=['DVARS']).to_csv(
                        output_path, index=False)
                    parameters_L.append(np.mean(sdvars))
                elif hemi == "R":
                    ccf = surface.load_surf_data(input_path).astype(np.float64)
                    ccf_v1 = ccf[indices_v1]
                    ccf_v2 = ccf[indices_v2]
                    right = np.vstack((ccf_v1, ccf_v2))
                    diff_data = np.diff(right, axis=1)

                    # Square the differences
                    squared_diff_data = diff_data ** 2
                    # Compute the voxel-wise variance of the time series
                    voxel_variance = np.var(diff_data, axis=1, ddof=1)

                    # To avoid division by zero, set zero variance to a small number (epsilon)
                    epsilon = 1e-8
                    voxel_variance[voxel_variance == 0] = epsilon

                    # Standardize the squared differences
                    standardized_squared_diff = squared_diff_data / \
                        voxel_variance[..., np.newaxis]

                    # Mean of standardized squared differences across voxels
                    mean_standardized_squared_diff = np.mean(
                        standardized_squared_diff, axis=0)

                    # Square root of mean standardized squared differences
                    sdvars = np.sqrt(mean_standardized_squared_diff)
                    # Save the DVARS values to a CSV file
                    output_path = f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_dvars.csv'
                    pd.DataFrame(sdvars, columns=['DVARS']).to_csv(
                        output_path, index=False)
                    parameters_R.append(np.mean(sdvars))
    df = {f"{group}_L": parameters_L, f"{group}_R": parameters_R}
# Check lengths of each array in 'groups_flatt'
for key, value in df.items():
    print(f"Length of {key}: {len(value)}")

# Determine the maximum length of the arrays
max_length = max(len(value) for value in df.values())

# Pad the arrays with None (or np.nan) to ensure all arrays have the same length
for key, value in df.items():
    if len(value) < max_length:
        # Using np.nan for padding
        padding = np.full(max_length - len(value), np.nan)
        df[key] = np.concatenate((value, padding))

results_df = pd.DataFrame(df)
results_df.to_csv(f'/data/p_02915/SPOT/02_averaged_dvars.csv')
