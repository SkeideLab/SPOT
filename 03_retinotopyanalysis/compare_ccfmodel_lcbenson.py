"""Compare the ccfmodel results in V2 to a local spatial correlations model and to the direct atlas values
"""

import numpy as np
from nilearn import surface

parameters_to_compare = ["eccentricity", "polar_angle"]


# COMPARING CCF PARAMS WITH BENSON TEMPLATE
def get_indices_roi(labels_area, visparc_array):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (numpy.array): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        np.logical_or(visparc_array == labels_area[0], visparc_array == labels_area[1])
    )[0]
    return indices_area


def calc_msd(param_area1, param_area2):
    # mean squared difference between two different model results in an area
    n1 = param_area1.size
    n2 = param_area2.size
    assert (
        n1 == n2
    ), "areas don't have the same number of vertices and cannot be compared!"

    return np.sum((param_area1 - param_area2) ** 2) / n1


root_dir = "/data/p_02495/dhcp_derivatives"
# file_ccf_v0 = "/data/p_02495/dhcp_derivatives/ccfmodel/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}_desc-test_v0i.gii"
file_eccentricitybenson = "/data/p_02495/dhcp_derivatives/dhcp_surface/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-eccentretinotbenson2014_seg.shape.gii"
file_anglebenson = "/data/p_02495/dhcp_derivatives/dhcp_surface/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-angleretinotbenson2014_seg.shape.gii"
file_visparc = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
file_eccccf = "/data/p_02495/dhcp_derivatives/ccfmodel/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}_label-eccentricity_desc-ccf_roi-v2_metric.gii"
file_angleccf = "/data/p_02495/dhcp_derivatives/ccfmodel/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}_label-angle_desc-ccf_roi-v2_metric.gii"
file_ecclc = "/data/p_02495/dhcp_derivatives/ccfmodel/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}_label-eccentricity_desc-lc_roi-v2_metric.gii"
file_anglelc = "/data/p_02495/dhcp_derivatives/ccfmodel/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}_label-angle_desc-lc_roi-v2_metric.gii"
sub = "CC00058XX09"
ses = "11300"

labels_v1 = [1, 2]  # see `ROIfiles_Labeling.txt` in this directory
labels_v2 = [3, 4]

for hemi in ["L"]:

    # load visual parcellation to define roi
    visparc = surface.load_surf_data(
        file_visparc.format(root_dir=root_dir, sub=sub, ses=ses, hemi=hemi)
    )
    indices_v1 = get_indices_roi(labels_v1, visparc)
    indices_v2 = get_indices_roi(labels_v2, visparc)

    # retinotopy data from benson template, ccf model and lc simulation
    ecc_benson = surface.load_surf_data(
        file_eccentricitybenson.format(sub=sub, ses=ses, hemi=hemi)
    )
    ecc_ccf = surface.load_surf_data(file_eccccf.format(sub=sub, ses=ses, hemi=hemi))
    ecc_lc = surface.load_surf_data(file_ecclc.format(sub=sub, ses=ses, hemi=hemi))
    angle_benson = surface.load_surf_data(
        file_anglebenson.format(sub=sub, ses=ses, hemi=hemi)
    )
    angle_ccf = surface.load_surf_data(
        file_angleccf.format(sub=sub, ses=ses, hemi=hemi)
    )
    angle_lc = surface.load_surf_data(file_anglelc.format(sub=sub, ses=ses, hemi=hemi))

    # calculate mean square differences between models for eccentricity and angle
    msd_ecc_bensonccf = calc_msd(ecc_benson[indices_v2], ecc_ccf[indices_v2])
    msd_angle_bensonccf = calc_msd(angle_benson[indices_v2], angle_ccf[indices_v2])

    msd_ecc_ccflc = calc_msd(ecc_lc[indices_v2], ecc_ccf[indices_v2])
    msd_angle_ccflc = calc_msd(angle_lc[indices_v2], angle_ccf[indices_v2])

    print(f"MSD for eccentricity: comparison ccf vs benson: {msd_ecc_bensonccf}")
    print(f"MSD for eccentricity: comparison ccf vs lc: {msd_ecc_ccflc}")
    print(f"MSD for polar angle: comparison ccf vs benson: {msd_angle_bensonccf}")
    print(f"MSD for polar angle: comparison ccf vs lc: {msd_angle_bensonccf}")
