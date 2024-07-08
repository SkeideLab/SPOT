import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from neuromaps.datasets import fetch_fsaverage
import matplotlib.pyplot as plt
from nilearn import plotting
import seaborn as sns
from matplotlib import cm


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


VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")
LABELS_V2 = (3, 4)

df = pd.DataFrame()
for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]:
    if group == "neonates<37":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "neonates>37":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal<29":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal>29":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "12-16y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "18-21y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
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
        data_flatten = data.flatten()
        data_r = pd.DataFrame({'ecc': data_flatten, 'Group': group_name})
        df = pd.concat([df, data_r], ignore_index=True)

bin = np.linspace(0, 23, 24)
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
rh_sulc_map_binary = np.where(rh_sulc < 0.0, 0.25, 0.6)


# Assuming `bin` is defined and `df` contains your data
print(bin)
fig = plt.figure(figsize=(11, 8))
groups = ["fetal<29", "fetal>29", "neonates<37",
          "neonates>37", '12-16y', '18-21y']
colors = {'neonates<37': 'limegreen', 'neonates>37': 'mediumaquamarine', 'fetal<29': 'gold',
          'fetal>29': 'yellowgreen', '12-16y': 'cadetblue', '18-21y': 'steelblue'}
text_label = {'neonates<37': 'Preterm neonatal', 'neonates>37': 'Full-term neonatal',
              'fetal<29': 'Fetal second trimester', 'fetal>29': 'Fetal third trimester', '12-16y': 'Adolescent', '18-21y': 'Adult'}
left_label = {'neonates<37': 'c', 'neonates>37': 'd',
              'fetal<29': 'a', 'fetal>29': 'b', '12-16y': 'e', '18-21y': 'f'}
for i, group in enumerate(groups):
    if group == "neonates<37":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_L_label-eccentricity_desc-real_roi-v2th00_metric_less_37.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_R_label-eccentricity_desc-real_roi-v2th00_metric_less_37.gii"
    elif group == "neonates>37":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_L_label-eccentricity_desc-real_roi-v2th00_metric_over_37.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_R_label-eccentricity_desc-real_roi-v2th00_metric_over_37.gii"
    elif group == "fetal<29":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_younger_fetal_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_younger_fetal_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "fetal>29":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_older_fetal_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_older_fetal_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "12-16y":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_L_label-eccentricity_desc-real_roi-v2th00_metric_young.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_R_label-eccentricity_desc-real_roi-v2th00_metric_young.gii"
    elif group == "18-21y":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_L_label-eccentricity_desc-real_roi-v2th00_metric_old.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_R_label-eccentricity_desc-real_roi-v2th00_metric_old.gii"

    color_range = (0, 20)

    fig_L = plotting.plot_surf_stat_map(
        lh, average_left, bg_map=lh_sulc_map_binary,
        hemi='left', view=(-10, -70.0), darkness=1.0,
        colorbar=False, cmap="jet", vmax=20,
        threshold=0.1)
    fig_R = plotting.plot_surf_stat_map(
        rh, average_right, bg_map=rh_sulc_map_binary,
        hemi='right', view=(-10, -110.0), darkness=1.0,
        colorbar=False, cmap="jet", vmax=20,
        threshold=0.1)
    # Adjust subplot parameters to minimize white space
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    fig_L.savefig(
        f'/data/p_02915/dhcp_derivatives_SPOT/Figures_v2/{group}_L_eccentricity.pdf')
    fig_R.savefig(
        f'/data/p_02915/dhcp_derivatives_SPOT/Figures_v2/{group}_R_eccentricity.pdf')
