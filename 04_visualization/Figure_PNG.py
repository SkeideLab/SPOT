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
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "neonates>37":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal<29":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal>29":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "12-16y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "18-21y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric.gii"
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

bin = np.linspace(0, 180, 37)
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
    data_L = df[df['Group'] == f"{group}_L"]
    data_R = df[df['Group'] == f"{group}_R"]

    # Histogram positions
    hist_L_pos = [0.1, 1 - i*0.15, 0.12, 0.1]
    # Adjusted position for the second histogram
    hist_R_pos = [0.61, 1 - i*0.15, 0.12, 0.1]

    ax_hist_L = fig.add_axes(hist_L_pos)
    ax_hist_R = fig.add_axes(hist_R_pos)

    sns.histplot(data=data_L, x='ecc', ax=ax_hist_L,
                 color=colors[group], bins=bin, stat='density')
    ax_hist_L.set_xlim(0, 180)
    ax_hist_L.set_ylim(0, 0.055)
    ax_hist_L.set_xlabel("Polar angle" if i == 5 else "", labelpad=1)
    ax_hist_L.set_ylabel("Probability density" if i == 5 else "",)

    sns.histplot(data=data_R, x='ecc', ax=ax_hist_R,
                 color=colors[group], bins=bin, stat='density')
    ax_hist_R.set_xlim(0, 180)
    ax_hist_R.set_ylim(0, 0.055)
    ax_hist_R.set_xlabel("Polar angle" if i == 5 else "", labelpad=1)
    ax_hist_R.set_ylabel("Probability density" if i == 5 else "",)

    fig.text(0.08, 1 - (i * 0.15) + 0.11,
             left_label[group], ha='center', va='center', fontsize=12, weight='bold')

    if group == "neonates<37":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_L_label-polarangle_desc-real_roi-v2th00_metric_less_37.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_R_label-polarangle_desc-real_roi-v2th00_metric_less_37.gii"
    elif group == "neonates>37":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_L_label-polarangle_desc-real_roi-v2th00_metric_over_37.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_R_label-polarangle_desc-real_roi-v2th00_metric_over_37.gii"
    elif group == "fetal<29":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_younger_fetal_L_label-polarangle_desc-real_roi-v2th00_metric.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_younger_fetal_R_label-polarangle_desc-real_roi-v2th00_metric.gii"
    elif group == "fetal>29":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_older_fetal_L_label-polarangle_desc-real_roi-v2th00_metric.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_older_fetal_R_label-polarangle_desc-real_roi-v2th00_metric.gii"
    elif group == "12-16y":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_L_label-polarangle_desc-real_roi-v2th00_metric_young.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_R_label-polarangle_desc-real_roi-v2th00_metric_young.gii"
    elif group == "18-21y":
        average_left = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_L_label-polarangle_desc-real_roi-v2th00_metric_old.gii"
        average_right = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_R_label-polarangle_desc-real_roi-v2th00_metric_old.gii"

    color_range = (0, 20)

    fig.text(0.395, 1 - (i * 0.15) + 0.115,
             text_label[group], ha='center', va='center', fontsize=14)

    # Create 3D axes for brain surface images and adjust positions to overlap by one-third
    ax1 = fig.add_axes([0.205, 1 - i*0.15-0.01, 0.12, 0.11], projection='3d')
    ax2 = fig.add_axes([0.29, 1 - i*0.15, 0.12, 0.1], projection='3d')
    ax3 = fig.add_axes([0.375, 1 - i*0.15, 0.12, 0.1], projection='3d')
    ax4 = fig.add_axes([0.455, 1 - i*0.15-0.01, 0.12, 0.11], projection='3d')

    plotting.plot_surf_stat_map(
        lh, average_left, bg_map=lh_sulc_map_binary,
        hemi='left', view=(-10, -70.0), axes=ax1, darkness=1.0,
        colorbar=False, cmap="jet", vmax=180,
        threshold=0.1)
    plotting.plot_surf_stat_map(
        lh, average_left, bg_map=lh_sulc_map_binary,
        hemi='left', view=(-10, -70.0), axes=ax2, darkness=1.0,
        colorbar=False, cmap="jet", vmax=180,
        threshold=0.1)
    plotting.plot_surf_stat_map(
        rh, average_right, bg_map=rh_sulc_map_binary,
        hemi='right', view=(-10, -110.0), axes=ax3, darkness=1.0,
        colorbar=False, cmap="jet", vmax=180,
        threshold=0.1)
    plotting.plot_surf_stat_map(
        rh, average_right, bg_map=rh_sulc_map_binary,
        hemi='right', view=(-10, -110.0), axes=ax4, darkness=1.0,
        colorbar=False, cmap="jet", vmax=180,
        threshold=0.1)

    # Remove axis ticks and labels for brain surface plots
    ax1.set_axis_off()
    ax2.set_axis_off()
    ax3.set_axis_off()
    ax4.set_axis_off()

    ax1.set_xlim([-50, 50])
    ax1.set_ylim([-50, 50])
    ax1.set_zlim([-50, 50])
    ax2.set_xlim([-40, 20])
    ax2.set_ylim([-30, 30])
    ax2.set_zlim([-30, 30])
    ax3.set_xlim([-20, 40])
    ax3.set_ylim([-30, 30])
    ax3.set_zlim([-30, 30])
    ax4.set_xlim([-50, 50])
    ax4.set_ylim([-50, 50])
    ax4.set_zlim([-50, 50])

# fig.tight_layout()

# Adjust layout to leave space for the colorbar at the bottom
plt.subplots_adjust(bottom=0.15, hspace=0.5)

# Add colorbar
# Adjust the position and height of the colorbar axis
cbar_ax = fig.add_axes([0.24, 0.18, 0.3, 0.02])
cbar = plt.colorbar(cm.ScalarMappable(cmap='jet', norm=plt.Normalize(
    vmin=0, vmax=180)), cax=cbar_ax, orientation='horizontal', fraction=.1)
cbar.set_label('Polar angle')

# Save the figure
fig.savefig("/data/p_02915/dhcp_derivatives_SPOT/Figures_v2/PNG_all.png",
            dpi=300, bbox_inches='tight')
