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


df = pd.DataFrame()
for group in ["preterm", "fullterm", "2nd", "3rd", "adolescent", "adult", "benson"]:
    if group == "preterm":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_v2.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fullterm":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_v2.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "2nd":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "3rd":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "adolescent":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "adult":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
            "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        )
        subject_info = pd.read_csv(
            '/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "benson":
        PREFIX_MODEL=("/data/p_02915/SPOT/01_dataprep/retinotopy/"
                      "templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-eccentretinotbenson2014_seg.shape.gii")

    for hemi in ["L", "R"]:
        visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
        indices_v2 = get_indices_roi(LABELS_V2, visparc)
        group_name = f"{group}_{hemi}"
        print(group_name)
        parameters = []
    # FORMAT PATHS FOR INPUT AND OUTPUT
        if group == "adolescent" or group == "adult":
            for index, row in subject_info.iterrows():            
                sub_id = subject_info.at[index, "sub_id"]
                prefix_model = PREFIX_MODEL.format(
                    sub=sub_id,
                    hemi=hemi,
                )
                input_path = prefix_model

                ccf = surface.load_surf_data(input_path)
                ccf_v0 = ccf[indices_v2].astype(np.float64)
                parameters.append(ccf_v0)            
        elif group == "benson":
            prefix_model = PREFIX_MODEL.format(
                    hemi=hemi,
                )
            input_path = prefix_model
            ccf = surface.load_surf_data(input_path)
            ccf_v0 = ccf[indices_v2].astype(np.float64)
            parameters.append(ccf_v0)
            if hemi == "L":
                benson_L = np.array(ccf_v0)                
            elif hemi == "R":
                benson_R = np.array(ccf_v0)
        else:
            for index, row in subject_info.iterrows():            
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

        if group == "preterm" and hemi == "L":
            neonates_y_L = np.array(parameters)
        elif group == "fullterm" and hemi == "L":
            neonates_o_L = np.array(parameters)
        elif group == "2nd" and hemi == "L":
            fetal_y_L = np.array(parameters)
        elif group == "3rd" and hemi == "L":
            fetal_o_L = np.array(parameters)
        elif group == "adolescent" and hemi == "L":
            hcp_y_L = np.array(parameters)
        elif group == "adult" and hemi == "L":
            hcp_o_L = np.array(parameters)
        elif group == "preterm" and hemi == "R":
            neonates_y_R = np.array(parameters)
        elif group == "fullterm" and hemi == "R":
            neonates_o_R = np.array(parameters)
        elif group == "2nd" and hemi == "R":
            fetal_y_R = np.array(parameters)
        elif group == "3rd" and hemi == "R":
            fetal_o_R = np.array(parameters)
        elif group == "adolescent" and hemi == "R":
            hcp_y_R = np.array(parameters)
        elif group == "adult" and hemi == "R":
            hcp_o_R = np.array(parameters)

        data = np.array(parameters).reshape(-1)
        data_flatten = data.flatten()
        data_r = pd.DataFrame({'ecc': data_flatten, 'Group': group_name})
        df = pd.concat([df, data_r], ignore_index=True)

bin = np.linspace(0, 90, 19)
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
left_temp = nib.load(VISPARC_PATH.format(hemi="L"))
indices_left = get_indices_roi([1, 2, 3], left_temp)
right_temp = nib.load(VISPARC_PATH.format(hemi="R"))
indices_right = get_indices_roi([1, 2, 3], right_temp)

# Assuming `bin` is defined and `df` contains your data
print(bin)
fig = plt.figure(figsize=(11, 8))
groups = ["2nd", "3rd", "preterm",
          "fullterm", 'adolescent', 'adult', "benson"]
text_label = {'preterm': 'preterm neonatal', 'fullterm': 'full-term neonatal',
              '2nd': 'tri-2nd prenatal', '3rd': 'tri-3rd prenatal', 'adolescent': 'adolescent', 'adult': 'adult', "benson": "template"}
left_label = ["a", "b", "c", "d", "e", "f", "g"]
colors = {'preterm': 'limegreen', 'fullterm': 'mediumaquamarine', '2nd': 'gold',
          '3rd': 'yellowgreen', 'adolescent': 'cadetblue', 'adult': 'steelblue', "benson":"#656364"}
vmax = 20
# Load the GIFTI file (atlas file)
atlas_lh = nib.load('/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-L_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii')
atlas_rh = nib.load('/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-R_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii')
template_lh = nib.load(lh)  # Load the left hemisphere surface
template_rh = nib.load(rh)

# Extract the label data
atlas_data_l = atlas_lh.darrays[0].data
atlas_data_r = atlas_rh.darrays[0].data  # Extract the data array
# Specify the labels of interest
labels_of_interest = [1, 2, 3]


for i, group in enumerate(groups):
    data_L = df[df['Group'] == f"{group}_L"]
    data_R = df[df['Group'] == f"{group}_R"]

    # Histogram positions
    hist_L_pos = [0.1, 1 - i*0.13, 0.12, 0.1]
    # Adjusted position for the second histogram
    hist_R_pos = [0.62, 1 - i*0.13, 0.12, 0.1]

    ax_hist_L = fig.add_axes(hist_L_pos)
    ax_hist_R = fig.add_axes(hist_R_pos)

    sns.histplot(data=data_L, x='ecc', ax=ax_hist_L,
                 color=colors[group], bins=bin, stat='density')
    ax_hist_L.set_xlim(0, 90)
    ax_hist_L.set_ylim(0, 0.12)
    ax_hist_L.set_xlabel("Eccentricity" if i == 6 else "", labelpad=1)
    ax_hist_L.set_ylabel("Probability\ndensity" if i == 6 else "",)
    ax_hist_L.set_xticks([0, 45, 90])

    sns.histplot(data=data_R, x='ecc', ax=ax_hist_R,
                 color=colors[group], bins=bin, stat='density')
    ax_hist_R.set_xlim(0, 90)
    ax_hist_R.set_ylim(0, 0.12)
    ax_hist_R.set_xlabel("Eccentricity" if i == 6 else "", labelpad=1)
    ax_hist_R.set_ylabel("Probability\ndensity" if i == 6 else "",)
    ax_hist_R.set_xticks([0, 45, 90])

    
    if group == "preterm":
        average_left = f"/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel/Averaged_younger_n_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        average_right = f"/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel/Averaged_younger_n_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "fullterm":
        average_left = f"/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel/Averaged_older_n_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        average_right = f"/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel/Averaged_older_n_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "2nd":
        average_left = f"/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_younger_fetal_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        average_right = f"/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_younger_fetal_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "3rd":
        average_left = f"/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_older_fetal_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        average_right = f"/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_older_fetal_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "adolescent":
        average_left = f"/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_young_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        average_right = f"/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_young_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "adult":
        average_left = f"/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_old_L_label-eccentricity_desc-real_roi-v2th00_metric.gii"
        average_right = f"/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_old_R_label-eccentricity_desc-real_roi-v2th00_metric.gii"
    elif group == "benson":
        template_left = "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-L_space-fsaverage_dens-164k_desc-eccentretinotbenson2014_seg.shape.gii"
        template_right = "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-R_space-fsaverage_dens-164k_desc-eccentretinotbenson2014_seg.shape.gii"
        ccf_l = surface.load_surf_data(template_left)
        average_left = np.zeros(ccf_l.shape)
        average_left[indices_left]= ccf_l[indices_left].astype(np.float64)
        ccf_r = surface.load_surf_data(template_right)
        average_right = np.zeros(ccf_r.shape)
        average_right[indices_right]= ccf_r[indices_right].astype(np.float64)

    color_range = (0, vmax)

    fig.text(0.395, 1 - (i * 0.13) + 0.11,
             text_label[group], ha='center', va='center', fontsize=14)

    # Create 3D axes for brain surface images and adjust positions to overlap by one-third
    ax1 = fig.add_axes([0.205, 1 - i*0.13-0.01, 0.12, 0.11], projection='3d')
    ax2 = fig.add_axes([0.29, 1 - i*0.13, 0.12, 0.1], projection='3d')
    ax3 = fig.add_axes([0.375, 1 - i*0.13, 0.12, 0.1], projection='3d')
    ax4 = fig.add_axes([0.455, 1 - i*0.13-0.01, 0.12, 0.11], projection='3d')
    
    plotting.plot_surf_stat_map(
        lh, average_left, bg_map=lh_sulc_map_binary,
        hemi='left', view=(0, 290), axes=ax1, darkness=1.0,
        colorbar=False, cmap="jet", vmax=vmax,
        threshold=0.000001)    
    
    plotting.plot_surf_contours(
        surf_mesh=lh, roi_map=atlas_data_l, hemi='left',
        view=(0, 290), axes=ax1, levels=[2, 1, 3], colors=['w','w', 'w'],
    )
    plotting.plot_surf_stat_map(
        lh, average_left, bg_map=lh_sulc_map_binary,
        hemi='left', view=(-10, -70.0), axes=ax2, darkness=1.0,
        colorbar=False, cmap="jet", vmax=vmax,
        threshold=0.000001)
    plotting.plot_surf_contours(
        surf_mesh=lh, roi_map=atlas_data_l, hemi='left',
        view=(-10, -70), axes=ax2, levels=[2, 1, 3], colors=["w", 'w', 'w'],
    )
    plotting.plot_surf_stat_map(
        rh, average_right, bg_map=rh_sulc_map_binary,
        hemi='right', view=(-10, -110.0), axes=ax3, darkness=1.0,
        colorbar=False, cmap="jet", vmax=vmax,
        threshold=0.000001)
    plotting.plot_surf_contours(
        surf_mesh=rh, roi_map=atlas_data_r, hemi='right',
        view=(-10, -110), axes=ax3, levels=[2, 1, 3], colors=["w", 'w', 'w'],
    )
    plotting.plot_surf_stat_map(
        rh, average_right, bg_map=rh_sulc_map_binary,
        hemi='right', view=(0, 250), axes=ax4, darkness=1.0,
        colorbar=False, cmap="jet", vmax=vmax,
        threshold=0.000001)
    plotting.plot_surf_contours(
        surf_mesh=rh, roi_map=atlas_data_r, hemi='right',
        view=(0, 250), axes=ax4, levels=[2, 1, 3], colors=["w", 'w', 'w'],
    )
    
    # Remove axis ticks and labels for brain surface plots
    ax1.set_axis_off()
    ax2.set_axis_off()
    ax3.set_axis_off()
    ax4.set_axis_off()

    ax1.set_xlim([-50, 50])
    ax1.set_ylim([-50, 50])
    ax1.set_zlim([-50, 50])
    ax2.set_xlim([-50, 10])
    ax2.set_ylim([-35, 35])
    ax2.set_zlim([-30, 35])
    ax3.set_xlim([-10, 50])
    ax3.set_ylim([-35, 35])
    ax3.set_zlim([-30, 35])
    ax4.set_xlim([-50, 50])
    ax4.set_ylim([-50, 50])
    ax4.set_zlim([-50, 50])
        

# Adjust layout to leave space for the colorbar at the bottom
plt.subplots_adjust(bottom=0.15, hspace=0.5)
#for i in range(7):
#        fig.text(0.085, 1 - (i * 0.13)+0.108,
#                left_label[i], ha='center', va='center', fontsize=12, weight='bold')
        
fig.text(0.24, 1.09, "L",  ha='center', va='center', fontsize=12)
fig.text(0.54, 1.09, "R",  ha='center', va='center', fontsize=12)
fig.text(0.392, 1.09, "———V3———",  ha='center', va='center', fontsize=5,)
fig.text(0.392, 1.075, "———V2———",  ha='center', va='center', fontsize=5,)
fig.text(0.392, 0.31, "———V3———",  ha='center', va='center', fontsize=5,)
fig.text(0.392, 0.295, "———V2———",  ha='center', va='center', fontsize=5,)
fig.text(0.392, 0.27, "———V1———",  ha='center', va='center', fontsize=5, )
# Add colorbar
# Adjust the position and height of the colorbar axis
cbar_ax = fig.add_axes([0.24, 0.18, 0.3, 0.02])
cbar = plt.colorbar(cm.ScalarMappable(cmap='jet', norm=plt.Normalize(
    vmin=0, vmax=vmax)), cax=cbar_ax, orientation='horizontal', fraction=.1)
cbar.set_label('Eccentricity')
# Set the ticks of the colorbar to include the min and max values
cbar.set_ticks([0, vmax/2, vmax])

# Set the tick labels to show min and max values (you can format them as needed)
cbar.set_ticklabels(['0', f'{vmax/2:.0f}', f'{vmax:.0f}'])

# Save the figure
fig.savefig("/data/p_02915/dhcp_derivatives_SPOT/Figures_v3/ECC_all_2.png",
            dpi=300, bbox_inches='tight')
