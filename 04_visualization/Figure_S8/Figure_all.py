import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from neuromaps.datasets import fetch_fsaverage
import matplotlib.pyplot as plt
from nilearn import plotting
import seaborn as sns
from matplotlib import cm
import os

# Get the full path of the script
file_path = os.path.dirname(os.path.abspath(__file__))


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
    "{file_path}/hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii")
LABELS_V2 = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

surfaces = fetch_fsaverage(density="164k")
lh, rh = surfaces["inflated"]
# Load sulcal depth data
# Load the left hemisphere sulcal depth data
lh_sulc = surface.load_surf_data(
    f"{file_path}/lh.curv")
# Load the right hemisphere sulcal depth data
rh_sulc = surface.load_surf_data(
    f"{file_path}/rh.curv")
lh_sulc_map_binary = np.where(lh_sulc < 0.0, 0.25, 0.6)
rh_sulc_map_binary = np.where(rh_sulc < 0.0, 0.25, 0.6)


figs, axs = plt.subplots(7, 4, figsize=(4, 8), subplot_kw={'projection': '3d'})
groups = ["2nd", "3rd", "preterm",
          "fullterm", 'adolescent', 'adult', "benson"]
text_label = {'preterm': 'preterm neonatal', 'fullterm': 'full-term neonatal',
              '2nd': 'tri-2nd prenatal', '3rd': 'tri-3rd prenatal', 'adolescent': 'adolescent', 'adult': 'adult', "benson": "template"}
left_label = {'preterm': 'c', 'fullterm': 'd',
              '2nd': 'a', '3rd': 'b', 'adolescent': 'e', 'adult': 'f', "benson": "g"}
left_temp = nib.load(VISPARC_PATH.format(file_path=file_path, hemi="L"))
indices_left = get_indices_roi(LABELS_V2, left_temp)
indices_left_v1 = get_indices_roi([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], left_temp)
right_temp = nib.load(VISPARC_PATH.format(file_path=file_path, hemi="R"))
indices_right = get_indices_roi(LABELS_V2, right_temp)
indices_right_v1 = get_indices_roi([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], right_temp)

for index, param in enumerate(["eccentricity", "polarangle"]):
    
    for i, group in enumerate(groups):
        if group == "preterm":
            average_left = f"{file_path}/Averaged_less_n_L_label-{param}_desc-real_roi-v2th00_metric.gii"
            average_right = f"{file_path}/Averaged_less_n_R_label-{param}_desc-real_roi-v2th00_metric.gii"
        elif group == "fullterm":
            average_left = f"{file_path}/Averaged_older_n_L_label-{param}_desc-real_roi-v2th00_metric.gii"
            average_right = f"{file_path}/Averaged_older_n_R_label-{param}_desc-real_roi-v2th00_metric.gii"
        elif group == "2nd":
            average_left = f"{file_path}/Averaged_younger_fetal_L_label-{param}_desc-real_roi-v2th00_metric.gii"
            average_right = f"{file_path}/Averaged_younger_fetal_R_label-{param}_desc-real_roi-v2th00_metric.gii"
        elif group == "3rd":
            average_left = f"{file_path}/Averaged_older_fetal_L_label-{param}_desc-real_roi-v2th00_metric.gii"
            average_right = f"{file_path}/Averaged_older_fetal_R_label-{param}_desc-real_roi-v2th00_metric.gii"
        elif group == "adolescent":
            average_left = f"{file_path}/Averaged_young_L_label-{param}_desc-real_roi-v2th00_metric.gii"
            average_right = f"{file_path}/Averaged_young_R_label-{param}_desc-real_roi-v2th00_metric.gii"
        elif group == "adult":
            average_left = f"{file_path}/Averaged_old_L_label-{param}_desc-real_roi-v2th00_metric.gii"
            average_right = f"{file_path}/Averaged_old_R_label-{param}_desc-real_roi-v2th00_metric.gii"
        elif group == "benson":           
            if param == "eccentricity":
                template_left = f"{file_path}/hemi-L_space-fsaverage_dens-164k_desc-eccentretinotbenson2014_seg.shape.gii"
                template_right = f"{file_path}/hemi-R_space-fsaverage_dens-164k_desc-eccentretinotbenson2014_seg.shape.gii"
            elif param == "polarangle":
                template_left = f"{file_path}/hemi-L_space-fsaverage_dens-164k_desc-angleretinotbenson2014_seg.shape.gii"
                template_right = f"{file_path}/hemi-R_space-fsaverage_dens-164k_desc-angleretinotbenson2014_seg.shape.gii"
            ccf_l = surface.load_surf_data(template_left)
            average_left = np.zeros(ccf_l.shape)
            average_left[indices_left_v1]= ccf_l[indices_left_v1].astype(np.float64)
            ccf_r = surface.load_surf_data(template_right)
            average_right = np.zeros(ccf_r.shape)
            average_right[indices_right_v1]= ccf_r[indices_right_v1].astype(np.float64)
        
        left_param = surface.load_surf_data(average_left)
        right_param = surface.load_surf_data(average_right)

        if param == "eccentricity":
            vmax = np.log10(20)
            ax1 = axs[i,0]
            ax2 = axs[i,1]
            average_left = np.log10(np.clip(left_param , 1e-6, None))  # Replace values <1e-6 to avoid log errors
            average_right = np.log10(np.clip(right_param, 1e-6, None))
            empty_left = np.full(left_param.shape, np.nan)
            empty_right = np.full(right_param.shape, np.nan)

            if i == 6:
                empty_left[indices_left_v1]= average_left[indices_left_v1].astype(np.float64)
                empty_right[indices_right_v1]= average_right[indices_right_v1].astype(np.float64)
            else:
                empty_left[indices_left]= average_left[indices_left].astype(np.float64)
                empty_right[indices_right]= average_right[indices_right].astype(np.float64)

            average_left_log = empty_left
            average_right_log = empty_right
            plotting.plot_surf_stat_map(
                lh, average_left_log, bg_map=lh_sulc_map_binary,
                hemi='left', view=(0, 290), axes=ax1, darkness=1.0,
                colorbar=False, cmap="jet", vmin=0, vmax=vmax, 
                )        
            
            plotting.plot_surf_stat_map(
                rh, average_right_log, bg_map=rh_sulc_map_binary,
                hemi='right', view=(0, 250), axes=ax2, darkness=1.0,
                colorbar=False, cmap="jet", vmin=0, vmax=vmax, 
                )
            


        elif param == "polarangle":
            figs.text(0.5, 0.88 - (i * 0.113),
             text_label[group], ha='center', va='center', fontsize=14)
            vmax = 180
            ax1 = axs[i,2]
            ax2 = axs[i,3]

            plotting.plot_surf_stat_map(
                lh, average_left, bg_map=lh_sulc_map_binary,
                hemi='left', view=(0, 290), axes=ax1, darkness=1.0,
                colorbar=False, cmap="jet", vmin=0, vmax=vmax, 
                threshold=1e-10)        
            
            plotting.plot_surf_stat_map(
                rh, average_right, bg_map=rh_sulc_map_binary,
                hemi='right', view=(0, 250), axes=ax2, darkness=1.0,
                colorbar=False, cmap="jet", vmin=0, vmax=vmax, 
                threshold=1e-10)
        
        ax1.set_xlim([-60, 10])
        ax1.set_ylim([-40, 40])
        ax1.set_zlim([-40, 40])
        ax2.set_xlim([-10, 60])
        ax2.set_ylim([-40, 40])
        ax2.set_zlim([-40, 40])
        
        

# Adjust layout to leave space for the colorbar at the bottom
plt.subplots_adjust(bottom=0.1, hspace=0.1)

# Add colorbar
# Adjust the position and height of the colorbar axis
cbar_ax = figs.add_axes([0.15, 0.08, 0.32, 0.02])
cbar = plt.colorbar(cm.ScalarMappable(cmap='jet', norm=plt.Normalize(
    vmin=0, vmax=np.log10(20))), cax=cbar_ax, orientation='horizontal', fraction=.1)
# Set the ticks of the colorbar to include the min and max values
cbar.set_ticks([0, np.log10(20)/2, np.log10(20)])

# Set the tick labels to show min and max values (you can format them as needed)
cbar.set_ticklabels([ '0', 'eccentricity',f'log20'])
# Add colorbar
# Adjust the position and height of the colorbar axis
cbar_ax = figs.add_axes([0.55, 0.08, 0.32, 0.02])
cbar = plt.colorbar(cm.ScalarMappable(cmap='jet', norm=plt.Normalize(
    vmin=0, vmax=180)), cax=cbar_ax, orientation='horizontal', fraction=.1)
cbar.set_ticks([0, 90, 180])
cbar.set_ticklabels([ 0, "polar angle",f'{vmax:.0f}'])

#plt.show()
figs.savefig(f"{file_path}/FigureS8.png",
                dpi=300, bbox_inches='tight')
