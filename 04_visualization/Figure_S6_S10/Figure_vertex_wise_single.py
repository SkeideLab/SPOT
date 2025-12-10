import matplotlib.pyplot as plt
from nilearn import surface
from neuromaps.datasets import fetch_fsaverage
from nilearn import plotting
import numpy as np
import nibabel as nib
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
    indices_area = np.nonzero(
        np.logical_or(
            visparc.agg_data() == labels_area[0], visparc.agg_data() == labels_area[1]
        )
    )[0]
    return indices_area

LABELS_V2 = (2, 3)

# Parameters for plot configuration
parameters = ["eccentricity","polarangle"]
num_cols = 6  # 2 columns for L and R hemispheres in both plots
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

text_label = {'preterm': 'preterm', 'fullterm': 'full-term',
              '2nd': 'tri-2nd', '3rd': 'tri-3rd', 'adolescent': 'adolescent', 'adult': 'adult'}
# Load the GIFTI file (atlas file)
atlas_lh = nib.load(f'{file_path}/hemi-L_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii')
atlas_rh = nib.load(f'{file_path}/hemi-R_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii')

indices_left = get_indices_roi(LABELS_V2, atlas_lh)
indices_right = get_indices_roi(LABELS_V2, atlas_rh)
# Extract the label data
atlas_data_l = atlas_lh.darrays[0].data
atlas_data_r = atlas_rh.darrays[0].data  # Extract the data array

# Create two figures, one for each parameter
for param in parameters:
    fig, axes = plt.subplots(5, ncols=num_cols, figsize=(8, 6), subplot_kw={'projection': '3d'})  # Adjust the figsize as needed

    # Flatten the axes for easy indexing
    group = ["2nd", "3rd", "preterm", "fullterm", "adolescent", "adult"]
    from itertools import combinations
    group_pairs = list(combinations(group, 2))

    for idx, (groupA_name, groupB_name) in enumerate(group_pairs):        
        # Load the corresponding GIFTI files for left and right hemispheres
        file_L = f'{file_path}/{groupA_name}_{groupB_name}_L_{param}_vertexwise.shape.gii'
        file_R = f'{file_path}/{groupA_name}_{groupB_name}_R_{param}_vertexwise.shape.gii'
        
        # Load data from GIFTI files
        data_L = surface.load_surf_data(file_L)
        data_R = surface.load_surf_data(file_R)
        print(f"{groupA_name} vs. {groupB_name}:")
        # Masking data: setting all values not equal to 1 to NaN
        data_L_masked = np.where(data_L == 1, data_L, 0)
        data_R_masked = np.where(data_R == 1, data_R, 0)
        print(f"Left {data_L_masked[data_L_masked == 1].shape[0]/indices_left.shape[0]*100:.2f}% Right {data_R_masked[data_R_masked == 1].shape[0]/indices_left.shape[0]*100:.2f}%")


        if idx < 5:
            ax_left = axes[idx, 0]
            ax_right = axes[idx, 1]
            fig.text(0.18, 1 - (idx * 0.195)-0.012,
             f'{text_label[groupA_name]} vs {text_label[groupB_name]}', ha='center', va='center', fontsize=12)
        elif idx >= 10:
            ax_left = axes[idx-10, 4]
            ax_right = axes[idx-10, 5]
            fig.text(0.82, 1 - ((idx-10) * 0.195)-0.012,
             f'{text_label[groupA_name]} vs {text_label[groupB_name]}', ha='center', va='center', fontsize=12)
        else:
            ax_left = axes[idx-5, 2]
            ax_right = axes[idx-5, 3]
            fig.text(0.5, 1 - ((idx-5) * 0.195)-0.012,
             f'{text_label[groupA_name]} vs {text_label[groupB_name]}', ha='center', va='center', fontsize=12)
            
        # Plot left hemisphere in the first column
        plotting.plot_surf_stat_map(
        lh, data_L_masked, bg_map=lh_sulc_map_binary,
        hemi='left', view=(0, 290), axes=ax_left, darkness=1.0,
        threshold=0.1, vmin=0, vmax=10, colorbar=False, cmap="autumn", avg_method="max")
        plotting.plot_surf_contours(
        surf_mesh=lh, roi_map=atlas_data_l, hemi='left',
        view=(0, 290), axes=ax_left, levels=[2, 1, 3], colors=["k", 'k', 'k'],
        )

        plotting.plot_surf_stat_map(
        rh, data_R_masked, bg_map=rh_sulc_map_binary,
        hemi='right', view=(0, 250), axes=ax_right, darkness=1.0,
        threshold=0.1, vmin=0, vmax=10, colorbar=False, cmap="autumn", avg_method="max")
        plotting.plot_surf_contours(
        surf_mesh=rh, roi_map=atlas_data_r, hemi='right',
        view=(0, 250), axes=ax_right, levels=[2, 1, 3], colors=["k", 'k', 'k'],
        )

        ax_left.set_xlim([-50, 50])
        ax_left.set_ylim([-50, 50])
        ax_left.set_zlim([-50, 50])   
        ax_right.set_xlim([-50, 50])
        ax_right.set_ylim([-50, 50])
        ax_right.set_zlim([-50, 50])   

    fig.text(0.05, 1 - 0.05,
             "L", ha='center', va='center', fontsize=12)
    fig.text(0.10, 1 - 0.09,
             "-V3", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))    
    fig.text(0.105, 1 - 0.105,
             "-V2", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))    
    fig.text(0.30, 1-0.05,
             "R", ha='center', va='center', fontsize=12)
    fig.text(0.255, 1 - 0.09,
             "V3-", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))    
    fig.text(0.25, 1 - 0.103,
             "V2-", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))
    
    fig.text(0.37, 1 - 0.05,
             "L", ha='center', va='center', fontsize=12)    
    fig.text(0.425, 1 - 0.09,
             "-V3", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))    
    fig.text(0.43, 1 - 0.105,
             "-V2", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))
    
    fig.text(0.62, 1-0.05,
             "R", ha='center', va='center', fontsize=12)    
    fig.text(0.575, 1 - 0.09,
             "V3-", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))    
    fig.text(0.57, 1 - 0.103,
             "V2-", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))
    
    fig.text(0.70, 1 - 0.05,
             "L", ha='center', va='center', fontsize=12)    
    fig.text(0.75, 1 - 0.09,
             "-V3", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))    
    fig.text(0.755, 1 - 0.105,
             "-V2", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))
    
    fig.text(0.95, 1-0.05,
             "R", ha='center', va='center', fontsize=12)
    fig.text(0.9, 1 - 0.09,
             "V3-", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))    
    fig.text(0.895, 1 - 0.103, 
             "V2-", ha='center', va='center', fontsize=6, backgroundcolor = "w", bbox=dict(facecolor='w', edgecolor='none', pad=0.1))

 
    if param =="eccentricity":
        fig_name = "S6"
    elif param == "polarangle":
        fig_name = "S10"
    # Adjust layout and save the plot
    plt.tight_layout()
    plt.savefig(f'{file_path}/Figure{fig_name}.png', dpi=600)
    print(param)
    #plt.show()