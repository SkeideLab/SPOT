import matplotlib.pyplot as plt
from nilearn import plotting, surface, datasets
import numpy as np
from matplotlib import cm
import os

# Get the full path of the script
file_path = os.path.dirname(os.path.abspath(__file__))

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation 
    about the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    return np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                     [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                     [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])

fsaverage = datasets.fetch_surf_fsaverage(mesh='fsaverage')
fs_sulc_l = surface.load_surf_data(fsaverage.sulc_left)
fs_sulc_r = surface.load_surf_data(fsaverage.sulc_right)
# Step 1: Define file paths for native surfaces, sulcal maps, and atlases for each individual
individuals = [
    {'surf_lh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-L_space-fslr_sphere.rot.surf.gii', 
     'sulc_lh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-left_sulc.shape.gii', 
     'atlas_lh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-L_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii',
     'surf_rh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-R_space-fslr_sphere.rot.surf.gii', 
     'sulc_rh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-right_sulc.shape.gii',
     'atlas_rh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-R_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii'},
    {'surf_lh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-L_space-fslr_sphere.rot.surf.gii', 
     'sulc_lh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-left_sulc.shape.gii', 
     'atlas_lh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-L_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii',
     'surf_rh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-R_space-fslr_sphere.rot.surf.gii', 
     'sulc_rh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-right_sulc.shape.gii', 
     'atlas_rh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-R_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii'},
    {'surf_lh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-left_sphere.surf.gii', 
     'sulc_lh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-left_sulc.shape.gii', 
     'atlas_lh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-L_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii',
     'surf_rh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-right_sphere.surf.gii', 
     'sulc_rh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-right_sulc.shape.gii', 
     'atlas_rh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-R_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii'},
    {'surf_lh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-left_sphere.surf.gii', 
     'sulc_lh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-left_sulc.shape.gii', 
     'atlas_lh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-L_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii',
     'surf_rh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-right_sphere.surf.gii', 
     'sulc_rh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-right_sulc.shape.gii', 
     'atlas_rh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-R_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii'},
    {'surf_lh': f'{file_path}/HCD0478356_V1_MR.L.sphere.rot.native.surf.gii', 
     'sulc_lh': f'{file_path}/HCD0478356_V1_MR.L.sulc.native.shape.gii', 
     'atlas_lh': f'{file_path}/HCD0478356_V1_MR_hemi-L_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii',
     'surf_rh': f'{file_path}/HCD0478356_V1_MR.R.sphere.rot.native.surf.gii', 
     'sulc_rh': f'{file_path}/HCD0478356_V1_MR.R.sulc.native.shape.gii', 
     'atlas_rh': f'{file_path}/HCD0478356_V1_MR_hemi-R_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii'},
     {'surf_lh': f'{file_path}/HCD2336346_V1_MR.L.sphere.rot.native.surf.gii', 
     'sulc_lh': f'{file_path}/HCD2336346_V1_MR.L.sulc.native.shape.gii', 
     'atlas_lh': f'{file_path}/HCD2336346_V1_MR_hemi-L_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii',
     'surf_rh': f'{file_path}/HCD2336346_V1_MR.R.sphere.rot.native.surf.gii', 
     'sulc_rh': f'{file_path}/HCD2336346_V1_MR.R.sulc.native.shape.gii', 
     'atlas_rh': f'{file_path}/HCD2336346_V1_MR_hemi-R_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii'},
    {'surf_lh': fsaverage.sphere_left, 
     'sulc_lh': fs_sulc_l, 
     'atlas_lh': f'{file_path}/hemi-L_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii',
     'surf_rh': fsaverage.sphere_right, 
     'sulc_rh': fs_sulc_r, 
     'atlas_rh': f'{file_path}/hemi-R_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii'}
]

# Step 2: Labels you want to highlight in the atlas (e.g., 1, 2, and 3)
labels_to_display = [1, 2, 3]
# Step 3: Create a 2x6 grid for plotting
fig, axes = plt.subplots(7, 2, subplot_kw={'projection': '3d'}, figsize=(4, 12))
left_view = [(0, 320),
             (0, 310),
             (0, 300),
             (0, 300),
             (10, 330),
             (10, 330),
             (0, 300),
]
right_view = [(10, 240),
             (30, 240),
             (0, 240),
             (0, 250),
             (10, 220),
             (10, 220),
             (0, 240),
]
lim = [
    (-85, 85),
    (-85, 85),
    (-0.85, 0.85),
    (-0.85, 0.85),
    (-85, 85),
    (-85, 85),
    (-85, 85),
]
group = [
    "tri-2 prenatal (26w)",
    "tri-3 prenatal (32w)",
    "preterm neonatal (33w/35w)",
    "full-term neonatal (40w/40w)",
    "adolescent (13.9y)",
    "adult (20.3y)",
    "template space"
]
# Step 4: Loop through all individuals and plot their native surfaces with selected atlas labels
for i, ind in enumerate(individuals):
    # Plot left hemisphere in the first row
    ax_lh = axes[i, 0]  # Get the subplot axis for left hemisphere
    surf_lh = surface.load_surf_mesh(ind['surf_lh'])
    sulc_lh = surface.load_surf_data(ind['sulc_lh'])
    atlas_lh = surface.load_surf_data(ind['atlas_lh'])
    if i == 6:
        lh_sulc_map_binary = np.where(sulc_lh < 0.0, 0.25, 0.6)
    else:    
        lh_sulc_map_binary = np.where(sulc_lh > 0.0, 0.25, 0.6)

    # Create a mask for left hemisphere labels
    mask_lh = np.isin(atlas_lh, labels_to_display)
    texture_lh = np.zeros_like(atlas_lh)
    texture_lh[mask_lh] = atlas_lh[mask_lh]
      
    # Plot left hemisphere surface
    plotting.plot_surf_roi(surf_lh, roi_map=texture_lh, hemi='left', view=left_view[i], bg_map=lh_sulc_map_binary, 
                           axes=ax_lh, figure=fig, cmap='viridis', colorbar=False, threshold=0.001)
    ax_lh.set_xlim(lim[i])
    ax_lh.set_ylim(lim[i])
    ax_lh.set_zlim(lim[i])
    #ax_lh.set_aspect("equal")
    #print(f'{group[i]}_{ax_lh.get_xlim()}')
    #print(ax_lh.get_ylim())
    #print(ax_lh.get_zlim())

    # Plot right hemisphere in the second row
    ax_rh = axes[i, 1]  # Get the subplot axis for right hemisphere
    surf_rh = surface.load_surf_mesh(ind['surf_rh'])
    sulc_rh = surface.load_surf_data(ind['sulc_rh'])
    atlas_rh = surface.load_surf_data(ind['atlas_rh'])
    if i == 6:
        rh_sulc_map_binary = np.where(sulc_rh < 0.0, 0.25, 0.6)
    else:    
        rh_sulc_map_binary = np.where(sulc_rh > 0.0, 0.25, 0.6)

    # Create a mask for right hemisphere labels
    mask_rh = np.isin(atlas_rh, labels_to_display)
    texture_rh = np.zeros_like(atlas_rh)
    texture_rh[mask_rh] = atlas_rh[mask_rh]

    # Plot right hemisphere surface
    plotting.plot_surf_roi(surf_rh, roi_map=texture_rh, hemi='right', view=right_view[i], bg_map=rh_sulc_map_binary,
                           axes=ax_rh, figure=fig, cmap='viridis', colorbar=False, threshold=0.001)
    ax_rh.set_xlim(lim[i])
    ax_rh.set_ylim(lim[i])
    ax_rh.set_zlim(lim[i])
    #ax_rh.set_aspect("equal")

    #print(ax_rh.get_xlim())
    #print(ax_rh.get_ylim())
    #print(ax_rh.get_zlim())

plt.tight_layout()   

for i in range(0,7):
    if i == 0:
        fig.text(0.5, 1 - 0.01,
            f'{group[i]}', ha='center', va='center', fontsize=16)
    else:
        fig.text(0.5, 1 - (i * 0.14) - 0.015,
            f'{group[i]}', ha='center', va='center', fontsize=16)
    
# Add "Left" and "Right" text labels
fig.text(0.08, 0.965, 'L', fontsize=14, va='center',)
fig.text(0.90, 0.965, 'R', fontsize=14, va='center',)

#plt.tight_layout()
# Step 5: Adjust layout and display
#plt.show()
fig.savefig(f"{file_path}/FigureS2.png",
                dpi=300, bbox_inches='tight')