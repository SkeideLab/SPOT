import matplotlib.pyplot as plt
from nilearn import plotting, surface, datasets
import numpy as np
from matplotlib import cm
import os

# Get the full path of the script
file_path = os.path.dirname(os.path.abspath(__file__))


fsaverage = datasets.fetch_surf_fsaverage(mesh='fsaverage')
fs_sulc_l = surface.load_surf_data(fsaverage.sulc_left)
fs_sulc_r = surface.load_surf_data(fsaverage.sulc_right)
# Step 1: Define file paths for native surfaces, sulcal maps, and atlases for each individual sub-CC01147XX18_ses-101330_hemi-L_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii
individuals = [
    {'surf_lh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-L_space-fslr_sphere.rot.surf.gii', 
     'sulc_lh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-left_sulc.shape.gii', 
     'atlas_lh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-L_mesh-native_dens-native_desc-real_r.gii',
     'ecc_lh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-L_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_lh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-L_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',
     'surf_rh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-R_space-fslr_sphere.rot.surf.gii', 
     'sulc_rh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-right_sulc.shape.gii',
     'atlas_rh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-R_mesh-native_dens-native_desc-simulated_r.gii',
     'ecc_rh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-R_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_rh': f'{file_path}/sub-CC01008XX10_ses-47430_hemi-R_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',},
    {'surf_lh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-L_space-fslr_sphere.rot.surf.gii', 
     'sulc_lh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-left_sulc.shape.gii', 
     'atlas_lh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-L_mesh-native_dens-native_desc-real_r.gii',
     'ecc_lh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-L_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_lh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-L_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',
     'surf_rh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-R_space-fslr_sphere.rot.surf.gii', 
     'sulc_rh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-right_sulc.shape.gii', 
     'atlas_rh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-R_mesh-native_dens-native_desc-real_r.gii',
     'ecc_rh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-R_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_rh': f'{file_path}/sub-CC01186XX16_ses-126230_hemi-R_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',},
    {'surf_lh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-left_sphere.surf.gii', 
     'sulc_lh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-left_sulc.shape.gii', 
     'atlas_lh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-L_mesh-native_dens-native_desc-real_r.gii',
     'ecc_lh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-L_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_lh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-L_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',
     'surf_rh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-right_sphere.surf.gii', 
     'sulc_rh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-right_sulc.shape.gii', 
     'atlas_rh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-R_mesh-native_dens-native_desc-real_r.gii',
     'ecc_rh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-R_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_rh': f'{file_path}/sub-CC00177XX13_ses-58500_hemi-R_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',},
    {'surf_lh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-left_sphere.surf.gii', 
     'sulc_lh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-left_sulc.shape.gii', 
     'atlas_lh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-L_mesh-native_dens-native_desc-real_r.gii',
     'ecc_lh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-L_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_lh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-L_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',
     'surf_rh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-right_sphere.surf.gii', 
     'sulc_rh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-right_sulc.shape.gii', 
     'atlas_rh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-R_mesh-native_dens-native_desc-real_r.gii',
     'ecc_rh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-R_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_rh': f'{file_path}/sub-CC00320XX07_ses-102300_hemi-R_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',},
    {'surf_lh': f'{file_path}/fs_LR.32k.L.sphere.surf.gii', 
     'sulc_lh': f'{file_path}/fs_LR.32k.L.sulc.shape.gii', 
     'atlas_lh': f'{file_path}/HCD0478356_V1_MR_hemi-L_mesh-native_dens-native_desc-real_r.gii',
     'ecc_lh': f'{file_path}/HCD0478356_V1_MR_hemi-L_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_lh': f'{file_path}/HCD0478356_V1_MR_hemi-L_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',
     'surf_rh': f'{file_path}/fs_LR.32k.R.sphere.surf.gii', 
     'sulc_rh': f'{file_path}/fs_LR.32k.R.sulc.shape.gii', 
     'atlas_rh': f'{file_path}/HCD0478356_V1_MR_hemi-R_mesh-native_dens-native_desc-real_r.gii',
     'ecc_rh': f'{file_path}/HCD0478356_V1_MR_hemi-R_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_rh': f'{file_path}/HCD0478356_V1_MR_hemi-R_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',},
    {'surf_lh': f'{file_path}/fs_LR.32k.L.sphere.surf.gii', 
     'sulc_lh': f'{file_path}/fs_LR.32k.L.sulc.shape.gii', 
     'atlas_lh': f'{file_path}/HCD2336346_V1_MR_hemi-L_mesh-native_dens-native_desc-real_r.gii',
     'ecc_lh': f'{file_path}/HCD2336346_V1_MR_hemi-L_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_lh': f'{file_path}/HCD2336346_V1_MR_hemi-L_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',
     'surf_rh': f'{file_path}/fs_LR.32k.R.sphere.surf.gii', 
     'sulc_rh': f'{file_path}/fs_LR.32k.R.sulc.shape.gii', 
     'atlas_rh': f'{file_path}/HCD2336346_V1_MR_hemi-R_mesh-native_dens-native_desc-real_r.gii',
     'ecc_rh': f'{file_path}/HCD2336346_V1_MR_hemi-R_mesh-native_dens-native_label-eccentricity_desc-real_roi-v2th00_metric.gii',
     'png_rh': f'{file_path}/HCD2336346_V1_MR_hemi-R_mesh-native_dens-native_label-polarangle_desc-real_roi-v2th00_metric.gii',}
]

# Step 2: Labels you want to highlight in the atlas (e.g., 1, 2, and 3)
labels_to_display = [1, 2, 3]
lim = [
    (-85, 85),
    (-85, 85),
    (-0.85, 0.85),
    (-0.85, 0.85),
    (-85, 85),
    (-85, 85),
    (-85, 85),
]
# Step 3: Create a 2x6 grid for plotting
for param in ["eccentricity", "polarangle"]:
    fig, axes = plt.subplots(6, 2, subplot_kw={'projection': '3d'}, figsize=(2.5, 8))
    left_view = [(0, 320),
                (0, 310),
                (0, 300),
                (0, 300),
                (10, 330),
                (10, 330),
    ]
    right_view = [(10, 240),
                (30, 240),
                (0, 240),
                (0, 250),
                (10, 220),
                (10, 220),
    ]
    
    name = ["tri-2 prenatal (26w)",
        "tri-3 prenatal (32w)",
        "preterm neonatal (33w/35w)",
        "full-term neonatal (40w/40w)",
        "adolescent (13.9y)",
        "adult (20.3y)"
    ]
    # Step 4: Loop through all individuals and plot their native surfaces with selected atlas labels
    for i, ind in enumerate(individuals):
        if param == "eccentricity":
            param_lh = ind['ecc_lh']
            param_rh = ind['ecc_rh']
            vmax=90
            label = "Eccentricity"
        elif param == "polarangle":
            param_lh = ind['png_lh']
            param_rh = ind['png_rh']    
            vmax=180
            label = "Polar angle"
        # Plot left hemisphere in the first row
        print(param_lh)
        ax_left = axes[i,0]
        ax_right = axes[i,1]
        surf_lh = surface.load_surf_mesh(ind['surf_lh'])
        sulc_lh = surface.load_surf_data(ind['sulc_lh'])
        brain_left = surface.load_surf_data(param_lh)
        lh_sulc_map_binary = np.where(sulc_lh > 0.0, 0.25, 0.6)
        surf_rh = surface.load_surf_mesh(ind['surf_rh'])
        sulc_rh = surface.load_surf_data(ind['sulc_rh'])
        brain_right = surface.load_surf_data(param_rh)
        rh_sulc_map_binary = np.where(sulc_rh > 0.0, 0.25, 0.6)

        # Plot eccentricity data (no individual colorbar)
        left_plot = plotting.plot_surf_stat_map(surf_lh, brain_left, hemi='left', view=left_view[i], 
                                            bg_map=lh_sulc_map_binary, axes=ax_left, figure=fig, 
                                            darkness=1.0,
                                            colorbar=False, cmap="jet", vmax=vmax, threshold=0.001)
        
        ax_left.set_xlim(lim[i])
        ax_left.set_ylim(lim[i])
        ax_left.set_zlim(lim[i])
        # Plot eccentricity data (no individual colorbar)
        right_plot = plotting.plot_surf_stat_map(surf_rh, brain_right, hemi='right', view=right_view[i], 
                                            bg_map=rh_sulc_map_binary, axes=ax_right, figure=fig, 
                                            darkness=1.0,
                                            colorbar=False, cmap="jet", vmax=vmax, threshold=0.001)
        
        ax_right.set_xlim(lim[i])
        ax_right.set_ylim(lim[i])
        ax_right.set_zlim(lim[i])

        if i == 0:
            fig.text(0.5, 1 - 0.01,
                f'{name[i]}', ha='center', va='center', fontsize=13)
        else:
            fig.text(0.5, 1 - (i * 0.145)-0.0215,
                f'{name[i]}', ha='center', va='center', fontsize=13)
            
    # Step 5: Adjust layout and display
    plt.tight_layout(rect=[0.01, 0.1, 1, 1])  # Leave space on the right for colorbars
        
    fig.text(0.06, 0.965, 'L', fontsize=12, va='center',)
    fig.text(0.92, 0.965, 'R', fontsize=12, va='center',)        
    cbar_ax = fig.add_axes([0.1, 0.08, 0.85, 0.02])
    cbar = plt.colorbar(cm.ScalarMappable(cmap='jet', norm=plt.Normalize(
        vmin=0, vmax=vmax)), cax=cbar_ax, orientation='horizontal', fraction=.1)
    cbar.set_label(f'{label}')

    # Set the ticks of the colorbar to include the min and max values
    cbar.set_ticks([0, vmax/2, vmax])

    # Set the tick labels to show min and max values (you can format them as needed)
    cbar.set_ticklabels(['0', f'{vmax/2:.0f}', f'{vmax:.0f}'])
    #plt.show()
    if param =="eccentricity":
        fig_name = "S6"
    elif param == "polarangle":
        fig_name = "S10"
    fig.savefig(f"{file_path}/Figure{fig_name}.png",
                    dpi=300, bbox_inches='tight')