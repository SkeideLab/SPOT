from surfplot import Plot
from neuromaps.datasets import fetch_fsaverage
from nilearn import surface
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from nilearn import plotting
surfaces = fetch_fsaverage(density="164k")
lh, rh = surfaces["inflated"]
#p = Plot(lh, rh, zoom=1.4,  views=['medial', 'posterior'], size=(900,600))

# Load sulcal depth data
lh_sulc = surface.load_surf_data("/data/p_02915/templates/template_fsaverage/fsaverage/surf/lh.curv")  # Load the left hemisphere sulcal depth data
rh_sulc = surface.load_surf_data("/data/p_02915/templates/template_fsaverage/fsaverage/surf/rh.curv")  # Load the right hemisphere sulcal depth data
lh_sulc_map_binary = np.where(lh_sulc < 0.0, 0.25, 0.6)
rh_sulc_map_binary = np.where(rh_sulc > 0.0, 0.25, 0.4)
#p.add_layer({'left': lh_sulc_map_binary, 'right': rh_sulc_map_binary}, cmap='binary_r', cbar=False, color_range=(0, 1))

color_range=(0,23)

fetal_old_left = surface.load_surf_data("/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_older_fetal_L_label-eccentricity_desc-real_roi-v2th00_metric.gii")
fetal_old_right = surface.load_surf_data("/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/Averaged_older_fetal_R_label-eccentricity_desc-real_roi-v2th00_metric.gii")
#p.add_layer({"left": fetal_old_left, "right": fetal_old_right}, cmap="jet", cbar=True )

cbar_kws = dict(outer_labels_only=True, pad=.02, n_ticks=2, decimals=0)
#fig = p.build(cbar_kws=cbar_kws)
# add units to colorbar

#fig.savefig("/data/p_02915/dhcp_derivatives_SPOT/Figures_v2/Averaged_older_fetal_eccentricity.png", dpi=300)
#plt.show()

fig = plotting.plot_surf_stat_map(
    lh, fetal_old_left, hemi="left", view=(-10, -70.0), darkness=1.0,
    colorbar=True, cmap="jet",
    threshold=0.1, bg_map=lh_sulc_map_binary,
)
fig.savefig("/data/p_02915/dhcp_derivatives_SPOT/Figures_v2/Averaged_older_fetal_eccentricity_V2.png", dpi=300)
plt.show()