import numpy as np
import nibabel as nib
from nilearn import plotting
from nilearn import surface
import plotly.io as pio


mesh = surface.load_surf_mesh("/data/p_02915/templates/template_corticalsurfaceneonatessym_williams2023_dhcp/dhcpSym_template/week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii")
data = surface.load_surf_data("/data/p_02915/dhcp_derivatives_SPOT/dhcp_surface/sub-CC00056XX07/ses-10700/anat/sub-CC00056XX07_ses-10700_hemi-L_mesh-native_dens-native_desc-eccentretinotbenson2014_seg.shape.gii")
bg_map = surface.load_surf_data("/data/p_02915/templates/template_corticalsurfaceneonatessym_williams2023_dhcp/dhcpSym_template/week-40_hemi-left_space-dhcpSym_dens-32k_sulc.shape.gii")
#curv_right = surface.load_surf_data("/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_anat_pipeline/sub-CC00056XX07/ses-10700/anat/sub-CC00056XX07_ses-10700_hemi-left_curv.shape.gii")
#curv_right_sign = np.sign(curv_right)

fig = plotting.plot_surf(
        mesh,
        surf_map=bg_map,
        view="medial",
        engine='plotly',
        symmetric_cmap="True",        
        cmap='gray',    
        title="Left hemisphere",
        vmax= 7,
        )
plotting.displays.PlotlySurfaceFigure.savefig(fig, output_file='/data/p_02915/SPOT/plotly_plot.png')
#pio.write_html(fig.figure, file='/data/p_02915/SPOT/plotly_plot.html', auto_open=False)
# https://nilearn.github.io/dev/auto_examples/01_plotting/plot_3d_map_to_surface_projection.html#making-a-surface-plot-of-a-3d-statistical-map
# https://nilearn.github.io/dev/plotting/index.html#surface-plotting