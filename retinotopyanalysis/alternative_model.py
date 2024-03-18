# alternative model of local spatial correlations
# in subject space for now
from nilearn import surface
from numpy import random
import nibabel as nib
import numpy as np
import subprocess
from pathlib import Path

sub = "CC00058XX09"
ses = "11300"
hemi = "L"
hemi_long = "left"
derivatives_path = "/data/p_02495/dhcp_derivatives"

surf_path = "{derivatives_path}/dhcp_anat_pipeline/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_space-T2w_pial.surf.gii"
func_path = "{derivatives_path}/dhcp_surface/sub-{sub}/ses-{ses}/space-func/func/func_hemi-{hemi_long}_mesh-native.func.gii"
sim_img_path = "{derivatives_path}/dhcp_surface/sub-{sub}/ses-{ses}/space-func/func/sub-{sub}_ses-{ses}_hemi-{hemi}_space-T2w_desc-simulated_bold.func.gii"

print("Simulating data on the surface...")

# load surface and functional data
path_thissim = sim_img_path.format(
    sub=sub, ses=ses, hemi=hemi, derivatives_path=derivatives_path
).replace("bold", "bold_tmp")
path_thissurf = surf_path.format(
    sub=sub, ses=ses, hemi=hemi, derivatives_path=derivatives_path
)
path_thisfunc = func_path.format(
    sub=sub, ses=ses, hemi_long=hemi_long, derivatives_path=derivatives_path
)
surf = surface.load_surf_mesh(path_thissurf)
func = surface.load_surf_data(path_thisfunc)


# extract number of nodes and timepoints from it
# everything we need for the simulation
n_nodes = surf.coordinates.shape[0]
n_timepoints = func.shape[1]

# simulate data with mean 0, sigma 1, for all timepoints and nodes
sim_data = random.default_rng().normal(size=(n_nodes, n_timepoints))

# save it as nifti with 1 data array per timepoint
sim_img = nib.gifti.GiftiImage(
    darrays=[
        nib.gifti.GiftiDataArray(
            np.float32(sim_data[:, i]),
            datatype="NIFTI_TYPE_FLOAT32",
        )
        for i in range(sim_data.shape[1])
    ]
)

nib.save(
    sim_img,
    path_thissim,
)

sigma = 30
fwhm = 2.354820 * sigma
# smooth surface with AFNI
cmd_surfsmooth = f"AFNI SurfSmooth -i {path_thissurf} -met HEAT_05 -input {path_thissim} -fwhm {fwhm:.0f}"
path_intermed = path_thissim.replace(".func", ".func_sm") + ".dset"
path_smoothed = path_thissim.replace("desc-simulated", "desc-simulatedsmoothed")
cmd_convertdset = (
    f"AFNI ConvertDset -o_gii -input {path_intermed} -prefix {path_smoothed}"
)

print("Smoothing the data over the surface...")
subprocess.run(cmd_surfsmooth, shell=True, check=True)
subprocess.run(cmd_convertdset, shell=True, check=True)

print("Adding random noise to data...")
# load smoothed image
sim_data = surface.load_surf_data(path_smoothed)

# add noise with amplitude 200
sim_data = sim_data + random.default_rng().normal(size=sim_data.shape) * 200

# save it as nifti with 1 data array per timepoint
sim_img = nib.gifti.GiftiImage(
    darrays=[
        nib.gifti.GiftiDataArray(
            np.float32(sim_data[:, i]),
            datatype="NIFTI_TYPE_FLOAT32",
        )
        for i in range(sim_data.shape[1])
    ]
)
nib.save(
    sim_img,
    path_thissim.replace("_tmp", ""),
)
print(f"Saving simulated image in {path_thissim.replace('_tmp', '')}.")

# delete intermediary images
for p in Path(path_thissim).parent.glob("*_tmp.*"):
    p.unlink()
