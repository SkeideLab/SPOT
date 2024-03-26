# alternative model of local spatial correlations
# in subject space for now
import subprocess
from argparse import ArgumentParser
from pathlib import Path

import nibabel as nib
import numpy as np
from numpy import random

from nilearn import surface

PREFIX_FUNC = "{derivatives_path}/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_hemi-{hemi}"
PATH_SURFPIAL = "{derivatives_path}/dhcp_anat_pipeline/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_space-T2w_pial.surf.gii"
PATH_FUNC = "{prefix_func}_mesh-native_bold.func.gii"
PATH_SIMIMG = "{prefix_func}_space-T2w_desc-simulated_bold.func.gii"


def parse_args():
    """Parses arguments from the command line."""

    parser = ArgumentParser()
    parser.add_argument(
        "-d",
        "--derivatives_directory",
        required=True,
        help="Superdirectory for top-level dhcp derivatives datasets ",
    )
    parser.add_argument(
        "-sub",
        required=True,
        help="Subject being processed",
    )
    parser.add_argument(
        "-ses",
        required=True,
        help="Session being processed",
    )
    args = parser.parse_args()
    return args


def save_gifti_timeseries(data, out_path):
    # n_nodes x n_timepoints
    # save it as nifti with 1 data array per timepoint
    gifti_arrays = [
        nib.gifti.GiftiDataArray(
            np.float32(data[:, i]),
            datatype="NIFTI_TYPE_FLOAT32",
        )
        for i in range(data.shape[1])
    ]
    img = nib.gifti.GiftiImage(darrays=gifti_arrays)
    nib.save(img, out_path)


def main():
    args = parse_args()

    for hemi in ["L", "R"]:

        prefix_func = PREFIX_FUNC.format(
            sub=args.sub,
            ses=args.ses,
            hemi=hemi,
            derivatives_path=args.derivatives_directory,
        )
        path_thissurf = PATH_SURFPIAL.format(
            sub=args.sub,
            ses=args.ses,
            hemi=hemi,
            derivatives_path=args.derivatives_directory,
        )
        path_thisfunc = PATH_FUNC.format(prefix_func)
        path_thissim = PATH_SIMIMG.format(prefix_func).replace("bold", "bold_tmp")

        # load surface and functional data
        surf = surface.load_surf_mesh(path_thissurf)
        func = surface.load_surf_data(path_thisfunc)

        # extract number of nodes and timepoints from it
        # everything we need for the simulation
        n_nodes = surf.coordinates.shape[0]
        n_timepoints = func.shape[1]

        #################################SIMULATING NOISE##############################
        print("Simulating data on the surface...")
        # simulate data with mean 0, sigma 1, for all timepoints and nodes
        sim_data = random.default_rng().normal(size=(n_nodes, n_timepoints))

        save_gifti_timeseries(sim_data, path_thissim)

        ##############################PRODUCING LOCAL CORRELATIONS #######################################################
        print("Smoothing the data over the surface to produce local correlations...")

        sigma = 20
        fwhm = 2.354820 * sigma
        path_intermed = path_thissim.replace(".func", ".func_smh7") + ".dset"
        path_smoothed = path_thissim.replace("desc-simulated", "desc-simulatedsmoothed")

        # smooth surface with AFNI
        cmd_surfsmooth = f"AFNI SurfSmooth -i {path_thissurf} -met HEAT_07 -input {path_thissim} -target_fwhm {fwhm:.0f}"
        cmd_convertdset = (
            f"AFNI ConvertDset -o_gii -input {path_intermed} -prefix {path_smoothed}"
        )
        subprocess.run(cmd_surfsmooth, shell=True, check=True)
        subprocess.run(cmd_convertdset, shell=True, check=True)

        ###################################ADDING NOISE#############################
        print("Adding random noise to data...")

        # load smoothed image
        sim_data = surface.load_surf_data(path_smoothed)

        # add noise with amplitude 200
        sim_data = sim_data + random.default_rng().normal(size=sim_data.shape) * 200

        # save it as nifti with 1 data array per timepoint
        save_gifti_timeseries(sim_data, path_thissim.replace("_tmp", ""))

        print(f"Saving simulated image in {path_thissim.replace('_tmp', '')}.")

        # delete intermediary images
        for p in Path(path_thissim).parent.glob("*_tmp.*"):
            p.unlink()
