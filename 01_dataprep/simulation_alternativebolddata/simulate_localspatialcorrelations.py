"""Generate artificial functional data on the surface (native space).
Simulated data will have the same number of timepoints as our resting-state data. 

Data are simulated based on:

1) Draw values from a random normal distribution for each volume/vertex datapoint
2) Smooth over the surface with a 30mm Gaussian kernel
3) Add random noise with amplitude 200x that of the random normal distribution

This results in data that show a pattern of local spatial correlations."""

import subprocess
from argparse import ArgumentParser
from pathlib import Path
import os
import nibabel as nib
import numpy as np
from numpy import random

from nilearn import surface

PREFIX_FUNC = (
    "{derivatives_path}/dhcp_surface/sub-{sub}/ses-{ses}/func/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}"
)
PATH_SURFPIAL = (
    "{derivatives_path}/sub-{sub}/ses-{ses}/anat/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_pial.surf.gii"
)
PATH_FUNC = "{prefix_func}_mesh-native_bold.func.gii"
PATH_SIMIMG = "{prefix_func}_mesh-native_desc-simulated_bold.func.gii"


def parse_args():
    """Parses arguments from the command line."""

    parser = ArgumentParser()
    parser.add_argument(
        "-da",
        "--anat_directory",
        required=True,
        help="Superdirectory for top-level dhcp derivatives datasets ",
    )
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
    """Save timeseries in gifti file format.

    Args:
        data (numpy.array): n_vertices x n_timepoints. Timeseries data.
        out_path (str): full path to output file for saving.
    """

    # save data as nifti with 1 data array with all vertices per timepoint
    # list of timepoint activation arrays
    timepoint_arrays = [data[:, i_timepoint] for i_timepoint in range(data.shape[1])]

    # convert normal arrays into gifti arrays and collect in gifti image
    img = nib.gifti.GiftiImage(
        darrays=[
            nib.gifti.GiftiDataArray(
                np.float32(timepoint_array), datatype="NIFTI_TYPE_FLOAT32"
            )
            for timepoint_array in timepoint_arrays
        ]
    )
    nib.save(img, out_path)


def main():
    args = parse_args()

    for hemi in ["L", "R"]:
        if hemi == "L":
            hemi_lower = "left"
        elif hemi =="R":
            hemi_lower = "right"

        prefix_func = PREFIX_FUNC.format(
            sub=args.sub,
            ses=args.ses,
            hemi=hemi,
            derivatives_path=args.derivatives_directory,
        )
        path_thissurf = PATH_SURFPIAL.format(
            sub=args.sub,
            ses=args.ses,
            hemi=hemi_lower,
            derivatives_path=args.anat_directory,
        )
        path_thisfunc = PATH_FUNC.format(prefix_func=prefix_func)

        # save intermediary files as 'tmp'
        path_thissim = PATH_SIMIMG.format(prefix_func=prefix_func).replace(
            "bold", "bold_tmp"
        )

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

        # TODO check if we should use afni flag -sigma directly? takes ~1day just to run this!
        # TODO check resulting sigma spread and correlation results of ccf model
        sigma = 20  # target sigma from Bock 2015
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

        # save final data
        save_gifti_timeseries(sim_data, path_thissim.replace("_tmp", ""))

        print(f"Saving simulated image in {path_thissim.replace('_tmp', '')}.")

        # delete intermediary images
        for p in Path(path_thissim).parent.glob("*_tmp.*"):
            p.unlink()


#if __name__ == "__main__":
#    main()
args = parse_args()
derivatives_path = args.derivatives_directory
sub = args.sub
ses = args. ses
 #Check if the files exist
if os.path.exists(f"{derivatives_path}/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_hemi-L_mesh-native_desc-simulated_bold.func.gii") and os.path.exists(f"{derivatives_path}/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_hemi-R_mesh-native_desc-simulated_bold.func.gii"):
    print("Files already exist. Skipping the process.")
else:
    # Process to be executed if the files don't exist
    print("Files do not exist. Proceeding with the process...")
    main()