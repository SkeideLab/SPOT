"""Generate artificial functional data on the surface (native space).
Simulated data will have the same number of timepoints as our resting-state data. 

Data are simulated based on:

1) Draw values from a random normal distribution for each volume/vertex datapoint
2) Smooth over the surface with a 30mm Gaussian kernel
3) Add random noise with amplitude 200x that of the random normal distribution

This results in data that show a pattern of local spatial correlations."""
# python /data/p_02915/SPOT/01_dataprep/simulation_alternativebolddata/simulated_model.py -da /data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_anat_pipeline -d /data/p_02915/dhcp_derivatives_SPOT -sub CC00060XX03 -ses 12501
import subprocess
from argparse import ArgumentParser
from pathlib import Path
import os
import nibabel as nib
import numpy as np
from numpy import random
from scipy.ndimage import gaussian_filter
from scipy.stats import norm
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
VISPARC_PATH = (
    "{derivatives_path}/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_"
    "desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
)

LABELS_V1 = (1, 2)
LABELS_V2 = (3, 4)

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
        "-da",
        "--anat_directory",
        required=True,
        help="Superdirectory for dhcp anatomical datasets ",
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

def main():
    args = parse_args()
    anat_directory=args.anat_directory
    derivatives_directory=args.derivatives_directory
    sub=args.sub
    ses=args.ses

    for hemi in ["L", "R"]:
        if hemi == "L":
            hemi_lower = "left"
        elif hemi =="R":
            hemi_lower = "right"

        prefix_func = PREFIX_FUNC.format(
            sub=sub,
            ses=ses,
            hemi=hemi,
            derivatives_path=derivatives_directory,
        )
        path_thissurf = PATH_SURFPIAL.format(
            sub=sub,
            ses=ses,
            hemi=hemi_lower,
            derivatives_path=anat_directory,
        )
        path_thisfunc = PATH_FUNC.format(prefix_func=prefix_func)

        # save intermediary files as 'tmp'
        path_thissim = PATH_SIMIMG.format(prefix_func=prefix_func)

        path_visparc=VISPARC_PATH.format(
            sub=sub,
            ses=ses,
            hemi=hemi,
            derivatives_path=derivatives_directory,
            )

        # load surface and functional data
        surf = surface.load_surf_mesh(path_thissurf)
        func_real = surface.load_surf_data(path_thisfunc)

        # extract number of nodes and timepoints from it
        # everything we need for the simulation
        n_nodes = surf.coordinates.shape[0]
        n_timepoints = func_real.shape[1]
        visparc = nib.load(path_visparc)
        indices_v1 = get_indices_roi(LABELS_V1, visparc)
        indices_v2 = get_indices_roi(LABELS_V2, visparc)

        #################################SIMULATING NOISE##############################
        print("Simulating data on the surface...")
        
        # simulate data with mean 0, sigma 1, for all timepoints and nodes
        sim_data = random.default_rng().normal(size=(n_nodes, n_timepoints))
        func_v2 = sim_data[indices_v2, :].astype(np.float64)
        func_v1 = sim_data[indices_v1, :].astype(np.float64)
        func = np.vstack((func_v1, func_v2))
        real_func = np.vstack((func_real[indices_v1, :].astype(np.float64), func_real[indices_v2, :].astype(np.float64)))    


        # Step 2: Compute the correlation matrix
        correlation_matrix = np.corrcoef(func)

        # Step 3: Apply Gaussian filter to the correlation matrix
        sigma = 30  # Standard deviation for Gaussian kernel (30 mm)
        smoothed_correlation_matrix = gaussian_filter(correlation_matrix, sigma=sigma) 
        smoothed_correlation_matrix = smoothed_correlation_matrix + np.eye(smoothed_correlation_matrix.shape[0]) * 1e-6

        # Step 4: Adjust correlations to match the desired mean and variance
        desired_mean = np.median(np.corrcoef(real_func))  # Desired mean correlation
        current_mean = np.median(smoothed_correlation_matrix)
        smoothed_correlation_matrix += (desired_mean - current_mean)

        # Step 6: Perform Cholesky decomposition
        cholesky_decomposition = np.linalg.cholesky(smoothed_correlation_matrix)

        # Step 7: Generate modified time courses
        modified_time_courses = cholesky_decomposition.T @ func

        # Step 8: Update sim_data with the modified time courses
        sim_data[indices_v1, :] = modified_time_courses[:len(indices_v1), :]
        sim_data[indices_v2, :] = modified_time_courses[len(indices_v1):, :]

        save_gifti_timeseries(sim_data, path_thissim)


if __name__ == "__main__":
    main()

# /data/p_02915/dhcp_derivatives_SPOT sub-CC00056XX07 ses-10700