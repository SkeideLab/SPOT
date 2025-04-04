"""Runs cortical connective field modelling. 
For exact model description, see SPOT/documentation/CorticalConnectiveFields.pdf.

# v_0: center voxel of connective field in V1
# sigma: SD (mm) along cortical surface

"""

import importlib
import sys
from argparse import ArgumentParser
from pathlib import Path

import nibabel as nib
import numpy as np

from nilearn import surface

# -d /data/p_02495/dhcp_derivatives -sub CC00069XX12 -ses 26300
# HACK
# ---------------------import workarounds to make repo structure work --------------------------#
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))

# Additionally remove the current file's directory from sys.path
try:
    sys.path.remove(str(parent))
except ValueError:  # Already removed
    pass

ccfanalysis = importlib.import_module("00_HCP")
perform_ccf_analysis = ccfanalysis.ccf_model.ccfprocedures.perform_ccf_analysis
calc_shortestpath = ccfanalysis.ccf_model.meshgraph.calc_shortestpath
# HACKEND
# ------------------------------------------------------------------------#

VISPARC_PATH = (
    "/data/p_02915/SPOT/00_HCP/template/"
    "hemi-{hemi}_mesh-native_dens-native_"
    "desc-retinotbenson2014_label-visarea_dparc.label.gii"
)
WM_PATH = (
    "/data/pt_02880/HCP_D/fmriresults01/{sub}/MNINonLinear/fsaverage_LR32k/"
    "{sub}.{hemi}.white_MSMAll.32k_fs_LR.surf.gii"
)

DISTANCEFILE_PATH = (
    "{root_dir}/ccfmodel/{sub}/"
    "{sub}_hemi-{hemi}_space-T2w_desc-V1_dijkstra.npy"
)
FUNC_PATH = (
    "{root_dir}/hcp_surface/{sub}/func/"
    "{sub}_hemi-{hemi}-Atlas-MSMAll_hp0_clean_bold{simulated}.func.gii"
)
OUTPUT_PREFIX = (
    "{root_dir}/ccfmodel/{sub}/"
    "{sub}_hemi-{hemi}_mesh-native_dens-native_desc-{datasource}"
)
LABELS_V1 = [1]
LABELS_V2 = [2, 3]


def make_percent_signal_change(func):
    """Compute activation in percent signal change from intensity.

    Args:
        func (numpy.array): n_vertices(n_voxels) x n_timepoints. Raw signal intensity.

    Returns:
        numpy.array: n_vertices(n_voxels) x n_timepoints.
            Activation in units of percent signal change.
    """
    return (func.T - np.mean(func, axis=1)).T * 100


def save_results(indices_v2, best_models, n_total_nodes, out_prefix):
    """Save model results to gifti files.

    Args:
        indices_v2 (numpy.array): indices of V2 vertices in whole hemisphere
        best_models (numpy.array): n_vertices_v2 x (v0i, sigma, rss, r)
        n_total_nodes (int): size of total hemispheres
        out_prefix (str): path prefix for result files
    """

    Path(out_prefix).parent.mkdir(exist_ok=True, parents=True)

    # empty mesh of whole hemi
    # will be filled at V2 indices
    param_full_mesh = np.zeros(n_total_nodes)
    parameters = ("v0i", "sigma", "rss", "r")

    for i, param in enumerate(parameters):

        param_full_mesh[indices_v2] = best_models[:, i]

        if param == "v0i":
            darray = nib.gifti.GiftiDataArray(np.int32(param_full_mesh))
        else:
            darray = nib.gifti.GiftiDataArray(np.float32(param_full_mesh))

        params_img = nib.gifti.GiftiImage(darrays=[darray])
        print(f"Saving results in {out_prefix}_{param}.gii...")
        nib.save(params_img, f"{out_prefix}_{param}.gii")


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
        "-nsigma",
        required=False,
        default=10,
        help="Number of different ccf spreads to try (between 3 and 25mm)",
    )
    parser.add_argument(
        "-th",
        "--threshold",
        required=False,
        default=0.1,
        help=(
            "Correlation threshold for acceptable models. Vertices that have a "
            "correlation >th with their best model after grid search will go on to "
            "optimization (if enabled). Other vertices' results will be discarded."
        ),
    )
    parser.add_argument(
        "--debug",
        required=False,
        default=False,
        help=(
            "Use only a few vertices from the input dataset and 1 sigma value "
            "for quick execution during debug"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--no-optimization",
        dest="optimize",
        action="store_false",
        help="Don't run nonlinear search for best sigma spread value of each vertex ccf after grid search",
    )

    args = parser.parse_args()
    return args


def fetch_distancematrix(indices_v1, wm, distance_file):
    """Get the distances between nodes in V1 along the mesh. Either load it or compute it.

    Args:
        indices_v1 (numpy.array): indices to V1 nodes
        wm (namedtuple): keys coordinates and faces
        distance_file (str): path to distance file if it exists

    Returns:
        numpy.array: distances between nodes
    """

    if Path(distance_file).exists():

        with open(distance_file, "rb") as f:
            distances_along_mesh = np.load(f)

    else:

        print("Distances between surface vertices are being calculated...")
        distances_along_mesh = calc_shortestpath(wm, indices_v1)

        Path(distance_file).parent.mkdir(exist_ok=True, parents=True)
        with open(distance_file, "wb") as f:
            np.save(f, distances_along_mesh)

    return distances_along_mesh


def main():

    args = parse_args()

    optimize_threshold = args.threshold  # model rsquared threshold for optimization
    # ccf spreads to try in grid search
    sigmas = np.linspace(3, 25, num=args.nsigma)

    # model both data sources (resting-state and simulated) per hemisphere
    for hemi in ["L", "R"]:

        ids = {
            "sub": args.sub,
            "hemi": hemi,
            "root_dir": args.derivatives_directory,
        }

        # load retinotopy data
        visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
        indices_v1 = get_indices_roi(LABELS_V1, visparc)
        indices_v2 = get_indices_roi(LABELS_V2, visparc)

        # load surface mesh data
        wm = surface.load_surf_mesh(WM_PATH.format(sub=args.sub, hemi=hemi))

        # load distance between surface vertices
        distance_file = DISTANCEFILE_PATH.format(**ids)
        distances_along_mesh = fetch_distancematrix(
            indices_v1, wm, distance_file)

        # small toy example for debugging
        if args.debug:
            indices_v1 = indices_v1[:10]
            indices_v2 = indices_v2[:10]
            distances_along_mesh = distances_along_mesh[:10, :10]
            sigmas = sigmas[[0]]
            optimize_threshold = 0

        for datasource in ["real", "simulated"]:
            print(f"Modeling {datasource} data on hemisphere {hemi}.")

            # get bold data for V1 and V2
            simulated = "_simulated" if datasource == "simulated" else ""

            func = surface.load_surf_data(
                FUNC_PATH.format(**ids, simulated=simulated))

            func_v2 = func[indices_v2, :].astype(np.float64)
            func_v1 = func[indices_v1, :].astype(np.float64)
            # ----------------------MODELING-------------------------------------------------------
            # prepare functional data for timeseries computation
            func_v1 = make_percent_signal_change(func_v1)
            func_v2 = make_percent_signal_change(func_v2)
            connfields = perform_ccf_analysis(
                args.optimize,
                optimize_threshold,
                distances_along_mesh,
                func_v1,
                func_v2,
                sigmas,
            )

            save_results(
                indices_v2=indices_v2,
                best_models=connfields["best_models"],
                n_total_nodes=wm.coordinates.shape[0],
                out_prefix=OUTPUT_PREFIX.format(**ids, datasource=datasource),
            )


if __name__ == "__main__":
    main()
