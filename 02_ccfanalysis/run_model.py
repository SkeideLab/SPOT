import importlib
import sys
from argparse import ArgumentParser
from pathlib import Path

import nibabel as nib
import numpy as np

from nilearn import plotting, surface

# HACK
# ---------------------------import workarounds to make repo structure work --------------------------#
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))

# Additionally remove the current file's directory from sys.path
try:
    sys.path.remove(str(parent))
except ValueError:  # Already removed
    pass

ccfanalysis = importlib.import_module("02_ccfanalysis")
perform_ccf_analysis = ccfanalysis.ccf_model.ccfprocedures.perform_ccf_analysis
calc_shortestpath = ccfanalysis.ccf_model.meshgraph.calc_shortestpath
# HACKEND
# ------------------------------------------------------------------------#

VISPARC_PATH = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
WM_PATH = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_space-bold_wm.surf.gii"
CURV_PATH = "{root_dir}/dhcp_anat_pipeline/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_space-T2w_curv.shape.gii"
DISTANCEFILE_PATH = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}_space-T2w_desc-V1_dijkstra.npy"
func_path = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native{simulated}_bold.func.gii"
LABELS_V1 = (1, 2)
LABELS_V2 = (3, 4)


def visualize_connective_field(mesh, v1_indices, connective_fields, curv):
    voxel_to_plot = 1000
    sigma_to_plot = 0
    cf_voxelseries = connective_fields[voxel_to_plot, sigma_to_plot, :]
    cf_map = np.zeros(mesh.coordinates.shape[0])
    cf_map[v1_indices] = cf_voxelseries

    plotting.plot_surf(
        mesh,
        surf_map=cf_map,
        view="posterior",
        bg_map=curv,
        threshold=0.01,
        cmap="rainbow",
    )


def make_percent_signal_change(func):
    """Compute activation in percent signal change from intensity.

    Args:
        func (numpy.array): n_vertices(n_voxels) x n_timepoints. Raw signal intensity.

    Returns:
        numpy.array: n_vertices(n_voxels) x n_timepoints.
            Activation in units of percent signal change.
    """
    return (func.T - np.mean(func, axis=1)).T * 100


def save_results(
    root_dir, sub, ses, hemi, indices_v2, best_models, n_total_nodes, model
):
    save_path = (
        Path(root_dir)
        / "ccfmodel"
        / f"sub-{sub}"
        / f"ses-{ses}"
        / f"sub-{sub}_ses-{ses}_hemi-{hemi}_desc-{model}"
    )
    save_path.parent.mkdir(exist_ok=True, parents=True)

    param_full_mesh = np.zeros(n_total_nodes)
    parameters = ("v0i", "sigma", "rss", "rsquared")
    for i, param in enumerate(parameters):
        param_full_mesh[indices_v2] = best_models[:, i]
        if param == "v0i":
            params_img = nib.gifti.GiftiImage(
                darrays=[
                    nib.gifti.GiftiDataArray(
                        np.int32(param_full_mesh), intent="NIFTI_INTENT_LABEL"
                    )
                ]
            )
        else:
            params_img = nib.gifti.GiftiImage(
                darrays=[nib.gifti.GiftiDataArray(np.float32(param_full_mesh))]
            )
        nib.save(params_img, f"{save_path}_{param}.gii")


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
    )
    parser.add_argument(
        "--optimize",
        required=False,
        default=True,
        help="Run nonlinear search for best sigma spread value of each vertex ccf after grid search",
    )

    args = parser.parse_args()
    return args


def fetch_distancematrix(indices_v1, wm, distance_file):

    if Path(distance_file).exists():

        with open(distance_file, "rb") as f:
            distances_along_mesh = np.load(f)

    else:

        print("Distances between surface vertices are being calculated...")
        distances_along_mesh = calc_shortestpath(wm, indices_v1)
        with open(distance_file, "wb") as f:
            np.save(f, distances_along_mesh)

    return distances_along_mesh


def main():

    args = parse_args()

    # sub = "CC00069XX12"
    # ses = "26300"
    # root_dir = "/data/p_02495/dhcp_derivatives"  # bids derivatives directory

    optimize_threshold = args.threshold  # model rsquared threshold for optimization
    sigmas = np.linspace(3, 25, num=args.nsigma)  # ccf spreads to try in grid search

    # model both data sources (resting-state and simulated) per hemisphere
    for hemi in ["L", "R"]:

        ids = {
            "sub": args.sub,
            "ses": args.ses,
            "hemi": hemi,
            "root_dir": args.root_dir,
        }

        # load retinotopy data
        visparc = nib.load(VISPARC_PATH.format(**ids))
        indices_v1 = get_indices_roi(LABELS_V1, visparc)
        indices_v2 = get_indices_roi(LABELS_V2, visparc)

        # load surface mesh data
        wm = surface.load_surf_mesh(WM_PATH.format(**ids))
        curv = CURV_PATH.format(**ids)

        # load distance between surface vertices
        distance_file = DISTANCEFILE_PATH.format(**ids)
        distances_along_mesh = fetch_distancematrix(indices_v1, wm, distance_file)

        for model in ["ccf", "lc"]:

            # get bold data for V1 and V2
            simulated = "_desc-simulated" if model == "lc" else ""
            func = surface.load_surf_data(func_path.format(**ids, simulated=simulated))
            func_v2 = func[indices_v2, :].astype(np.float64)
            func_v1 = func[indices_v1, :].astype(np.float64)

            # small toy example for debugging
            if args.debug:
                func_v1 = func_v1[:10, :]
                func_v2 = func_v2[:10, :]
                indices_v1 = indices_v1[:10]
                indices_v2 = indices_v2[:10]
                distances_along_mesh = distances_along_mesh[:10, :10]
                sigmas = sigmas[[0]]
                optimize_threshold = 0

            # ----------------------MODELING------------------------------------------------------------
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

            if not args.debug:
                visualize_connective_field(
                    wm, indices_v1, connfields["connfield_weights"], curv
                )

            save_results(
                **ids,
                indices_v2=indices_v2,
                best_models=connfields["best_models"],
                n_total_nodes=wm.coordinates.shape[0],
                model=model,
            )


# TODO transform this info into info text string
# g: connective field of V2 seed voxel, v: voxel within V1,
# v_0: center voxel of connective field, d: shortest distance along the cortical surface mesh
# sigma: SD (mm) along cortical surface
# g(v_0,sigma) = exp - [d(v,v_0)²/ 2*sigma²]

# -------------------------Questions-----------------------------------
# what is SD along cortical surface?
# --> " For each v0 value, we created 10 basis functions, using σ values linearly spaced between 3 and 25 mm"
# which cortical surface mesh? midthickness, white, pial?
# --> grey/white matter border so maybe white
# how to understand the exp and sigma² in equation?
# What is v_0?
# --> "center v0 corresponding to the surface vertex closest to that V1 voxel at the gray/white matter border"
# how to get predicted time course?
# --> "taking the linear sum of all V1 time courses convolved with all possible basis functions"
# convolve connective field g(v_0,sigma) for different v0 and sigma with v timecourses + sum up linearly
# -------------------------Questions end------------------------------------


if __name__ == "__main__":
    main()
