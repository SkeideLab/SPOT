from pathlib import Path

import nibabel as nib
import numpy as np

from nilearn import plotting, surface

from ccfprocedures import perform_ccf_analysis
from meshgraph import calc_shortestpath


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


def save_results(root_dir, sub, ses, hemi, indices_v2, best_models, n_total_nodes):
    save_path = Path(root_dir) / "ccfmodel" / f"sub-{sub}" / f"ses-{ses}"
    save_path.mkdir(exist_ok=True, parents=True)
    if SIMULATED:
        # results from local correlations simulation
        save_file = save_path / f"sub-{sub}_ses-{ses}_hemi-{hemi}_desc-lc"
    else:
        save_file = save_path / f"sub-{sub}_ses-{ses}_hemi-{hemi}_desc-ccf"

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
        nib.save(params_img, f"{save_file}_{param}.gii")


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


def main():
    # -------------OPTIONS-----------------------------------
    global DEBUG, OPTIMIZE, VISUALIZE, SIMULATED
    DEBUG = False
    OPTIMIZE = True
    VISUALIZE = False
    # should we use original subjects data or simulated data for this subject based on
    # local spatial correlations
    SIMULATED = True
    sub = "CC00058XX09"
    ses = "11300"
    hemi = "left"
    root_dir = "/data/p_02495/dhcp_derivatives"  # bids derivatives directory
    # -------------------------------------------------------------------------

    # default values
    n_sigma = 10  # number of different sigma values for spread of CCF to try
    optimize_threshold = 0.1  # model rsquared threshold for optimization
    if SIMULATED:
        func_path = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_space-T2w_desc-simulated_bold.func.gii"
    else:
        func_path = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/func/func_hemi-{hemi}_mesh-native.func.gii"
    visparc_path = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native_desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
    wm_path = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_space-func_wm.surf.gii"
    curv_path = "{root_dir}/dhcp_anat_pipeline/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_space-T2w_curv.shape.gii"
    distancefile_path = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_space-T2w_desc-V1_dijkstra.npy"
    labels_v1 = [1, 2]  # see `ROIfiles_Labeling.txt` in this directory
    labels_v2 = [3, 4]

    # -----------------LOAD DATA-----------------------------------
    hemi_upper = hemi[0].upper()
    func = surface.load_surf_data(
        func_path.format(
            root_dir=root_dir, sub=sub, ses=ses, hemi=hemi, hemi_upper=hemi_upper
        )
    )
    visparc = nib.load(
        visparc_path.format(root_dir=root_dir, sub=sub, ses=ses, hemi_upper=hemi_upper)
    )
    wm = surface.load_surf_mesh(
        wm_path.format(root_dir=root_dir, sub=sub, ses=ses, hemi_upper=hemi_upper)
    )
    curv = curv_path.format(root_dir=root_dir, sub=sub, ses=ses, hemi_upper=hemi_upper)
    distance_file = distancefile_path.format(
        root_dir=root_dir, sub=sub, ses=ses, hemi_upper=hemi_upper
    )

    # ------------------------------PREPARE DATA---------------------------------------------------
    indices_v1 = get_indices_roi(labels_v1, visparc)
    func_v1 = func[indices_v1, :].astype(np.float64)

    indices_v2 = get_indices_roi(labels_v2, visparc)
    func_v2 = func[indices_v2, :].astype(np.float64)

    sigmas = np.linspace(3, 25, num=n_sigma)

    # get distances between nodes
    if Path(distance_file).exists():
        with open(distance_file, "rb") as f:
            distances_along_mesh = np.load(f)
    else:
        distances_along_mesh = calc_shortestpath(wm, indices_v1)
        with open(distance_file, "wb") as f:
            np.save(f, distances_along_mesh)

    # small toy example for debugging
    if DEBUG:
        func_v1 = func_v1[:10, :]
        func_v2 = func_v2[:10, :]
        indices_v1 = indices_v1[:10]
        indices_v2 = indices_v2[:10]
        distances_along_mesh = distances_along_mesh[:10, :10]
        sigmas = sigmas[[0]]
        optimize_threshold = 0

    # prepare functional data for timeseries computation
    func_v1 = make_percent_signal_change(func_v1)
    func_v2 = make_percent_signal_change(func_v2)

    # ----------------------MODELING------------------------------------------------------------
    connfields = perform_ccf_analysis(
        OPTIMIZE,
        optimize_threshold,
        distances_along_mesh,
        func_v1,
        func_v2,
        sigmas,
    )

    if VISUALIZE and not DEBUG:
        visualize_connective_field(
            wm, indices_v1, connfields["connfield_weights"], curv
        )

    save_results(
        root_dir,
        sub,
        ses,
        hemi_upper,
        indices_v2,
        connfields["best_models"],
        wm.coordinates.shape[0],
    )

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
