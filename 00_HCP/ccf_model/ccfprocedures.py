import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import shared_memory
from scipy.optimize import minimize
from .ccflogic import (
    create_cf_fit,
    define_connfield_candidates,
    create_cf_timecourse,
    calc_correlation,
)

# ======================================================================
# ------------------------ Shared Memory Utilities ---------------------
# ======================================================================
def _to_shared(arr):
    """Create shared memory for a numpy array and copy data into it."""
    # Note: arr.nbytes automatically calculates size based on dtype (now float64)
    shm = shared_memory.SharedMemory(create=True, size=arr.nbytes)
    shared_arr = np.ndarray(arr.shape, dtype=arr.dtype, buffer=shm.buf)
    shared_arr[:] = arr[:]
    return shm, arr.shape, arr.dtype


def _from_shared(name, shape, dtype):
    """Access a numpy array in existing shared memory."""
    existing_shm = shared_memory.SharedMemory(name=name)
    return np.ndarray(shape, dtype=dtype, buffer=existing_shm.buf), existing_shm


# ======================================================================
# -------------------------- Helper Functions --------------------------
# ======================================================================
def gridsearch_connfield(timeseries_targetvoxel, timeseries_v0):
    """Grid search for best-fitting connective field."""
    model_rss = [None, None, -np.inf, None]
    for v0_i in range(timeseries_v0.shape[0]):
        for sigma_i in range(timeseries_v0.shape[1]):
            modelfit = calc_correlation(timeseries_targetvoxel, timeseries_v0[v0_i, sigma_i, :])
            if modelfit > model_rss[2]:
                model_rss = [v0_i, sigma_i, modelfit, modelfit**2]
    return np.array(model_rss)


def optimize_connfield(v0_i, sigma_init, distances, timeseries_sources, timeseries_target):
    """Nonlinear optimization of sigma parameter for one target voxel."""
    res = minimize(
        create_cf_fit,
        (sigma_init,),
        args=(timeseries_target, distances[[int(v0_i)], :], timeseries_sources),
        bounds=[(0.001, None)]
    )
    sigma = res.x[0]
    rss = res.fun
    fit_results = create_cf_fit(
        (sigma,),
        timeseries_target,
        distances[[int(v0_i)], :],
        timeseries_sources,
        return_modelfit=True,
    )
    rsquared = fit_results["rsquared"]
    return (v0_i, sigma, rss, rsquared)


# ======================================================================
# ------------------ Worker: Standard (non-cross) CCF ------------------
# ======================================================================
def _fit_target_voxel_shared(
    target_idx,
    y_target,
    sigmas,
    optimize,
    optimize_threshold,
    shm_meta,
):
    """Worker for one target voxel using shared-memory arrays."""
    # Attach to shared arrays
    shm_names, shapes, dtypes = shm_meta["names"], shm_meta["shapes"], shm_meta["dtypes"]

    timeseries_v0, shm_v0 = _from_shared(shm_names["timeseries_v0"], shapes["timeseries_v0"], dtypes["timeseries_v0"])
    distances, shm_dist = _from_shared(shm_names["distances"], shapes["distances"], dtypes["distances"])
    sources, shm_sources = _from_shared(shm_names["sources"], shapes["sources"], dtypes["sources"])

    # Grid search
    best_model = gridsearch_connfield(y_target, timeseries_v0)
    best_model[1] = sigmas[int(best_model[1])]  # replace sigma index with value

    # Optional optimization
    if optimize:
        rsquared = best_model[3]
        # if rsquared > optimize_threshold:
        best_model = optimize_connfield(
            best_model[0],
            best_model[1],
            distances,
            sources,
            y_target,
        )

    # Close shared handles
    shm_v0.close()
    shm_dist.close()
    shm_sources.close()

    return target_idx, best_model


# ======================================================================
# --------------------------- MAIN ANALYSIS ----------------------------
# ======================================================================
def perform_ccf_analysis(
    OPTIMIZE,
    optimize_threshold,
    distances_along_mesh,
    timeseries_sources,
    timeseries_targets,
    sigmas,
    n_jobs=30,
):
    """
    Parallel + shared-memory cortical connective field analysis.
    UPDATED: Uses float64 (default) to prevent precision loss.
    """
    # --- 1. Candidate CFs
    connfield_weights = define_connfield_candidates(distances_along_mesh, sigmas)

    # --- 2. Precompute CF timecourses
    n_v0, n_sigmas, n_timepoints = connfield_weights.shape[0], connfield_weights.shape[1], timeseries_sources.shape[1]
    
    # CHANGED: dtype=np.float64 (was float32)
    timeseries_v0 = np.empty((n_v0, n_sigmas, n_timepoints), dtype=np.float64)
    
    for v0_i in range(n_v0):
        for sigma_i in range(n_sigmas):
            timeseries_v0[v0_i, sigma_i, :] = create_cf_timecourse(
                timeseries_sources, connfield_weights[v0_i, sigma_i, :]
            )

    # --- 3. Create shared memory
    # CHANGED: Removed .astype(np.float32). Now uses default (likely float64).
    shm_v0, shape_v0, dtype_v0 = _to_shared(timeseries_v0)
    shm_dist, shape_dist, dtype_dist = _to_shared(distances_along_mesh) 
    shm_sources, shape_sources, dtype_sources = _to_shared(timeseries_sources)

    shm_meta = {
        "names": {"timeseries_v0": shm_v0.name, "distances": shm_dist.name, "sources": shm_sources.name},
        "shapes": {"timeseries_v0": shape_v0, "distances": shape_dist, "sources": shape_sources},
        "dtypes": {"timeseries_v0": dtype_v0, "distances": dtype_dist, "sources": dtype_sources},
    }

    # --- 4. Parallel loop
    n_targets = timeseries_targets.shape[0]
    best_models = np.zeros((n_targets, 4))

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = [
            executor.submit(
                _fit_target_voxel_shared,
                i,
                timeseries_targets[i, :],
                sigmas,
                OPTIMIZE,
                optimize_threshold,
                shm_meta,
            )
            for i in range(n_targets)
        ]
        for f in as_completed(futures):
            idx, model = f.result()
            best_models[idx, :] = model

    # --- 5. Cleanup
    shm_v0.close(); shm_v0.unlink()
    shm_dist.close(); shm_dist.unlink()
    shm_sources.close(); shm_sources.unlink()

    return {"best_models": best_models, "connfield_weights": connfield_weights}


# ======================================================================
# ---------------- Worker: Cross-Validation Analysis -------------------
# ======================================================================
def _fit_target_voxel_cross_shared(
    target_idx,
    y_first,
    y_second,
    sigmas,
    optimize,
    optimize_threshold,
    shm_meta,
):
    """Worker for one target voxel in cross-validation CCF."""
    # Attach to shared arrays
    shm_names, shapes, dtypes = shm_meta["names"], shm_meta["shapes"], shm_meta["dtypes"]
    timeseries_v0_first, shm_v0 = _from_shared(shm_names["timeseries_v0_first"], shapes["timeseries_v0_first"], dtypes["timeseries_v0_first"])
    distances, shm_dist = _from_shared(shm_names["distances"], shapes["distances"], dtypes["distances"])
    sources_first, shm_sources_first = _from_shared(shm_names["sources_first"], shapes["sources_first"], dtypes["sources_first"])
    sources_second, shm_sources_second = _from_shared(shm_names["sources_second"], shapes["sources_second"], dtypes["sources_second"])

    # --- GRID SEARCH on first half ---
    v0_idx, sigma_idx, rss_first, r2_first = gridsearch_connfield(y_first, timeseries_v0_first)
    best_v0 = int(v0_idx)
    best_sigma = sigmas[int(sigma_idx)]

    # --- OPTIONAL OPTIMIZATION ---
    res = minimize(
        create_cf_fit,
        [best_sigma],
        args=(y_first, distances[[best_v0], :], sources_first),
        bounds=[(0.001, None)]
    )
    best_sigma = res.x[0]
    rss_first = res.fun
    fit_first = create_cf_fit(
        [best_sigma],
        y_first,
        distances[[best_v0], :],
        sources_first,
        return_modelfit=True,
    )
    r2_first = fit_first["rsquared"]

    # --- EVALUATE SECOND HALF ---
    fit_second = create_cf_fit(
        [best_sigma],
        y_second,
        distances[[best_v0], :],
        sources_second,
        return_modelfit=True,
    )
    r2_second = fit_second["rsquared"]

    # Cleanup
    shm_v0.close(); shm_dist.close()
    shm_sources_first.close(); shm_sources_second.close()

    return target_idx, [best_v0, best_sigma, rss_first, r2_first, r2_second]


# ======================================================================
# ------------------- Cross-Validation Main Function -------------------
# ======================================================================
def perform_ccf_analysis_cross(
    OPTIMIZE,
    optimize_threshold,
    distances_along_mesh,
    timeseries_sources_first,
    timeseries_sources_second,
    timeseries_targets_first,
    timeseries_targets_second,
    sigmas,
    n_jobs=30,
):
    """
    Parallel + shared-memory cross-validated CCF analysis.
    UPDATED: Uses float64 (default) to prevent precision loss.
    """
    # --- 1. Precompute candidate fields ---
    connfield_weights = define_connfield_candidates(distances_along_mesh, sigmas)

    # --- 2. Precompute CF timecourses for first half ---
    n_v0, n_sigmas, n_t_first = connfield_weights.shape[0], connfield_weights.shape[1], timeseries_sources_first.shape[1]
    
    # CHANGED: dtype=np.float64 (was float32)
    timeseries_v0_first = np.empty((n_v0, n_sigmas, n_t_first), dtype=np.float64)
    
    for v0_i in range(n_v0):
        for sigma_i in range(n_sigmas):
            timeseries_v0_first[v0_i, sigma_i, :] = create_cf_timecourse(
                timeseries_sources_first, connfield_weights[v0_i, sigma_i, :]
            )

    # --- 3. Create shared memory arrays ---
    # CHANGED: Removed .astype(np.float32). Now uses default (likely float64).
    shm_v0_first, shape_v0_first, dtype_v0_first = _to_shared(timeseries_v0_first)
    shm_dist, shape_dist, dtype_dist = _to_shared(distances_along_mesh)
    shm_sources_first, shape_sources_first, dtype_sources_first = _to_shared(timeseries_sources_first)
    shm_sources_second, shape_sources_second, dtype_sources_second = _to_shared(timeseries_sources_second)

    shm_meta = {
        "names": {
            "timeseries_v0_first": shm_v0_first.name,
            "distances": shm_dist.name,
            "sources_first": shm_sources_first.name,
            "sources_second": shm_sources_second.name,
        },
        "shapes": {
            "timeseries_v0_first": shape_v0_first,
            "distances": shape_dist,
            "sources_first": shape_sources_first,
            "sources_second": shape_sources_second,
        },
        "dtypes": {
            "timeseries_v0_first": dtype_v0_first,
            "distances": dtype_dist,
            "sources_first": dtype_sources_first,
            "sources_second": dtype_sources_second,
        },
    }

    # --- 4. Parallel processing ---
    n_targets = timeseries_targets_first.shape[0]
    best_models = np.zeros((n_targets, 5))

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = [
            executor.submit(
                _fit_target_voxel_cross_shared,
                i,
                timeseries_targets_first[i, :],
                timeseries_targets_second[i, :],
                sigmas,
                OPTIMIZE,
                optimize_threshold,
                shm_meta,
            )
            for i in range(n_targets)
        ]
        for f in as_completed(futures):
            idx, model = f.result()
            best_models[idx, :] = model

    # --- 5. Cleanup shared memory ---
    shm_v0_first.close(); shm_v0_first.unlink()
    shm_dist.close(); shm_dist.unlink()
    shm_sources_first.close(); shm_sources_first.unlink()
    shm_sources_second.close(); shm_sources_second.unlink()

    return {"best_models": best_models, "connfield_weights": connfield_weights}