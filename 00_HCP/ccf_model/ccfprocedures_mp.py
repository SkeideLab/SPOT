import numpy as np
from scipy.optimize import minimize
from .ccflogic import (
    create_cf_fit,
    define_connfield_candidates,
    create_cf_timecourse,
    calc_correlation,
)
import multiprocessing as mp

def parallel_create_cf_timecourse(args):
    """Wrapper function to create cf timecourse for multiprocessing."""
    v0_i, sigma_i, timeseries_sources, connfield_weights = args
    return v0_i, sigma_i, create_cf_timecourse(timeseries_sources, connfield_weights[v0_i, sigma_i, :])

def parallel_gridsearch(args):
    """Wrapper function to perform grid search for multiprocessing."""
    targetvoxel_i, timeseries_targetvoxel, timeseries_v0 = args
    return targetvoxel_i, gridsearch_connfield(timeseries_targetvoxel, timeseries_v0)

def parallel_optimize(args):
    """Wrapper function to perform optimization for multiprocessing."""
    targetvoxel_i, best_model, distances_along_mesh, timeseries_sources, timeseries_target = args
    return targetvoxel_i, optimize_connfield(
        best_model[0], best_model[1], distances_along_mesh, timeseries_sources, timeseries_target
    )

def perform_ccf_analysis(
    OPTIMIZE,
    optimize_threshold,
    distances_along_mesh,
    timeseries_sources,
    timeseries_targets,
    sigmas,
):
    """Runs cortical connective field analysis between source and target vertices using multiprocessing."""

    # Define connective field candidates
    connfield_weights = define_connfield_candidates(distances_along_mesh, sigmas)

    # Prepare multiprocessing for timeseries_v0 computation
    timeseries_v0 = np.empty(
        (
            connfield_weights.shape[0],  # n_nodes_v0
            connfield_weights.shape[1],  # n_sigmas
            timeseries_sources.shape[1],  # n_timepoints
        )
    )

    # Parallel processing for timeseries_v0
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(
            parallel_create_cf_timecourse,
            [
                (v0_i, sigma_i, timeseries_sources, connfield_weights)
                for v0_i in range(connfield_weights.shape[0])
                for sigma_i in range(connfield_weights.shape[1])
            ],
        )

    # Populate timeseries_v0 with results
    for v0_i, sigma_i, cf_timecourse in results:
        timeseries_v0[v0_i, sigma_i, :] = cf_timecourse

    # Parallel processing for grid search
    best_models = np.zeros((timeseries_targets.shape[0], 4))

    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = pool.map(
            parallel_gridsearch,
            [
                (targetvoxel_i, timeseries_targetvoxel, timeseries_v0)
                for targetvoxel_i, timeseries_targetvoxel in enumerate(timeseries_targets)
            ],
        )

    # Populate best_models with results
    for targetvoxel_i, best_model in results:
        best_models[targetvoxel_i, :] = best_model

    # Replace sigma index with sigma value
    best_models[:, 1] = [sigmas[int(sigma_i)] for sigma_i in best_models[:, 1]]

    # Optimization step
    if OPTIMIZE:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = pool.map(
                parallel_optimize,
                [
                    (targetvoxel_i, best_models[targetvoxel_i, :], distances_along_mesh, timeseries_sources, timeseries_targets[targetvoxel_i, :])
                    for targetvoxel_i in range(timeseries_targets.shape[0])
                    # Uncomment the following if you want to optimize only above threshold:
                    # if best_models[targetvoxel_i, 3] > optimize_threshold
                ],
            )

        # Populate best_models with optimized results
        for targetvoxel_i, optimized_model in results:
            best_models[targetvoxel_i, :] = optimized_model

    return {"best_models": best_models, "connfield_weights": connfield_weights}


def gridsearch_connfield(timeseries_targetvoxel, timeseries_v0):
    """Find the best fitting connective field for 1 targetvoxel
        out of a number of pre-computed candidates

    Args:
        timeseries_targetvoxel (numpy.array): n_timepoints. Activity timeseries of target voxel.
        timeseries_v0 (numpy.array): n_centervoxels x n_sigmas x n_timepoints.
            Timeseries of connective field candidates defined by center voxel and spread.

    Returns:
        numpy.array: (v0_i, sigma_i, RSS, rsquared) array of best fitting model. I
            ndex of center voxel, index of sigma, residual sum of squares, r squared.
    """

    # fill with v0_i, sigma_i, RSS, rsquared
    model_rss = [None, None, -np.inf, None]

    # try all possible models
    for v0_i in range(timeseries_v0.shape[0]):
        for sigma_i in range(timeseries_v0.shape[1]):
            modelfit = calc_correlation(
                timeseries_targetvoxel, timeseries_v0[v0_i, sigma_i, :]
            )
            if modelfit > model_rss[2]:
                model_rss = [v0_i, sigma_i, modelfit, modelfit**2]

    return np.array(model_rss)


def optimize_connfield(
    v0_i,
    sigma_init,
    distances,
    timeseries_sources,
    timeseries_target,
):
    """Nonlinear search for best-fitting connective field for 1 target vertex.
    Based on initial best connective field parameters from grid search.
    CF center voxel will not be optimized, only sigma.

    Args:
        v0_i (int): index of model center vertex of previously bestfitting model
        sigma_init (float): spread parameter of previously bestfitting model
        distances (numpy.array): n_vertices_source x n_vertices_source.
            Distances between vertices in graph, sum of shortest edge distances between vertices.
        timeseries_sources (numpy.array): n_vertices_source x n_timepoints.
            Timeseries of source vertices.
        timeseries_target (numpy.array):  n_timepoints.
            Timeseries of target vertex

    Returns:
        tuple: (v0_i, sigma, rss, rsquared) of optimized model
    """
    # minimize rss
    res = minimize(
        create_cf_fit,
        (sigma_init,),
        args=(
            timeseries_target,
            distances[[int(v0_i)], :],
            timeseries_sources,
        ),
        bounds = [(0.1, None)]
    )

    # extract optimization results
    sigma = res.x[0]
    rss = res.fun

    # calculate rsquared with optimized model
    rsquared = create_cf_fit(
        (sigma,),
        timeseries_target,
        distances[[int(v0_i)], :],
        timeseries_sources,
        return_modelfit=True,
    )

    return (v0_i, sigma, rss, rsquared)
