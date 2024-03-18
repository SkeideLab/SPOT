import numpy as np
from scipy.optimize import minimize
from ccflogic import (
    create_cf_fit,
    define_connfield_candidates,
    create_cf_timecourse,
    calc_correlation,
)


def perform_ccf_analysis(
    OPTIMIZE,
    optimize_threshold,
    distances_along_mesh,
    timeseries_sources,
    timeseries_targets,
    sigmas,
):
    """Runs cortical connective field analysis between source and target vertices
        (e.g. V1 and V2) on a surface mesh. Will run a grid search over potential
        cortial connective fields for each target vertex and search non-linearly to
        optimize solution if model is explaining some variance. Cortical connective
        fields are defined by their center node v0 and their spread sigma.

    Args:
        OPTIMIZE (bool): should model fits be optimized or only gridsearched
        optimize_threshold (float): 0 < optimize_threshold < 1.
            R squared thresholds for model optimization
        distances_along_mesh (numpy.array): n_vertices_source x n_vertices_source.
            Distances along mesh between vertices
        timeseries_sources (numpy.array): n_vertices_source x n_timepoints.
            Activity timeseries per vertex in source
        timeseries_targets (numpy.array): n_vertices_target x n_timepoints.
            Activity timeseries per vertex in target
        sigmas (numpy.array): n_sigmas. Spreads of cortical fields to try in grid search.

    Returns:
        dict: keys connfield_weights and best_models.
    """
    # make connective field for all nodes and all sigmas
    connfield_weights = define_connfield_candidates(distances_along_mesh, sigmas)

    # empty cf timeseries
    timeseries_v0 = np.empty(
        (
            connfield_weights.shape[0],  # n_nodes_v0 x
            connfield_weights.shape[1],  # n_sigmas x
            timeseries_sources.shape[1],  # n_timepoints
        )
    )
    for v0_i in range(connfield_weights.shape[0]):
        for sigma_i in range(connfield_weights.shape[1]):
            timeseries_v0[v0_i, sigma_i, :] = create_cf_timecourse(
                timeseries_sources,
                connfield_weights[v0_i, sigma_i, :],
            )

    # fill with v0_i, sigma, RSS, rsquared
    best_models = np.zeros((timeseries_targets.shape[0], 4))
    for targetvoxel_i, timeseries_targetvoxel in enumerate(timeseries_targets):

        # compute model fits for previously computed candidate CCFs and targetvoxel
        best_models[targetvoxel_i, :] = gridsearch_connfield(
            timeseries_targetvoxel, timeseries_v0
        )

    # replace sigma index with sigma value
    best_models[:, 1] = [sigmas[int(sigma_i)] for sigma_i in best_models[:, 1]]

    # ---------------------------------------------COMPUTE CONNECTIVE FIELDS-------------------
    if OPTIMIZE:
        for targetvoxel_i, rsquared in enumerate(best_models[:, 3]):
            # if model found some signal (rsquared threshold), then optimize
            if rsquared > optimize_threshold:
                best_models[targetvoxel_i, :] = optimize_connfield(
                    best_models[targetvoxel_i, 0],
                    best_models[targetvoxel_i, 1],
                    distances_along_mesh,
                    timeseries_sources,
                    timeseries_targets[targetvoxel_i, :],
                )

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
