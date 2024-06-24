import numpy as np
import statsmodels.api as sm


def calc_CF(distances, sigma):
    """Calculates spread of connective field and gives all cf nodes a weight for this cf.
    Will calculate one CF for each node in the first dimension of 'distances'.

    Args:
        distances (numpy.array): n_nodes_from x n_nodes_to.
            Shortest distance along the mesh between pairs of nodes.
            All nodes_from will have a cf centered on them.
        sigma (float): spread of cf along cortical mesh

    Returns:
        numpy.array: n_nodes_center x n_nodes_all. n_nodes_center=n_nodes_from,
            n_nodes_all=n_nodes_to.
            Cortical fields centered on all 'from' nodes with weights for all nodes.
    """
    return np.exp(-(distances**2 / (2 * sigma**2)))


def define_connfield_candidates(distances, sigmas):
    """Generates candidate connective field spreads for multiple center vertices and spreads.
        Calculates the weight of each vertex in that connective field.

    Args:
        distances (numpy.array): n_vertices x n_vertices. distances between all pairs of vertices
        sigmas (numpy.array): sigma values to use for the candidates.
            Sigma controls the spread of the connective field.

    Returns:
        numpy.array: n_vertices x n_sigmas x n_vertices.
            For each center vertex (dim 1) and each sigma (dim 2), the weights for all other
            vertices in the connective field (dim 3).
    """
    connfield_weights = np.zeros((distances.shape[0], len(sigmas), distances.shape[1]))
    for i, sigma in enumerate(sigmas):
        connfield_weights[:, i, :] = calc_CF(
            distances=distances,
            sigma=sigma,
        )
    return connfield_weights


def create_cf_timecourse(timeseries_sources, connfield):
    """Creates a timecourse of the whole connective field.
        The timeseries of each node in the field is weighed according to the CF and then
        the values of the contributing nodes at each timepoint are summed up.

    Args:
        timeseries_sources (numpy.array): n_vertices x n_timepoints.
            Contains the timeseries of each vertex in the CF
        connfield (numpy.array): n_vertices. Weights of all nodes in this connective field

    Returns:
        numpy.array: n_timepoints,. Linear weighted sum of the node activations per timepoint.
    """
    # n_timepoints x n_nodes
    connective_field_activation = timeseries_sources.T * connfield

    # sum along nodes
    return np.sum(connective_field_activation, axis=1)


def create_cf_fit(
    x, timeseries_targetvoxel, distances, timeseries_sources, return_modelfit=False
):
    """Create connective field weights based on center voxel and sigma and fit this
        model to a targetvoxel timeseries. To be used with optimize function.

    Args:
        x (list): contains only 1 value. sigma (spread of the conective field)
        timeseries_targetvoxel (numpy.array): n_timepoints.
            Activity timeseries of voxel to be modeled as Y.
        distances (numpy.array): 1 x n_voxels_sources. Distances
            from CF center voxel to all other voxels in the source roi.
        timeseries_sources (numpy.array): n_voxels_sources x n_timepoints.
            Timeseries of all voxels in the source roi.
        return_modelfit (bool, optional): Should the modelfit dict ('rss', 'rsquared')
            be returned (if True) or only the rss/residual sum of squares value (if False).
            Defaults to False.

    Returns:
        float/dict : value of model rss or dict consisting of rss and rsquared
    """
    sigma = x[0]
    weights_cf = define_connfield_candidates(distances, [sigma])
    timeseries_cf = create_cf_timecourse(timeseries_sources, np.squeeze(weights_cf))

    if return_modelfit:
        return calc_correlation(timeseries_targetvoxel, np.squeeze(timeseries_cf))
    return 1 - calc_correlation(timeseries_targetvoxel, np.squeeze(timeseries_cf))


def calc_correlation(timeseries_targetvoxel, timeseries_cf):
    r = np.corrcoef(timeseries_targetvoxel, timeseries_cf)[0, 1]
    return r
