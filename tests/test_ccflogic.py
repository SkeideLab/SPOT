# test all the functions of meshgraph module
import numpy as np
from collections import namedtuple
from ..ccfanalysis.ccflogic import (
    calc_CF,
    define_connfield_candidates,
    create_cf_timecourse,
)
from ..ccfanalysis.modelPCF import make_percent_signal_change
import pytest


@pytest.fixture
def example_distances():
    distances = np.array(
        [
            [
                0.0,
                0.80319059,
                2.91988581,
                2.18483418,
                2.47272342,
                1.51983225,
                2.16783452,
                2.2252509,
                1.7705344,
                1.20916122,
            ],
            [
                0.80319059,
                0.0,
                3.05637002,
                2.28789002,
                2.18919998,
                0.71664166,
                1.36464393,
                1.42206031,
                1.48701096,
                1.3194527,
            ],
            [
                2.91988581,
                3.05637002,
                0.0,
                0.76848,
                1.40686929,
                2.33972836,
                2.25074518,
                2.76781607,
                1.56935906,
                2.09526038,
            ],
            [
                2.18483418,
                2.28789002,
                0.76848,
                0.0,
                0.65075725,
                1.57124835,
                1.48226517,
                1.99933606,
                0.80087906,
                1.36020875,
            ],
            [
                2.47272342,
                2.18919998,
                1.40686929,
                0.65075725,
                0.0,
                1.47255832,
                1.38357514,
                1.90064603,
                0.70218903,
                1.2635622,
            ],
            [
                1.51983225,
                0.71664166,
                2.33972836,
                1.57124835,
                1.47255832,
                0.0,
                0.64800227,
                0.70541865,
                0.77036929,
                0.60281104,
            ],
            [
                2.16783452,
                1.36464393,
                2.25074518,
                1.48226517,
                1.38357514,
                0.64800227,
                0.0,
                0.51707089,
                0.68138611,
                1.24275929,
            ],
            [
                2.2252509,
                1.42206031,
                2.76781607,
                1.99933606,
                1.90064603,
                0.70541865,
                0.51707089,
                0.0,
                1.198457,
                1.30822968,
            ],
            [
                1.7705344,
                1.48701096,
                1.56935906,
                0.80087906,
                0.70218903,
                0.77036929,
                0.68138611,
                1.198457,
                0.0,
                0.56137317,
            ],
            [
                1.20916122,
                1.3194527,
                2.09526038,
                1.36020875,
                1.2635622,
                0.60281104,
                1.24275929,
                1.30822968,
                0.56137317,
                0.0,
            ],
        ]
    )
    return distances


def test_create_connfield_candidate_weights():
    # validated by visualization(?)
    # TODO validate with example
    define_connfield_candidates


def test_calculate_connective_field_batch(example_distances):
    # validated by visualization

    assert np.allclose(
        calc_CF(example_distances[0, :], 5), calc_CF(example_distances, 5)[0, :]
    )


def test_cf_diagonal_highest(example_distances):
    cf = calc_CF(example_distances, 3)
    assert np.diagonal(cf).all() == 1
    assert np.max(cf) == 1


def test_cf_timecourse_creation():
    v1 = np.array([0, 0, 1, 2, 1, 0, 0])
    v2 = np.array([0, 1, 2, 2, 0, 0, 0])
    timeseries_sources = np.stack([v1, v2])
    timeseries_sources = make_percent_signal_change(timeseries_sources)
    connfield = np.array([1, 0])
    assert np.allclose(
        create_cf_timecourse(timeseries_sources, connfield), timeseries_sources[0, :]
    )
    connfield = np.array([0, 1])
    assert np.allclose(
        create_cf_timecourse(timeseries_sources, connfield), timeseries_sources[1, :]
    )
    connfield = np.array([1, 1])
    assert np.allclose(
        create_cf_timecourse(timeseries_sources, connfield),
        np.sum(timeseries_sources, axis=0),
    )
    connfield = np.array([0, 0])
    assert create_cf_timecourse(timeseries_sources, connfield).all() == 0
