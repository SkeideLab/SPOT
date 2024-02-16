# test all the functions of meshgraph module
import numpy as np
from collections import namedtuple
from ..ccfanalysis.meshgraph import find_neighbours, calc_3d_dist, calc_shortestpath
import pytest


@pytest.fixture
def example_mesh():
    # this example validated by hand
    Mesh = namedtuple("mesh", ["coordinates", "faces"])
    return {
        "mesh": Mesh(
            coordinates=np.array(
                [[0, 1, 1], [0, 2, 1], [3, 1, 1], [0, 3, 1], [-10, 3.7, 5]]
            ),
            faces=np.array([[0, 1, 2], [3, 2, 1], [1, 4, 3]]),
        ),
        "neighbours": np.array(
            [
                [0, 1, 1, 0, 0],
                [1, 0, 1, 1, 1],
                [1, 1, 0, 1, 0],
                [0, 1, 1, 0, 1],
                [0, 1, 0, 1, 0],
            ]
        ),
        "neighbour_distances": np.array(
            [
                [0.0, 1.0, 3.0, 0.0, 0],
                [1.0, 0.0, 3.16227766, 1.0, 10.903669],
                [3.0, 3.16227766, 0.0, 3.60555128, 0],
                [0.0, 1.0, 3.60555128, 0.0, 10.7930533],
                [0, 10.903669, 0, 10.7930533, 0],
            ]
        ),
        "shortestpaths": np.array(
            [
                [0.0, 1.0, 3.0, 2.0, 11.90366911],
                [1.0, 0.0, 3.16227766, 1.0, 10.90366911],
                [3.0, 3.16227766, 0.0, 3.60555128, 14.06594677],
                [2.0, 1.0, 3.60555128, 0.0, 10.79305332],
                [11.90366911, 10.90366911, 14.06594677, 10.79305332, 0.0],
            ]
        ),
    }


def test_neighbours_symmetric(example_mesh):

    neighbours = find_neighbours(example_mesh["mesh"])
    assert np.allclose(neighbours, neighbours.T)


def test_neighbours_result(example_mesh):

    expected_result = example_mesh["neighbours"]
    neighbours = find_neighbours(example_mesh["mesh"])
    assert np.allclose(neighbours, expected_result)


def test_calculate_3d_distance_euclidean_validate(example_mesh):
    assert np.allclose(
        calc_3d_dist(example_mesh["mesh"], example_mesh["neighbours"]),
        example_mesh["neighbour_distances"],
    )


def test_shortestpath_calculation(example_mesh):
    mesh = example_mesh["mesh"]
    shortestpaths = calc_shortestpath(mesh, np.arange(mesh.coordinates.shape[0]))
    print(shortestpaths)
    assert np.allclose(shortestpaths, example_mesh["shortestpaths"])
