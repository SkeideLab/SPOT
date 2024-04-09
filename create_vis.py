# create graphs for model visualization
import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csgraph
from scipy.spatial import distance


def find_neighbours(coordinates, faces):
    """Find neighbours in a graph represented as namedtuple with fields faces and coordinates.
    Faces contains triangles of vertex indices. Coordinates contains the x,y,z coordinates in space
    per vertex. Checks which vertices are direct neighbours in a triangle. Undirected.

    Args:
        mesh (nilearn.surface.mesh): surface representation of brain as e.g. read
            from gifti file with nilearn.surface.load_surf_mesh.
            mesh.faces will be used to check which vertices are direct neighbours.

    Returns:
        numpy.array: array of the form n_vertices x n_vertices. Symmetric.
            1 if two vertices are direct neighbours,
            0 otherwise.
    """
    n_vertices = coordinates.shape[0]
    neighbours = np.zeros((n_vertices, n_vertices), dtype=np.int8)

    for face in faces:
        # a face consists of 3 edges and forms a triangle
        # note which vertices are direct neighbours in a triangle
        neighbours[face[0], face[1]] = 1
        neighbours[face[1], face[0]] = 1
        neighbours[face[1], face[2]] = 1
        neighbours[face[2], face[1]] = 1
        neighbours[face[0], face[2]] = 1
        neighbours[face[2], face[0]] = 1

    return neighbours


def calc_2d_dist(coordinates, neighbours):
    """Calculate euclidean distance in 3-dimensional space between
    all nodes that are direct neigbours in a mesh.

    Args:
        mesh (nilearn.surface.mesh): brain surface representation with fields coordinates and faces.
            mesh.coordinates will be used as position in space for vertices.
        neighbours (numpy.array): n_vertices x n_vertices.
            1 if two vertices are direct neighbours,
            0 otherwise.

    Returns:
        numpy.array: n_vertices x n_vertices. Symmetric. Undirected.
            Euclidean distance between vertices that are direct neighbours,
            0 otherwise.
    """
    neighbour_distances = np.zeros(neighbours.shape)
    neighbour_index_dim1, neighbour_index_dim2 = np.nonzero(neighbours)
    for i, j in zip(neighbour_index_dim1, neighbour_index_dim2):
        dist = distance.euclidean(coordinates[i], coordinates[j])
        neighbour_distances[i, j] = dist
        neighbour_distances[j, i] = dist

    return neighbour_distances


def calc_shortestpath(coordinates, faces, indices_of_interest):
    """Calculate the shortest path between all pairs of vertices in a graph mesh
    where at least one member is given in indices_of_interest.
    Distance will be a sum of all edge euclidean distances
    in the shortest path (found with dijkstra algorithm) in the graph between vertex a and vertex b.

    Args:
        mesh (nilearn.surface.mesh): brain surface representation as mesh.
            named tuple with fields coordinates and faces.
        indices_of_interest (numpy.array): indices of vertices that are of interest.
            Shortest paths will only be calculated starting from these vertices.

    Returns:
        numpy.array: n_vertices_interest x n_vertices_interest.
            Euclidean distances along surface mesh for all vertices of interest
    """
    neighbours = find_neighbours(coordinates, faces)
    neighbour_distances = calc_2d_dist(coordinates, neighbours)
    gr = csgraph.csgraph_from_dense(neighbour_distances)

    # compute shortest paths only to and from indices of interest
    shortest_paths_interest = csgraph.dijkstra(
        gr, directed=False, indices=indices_of_interest
    )

    return shortest_paths_interest[:, indices_of_interest]


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


def draw_gaussian_on_triangle_grid():
    # 1. gaussian on triangle grid
    plt.style.use("ggplot")

    n_points_per_dim = 20
    mesh = np.meshgrid(
        np.arange(0, n_points_per_dim, 1), np.arange(0, n_points_per_dim, 1)
    )
    x, y = [np.ravel(coords) for coords in mesh]

    t = tri.Triangulation(x, y)
    center_i = 250  # np.random.randint(len(x))

    distances = calc_shortestpath(
        coordinates=np.stack([x, y], axis=1),
        faces=t.triangles,
        indices_of_interest=np.arange(x.shape[0]),
    )

    cf = calc_CF(distances, sigma=5)
    plt.figure()
    tcp = plt.tripcolor(t, cf[center_i, :], edgecolors="w", linewidth=1)
    plt.axis("off")
    plt.show()


def make_two_random_vertex_timeseries():
    plt.style.use("fivethirtyeight")

    n_timepoints = 100
    ts_1 = np.random.normal(size=n_timepoints)
    ts_2 = np.random.normal(size=n_timepoints) / 2

    plt.figure()
    plt.ylim(-2.5, 2.5)
    plt.grid(False)
    plt.axis("off")
    plt.plot(np.arange(0, n_timepoints, 1), ts_1, color="#FDE725FF")
    plt.show(block=True)

    plt.figure()
    plt.ylim(-2.5, 2.5)
    plt.grid(False)
    plt.plot(np.arange(0, n_timepoints, 1), ts_2, color="#481567FF")
    plt.axis("off")
    plt.show()

    plt.figure()
    plt.ylim(-2.5, 2.5)
    plt.grid(False)
    plt.plot(np.arange(0, n_timepoints, 1), (ts_1 + ts_2) / 2, color="#95D840FF")
    plt.axis("off")
    plt.show()


def main():
    make_two_random_vertex_timeseries()
    x = 1


if __name__ == "__main__":
    main()
