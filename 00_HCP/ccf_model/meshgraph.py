import numpy as np
from scipy.sparse import csgraph
from scipy.spatial import distance


def find_neighbours(mesh):
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
    n_vertices = mesh.coordinates.shape[0]
    neighbours = np.zeros((n_vertices, n_vertices), dtype=np.int8)

    for face in mesh.faces:
        # a face consists of 3 edges and forms a triangle
        # note which vertices are direct neighbours in a triangle
        neighbours[face[0], face[1]] = 1
        neighbours[face[1], face[0]] = 1
        neighbours[face[1], face[2]] = 1
        neighbours[face[2], face[1]] = 1
        neighbours[face[0], face[2]] = 1
        neighbours[face[2], face[0]] = 1

    return neighbours


def calc_3d_dist(mesh, neighbours):
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
        dist = distance.euclidean(mesh.coordinates[i], mesh.coordinates[j])
        neighbour_distances[i, j] = dist
        neighbour_distances[j, i] = dist

    return neighbour_distances

 
def calc_shortestpath(mesh, indices_of_interest):
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
    neighbours = find_neighbours(mesh)
    neighbour_distances = calc_3d_dist(mesh, neighbours)
    gr = csgraph.csgraph_from_dense(neighbour_distances)

    # compute shortest paths only to and from indices of interest
    shortest_paths_interest = csgraph.dijkstra(
        gr, directed=False, indices=indices_of_interest
    )

    return shortest_paths_interest[:, indices_of_interest]
