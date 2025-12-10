import nibabel as nib
import numpy as np
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
import pandas as pd

def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest."""
    if not isinstance(labels_area, list):
        labels_area = [labels_area]  # Ensure it's a list

    indices_area = np.concatenate([
        np.nonzero(visparc.agg_data() == label)[0] for label in labels_area
    ])
    return indices_area

VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/"
    "hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii"
)
 
for hemi in ["L", "R"]:
    visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
    coords, faces = nib.freesurfer.io.read_geometry(
        f"/data/p_02915/templates/template_fsaverage/fsaverage/surf/{hemi.lower()}h.sphere"
    )

    # Get vertex indices for dorsal and ventral regions
    v1_dorsal_vertices = get_indices_roi(2, visparc)
    v1_ventral_vertices = get_indices_roi(1, visparc)
    v2_dorsal_vertices = get_indices_roi(4, visparc)
    v2_ventral_vertices = get_indices_roi(3, visparc)
    v3_dorsal_vertices = get_indices_roi(6, visparc)
    v3_ventral_vertices = get_indices_roi(5, visparc)
    out_dorsal_vertices = get_indices_roi([15, 16, 17], visparc)
    out_ventral_vertices = get_indices_roi([8, 9], visparc)

    # Create adjacency matrix from faces
    row_indices = np.concatenate([faces[:, 0], faces[:, 1], faces[:, 2]])
    col_indices = np.concatenate([faces[:, 1], faces[:, 2], faces[:, 0]])
    data = np.ones(len(row_indices), dtype=np.float32)

    adj_matrix = csr_matrix((data, (row_indices, col_indices)), shape=(coords.shape[0], coords.shape[0]))
    adj_matrix = adj_matrix + adj_matrix.T  # Make it symmetric

    def find_border_vertices(region_vertices, adjacent_regions):
        """Find vertices at the border between a region and its neighbors."""
        border_vertices = []
        for v in region_vertices:
            neighbors = adj_matrix[v].nonzero()[1]
            if any(n in np.concatenate(adjacent_regions) for n in neighbors):
                border_vertices.append(v)
        return np.array(border_vertices)

    # Find border vertices separately for dorsal and ventral
    border_vertices_dorsal_v2 = find_border_vertices(v2_dorsal_vertices, [v3_dorsal_vertices, v1_dorsal_vertices])
    border_vertices_dorsal_v3 = find_border_vertices(v3_dorsal_vertices, [out_dorsal_vertices, v2_dorsal_vertices])
    border_vertices_ventral_v2 = find_border_vertices(v2_ventral_vertices, [v3_ventral_vertices, v1_ventral_vertices])
    border_vertices_ventral_v3 = find_border_vertices(v3_ventral_vertices, [out_ventral_vertices, v2_ventral_vertices])

    # Compute adjacency matrix edge weights (Euclidean distances)
    row, col = adj_matrix.nonzero()
    distances = np.linalg.norm(coords[row] - coords[col], axis=1)
    weighted_adj_matrix = csr_matrix((distances, (row, col)), shape=adj_matrix.shape)

    # Run Dijkstra's algorithm separately for dorsal and ventral regions
    v2_vertices_dorsal = v2_dorsal_vertices
    v3_vertices_dorsal = v3_dorsal_vertices
    v2_vertices_ventral = v2_ventral_vertices
    v3_vertices_ventral = v3_ventral_vertices

    dist_matrix_dorsal_v2, _ = dijkstra(csgraph=weighted_adj_matrix, directed=False, indices=v2_vertices_dorsal, return_predecessors=True)
    dist_matrix_dorsal_v3, _ = dijkstra(csgraph=weighted_adj_matrix, directed=False, indices=v3_vertices_dorsal, return_predecessors=True)
    dist_matrix_ventral_v2, _ = dijkstra(csgraph=weighted_adj_matrix, directed=False, indices=v2_vertices_ventral, return_predecessors=True)
    dist_matrix_ventral_v3, _ = dijkstra(csgraph=weighted_adj_matrix, directed=False, indices=v3_vertices_ventral, return_predecessors=True)

    # Filter border vertices within bounds
    border_vertices_dorsal_v2 = border_vertices_dorsal_v2[border_vertices_dorsal_v2 < dist_matrix_dorsal_v2.shape[1]]
    border_vertices_dorsal_v3 = border_vertices_dorsal_v3[border_vertices_dorsal_v3 < dist_matrix_dorsal_v3.shape[1]]
    border_vertices_ventral_v2 = border_vertices_ventral_v2[border_vertices_ventral_v2 < dist_matrix_ventral_v2.shape[1]]
    border_vertices_ventral_v3 = border_vertices_ventral_v3[border_vertices_ventral_v3 < dist_matrix_ventral_v3.shape[1]]

    # Compute distances to borders separately
    border_distances_dorsal_v2 = dist_matrix_dorsal_v2[:, border_vertices_dorsal_v2].min(axis=1)
    border_distances_dorsal_v3 = dist_matrix_dorsal_v3[:, border_vertices_dorsal_v3].min(axis=1)
    border_distances_ventral_v2 = dist_matrix_ventral_v2[:, border_vertices_ventral_v2].min(axis=1)
    border_distances_ventral_v3 = dist_matrix_ventral_v3[:, border_vertices_ventral_v3].min(axis=1)

    # DataFrame for output
    output = pd.DataFrame({
        "Vertex": np.concatenate([v2_vertices_dorsal, v3_vertices_dorsal, v2_vertices_ventral, v3_vertices_ventral]),
        "Distance_to_Border": np.concatenate([border_distances_dorsal_v2, border_distances_dorsal_v3, border_distances_ventral_v2, border_distances_ventral_v3]),
    })
    output.to_csv(f"/data/p_02915/dhcp_derivatives_SPOT/{hemi}_eccentricity_distance_wang.csv", index=False)    
    output_2 = pd.DataFrame({"Vertex": np.concatenate([border_vertices_dorsal_v2, border_vertices_ventral_v2, border_vertices_dorsal_v3, border_vertices_ventral_v3])})
    output_2.to_csv(f"/data/p_02915/dhcp_derivatives_SPOT/{hemi}_eccentricity_border_wang.csv", index=False) 