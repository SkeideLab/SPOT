import nibabel as nib
import numpy as np
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
import pandas as pd

def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
     # Ensure labels_area is a list
    if not isinstance(labels_area, list):
        labels_area = [labels_area]
    
    # Collect indices for all labels in labels_area
    indices_area = np.concatenate([
        np.nonzero(visparc.agg_data() == label)[0]
        for label in labels_area
    ])

    return indices_area

VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")


for hemi in ["L", "R"]:
    # Load the annotation or label data
    visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
    # Load the surface geometry
    print(hemi.lower())
    coords, faces = nib.freesurfer.io.read_geometry(f"/data/p_02915/templates/template_fsaverage/fsaverage/surf/{hemi.lower()}h.sphere")

    v2_vertices = get_indices_roi([3, 4], visparc)
    v3_vertices = get_indices_roi([5, 6], visparc)
    # Step 3: Find border vertices between V2 and V3
    # Create a sparse adjacency matrix for the surface mesh
    # Create adjacency matrix from the faces
    row_indices = np.concatenate([faces[:, 0], faces[:, 1], faces[:, 2]])
    col_indices = np.concatenate([faces[:, 1], faces[:, 2], faces[:, 0]])
    data = np.ones(row_indices.shape[0], dtype=np.float32)  # All edges initially have weight 1

    # The adjacency matrix should use the `row_indices` and `col_indices` arrays for the connections
    adj_matrix = csr_matrix((data, (row_indices, col_indices)), shape=(coords.shape[0], coords.shape[0]))

    # Make the matrix symmetric (undirected graph)
    adj_matrix = adj_matrix + adj_matrix.T

    # Check adjacency of V2 and V3
    border_vertices = []
    for v in v2_vertices:
        neighbors = adj_matrix[v].nonzero()[1]
        if any(n in v3_vertices for n in neighbors):
            border_vertices.append(v)
    border_vertices = np.array(border_vertices)

    # Step 4: Calculate distances from V2 vertices to the border
    border_coords = coords[border_vertices]
    v2_coords = coords[v2_vertices]
    v3_coords = coords[v3_vertices]

    # Calculate the edge weights for the adjacency matrix
    # Edge weights are distances between connected vertices
    row, col = adj_matrix.nonzero()
    distances = np.linalg.norm(coords[row] - coords[col], axis=1)

    # Create a weighted adjacency matrix
    weighted_adj_matrix = csr_matrix((distances, (row, col)), shape=adj_matrix.shape)

    # Step 2: Run Dijkstra's algorithm from all V2 vertices
    # Initialize the distance matrix for shortest paths from V2 vertices
    # Make sure dist_matrix has the correct shape and contains distances for the vertices
    # from both V2 and V3 (i.e., the indices passed to Dijkstra's algorithm)

    # Concatenate the V2 and V3 vertices
    indices = np.concatenate([v2_vertices, v3_vertices])

    # Run Dijkstra's algorithm for the combined indices
    dist_matrix, predecessors = dijkstra(csgraph=weighted_adj_matrix, directed=False, indices=indices, return_predecessors=True)

    # Ensure the size of dist_matrix is correct, containing distances for both V2 and V3 vertices
    assert dist_matrix.shape[0] == len(indices), "Mismatch in dist_matrix size and indices length"

    # Calculate border_distances for the border_vertices
    # Check if border_vertices contains valid indices within the range of dist_matrix
    border_distances = dist_matrix[:, border_vertices]  # Ensure this slicing works correctly

    # If border_vertices contains out-of-bounds indices, filter them out
    border_vertices = [v for v in border_vertices if v < dist_matrix.shape[1]]

    # Now compute the minimum distances for the border
    min_distances = border_distances.min(axis=1)

    # Ensure proper alignment for results
    v2_indices = np.arange(len(v2_vertices))
    v3_indices = np.arange(len(v2_vertices), len(indices))

    # DataFrame for output
    output = pd.DataFrame({
        "Vertex": np.concatenate([v2_vertices, v3_vertices]),
        "Distance_to_Border": np.concatenate([min_distances[:len(v2_vertices)], -min_distances[len(v2_vertices):]]),
    })
    output.to_csv(f"/data/p_02915/dhcp_derivatives_SPOT/{hemi}_polarangle_distance_wang.csv", index=False)    


    # # Continue with the rest of the code
    # for index, row in subject_info.iloc[sub_num:sub_num + 1].iterrows():
    #     if group == "12-16y" or group == "18-21y":
    #         sub_id = subject_info.at[index, "sub_id"]
    #         prefix_model = PREFIX_MODEL.format(
    #             sub=sub_id,
    #             hemi=hemi,
    #         )
    #         input_path = prefix_model         
    #         param_value = "polarangle"
    #         output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub_id}_{hemi}_{param_value}'
    #     else:
    #         sub_id = subject_info.at[index, "sub_id"]
    #         sess_id = subject_info.at[index, "sess_id"]
    #         sub = sub_id.replace('sub-','')
    #         ses = sess_id.replace('ses-','')
    #         prefix_model = PREFIX_MODEL.format(
    #             sub=sub,
    #             ses=ses,
    #             hemi=hemi,
    #         )
    #         input_path = prefix_model
    #         param_value = "polarangle"
    #         output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_{param_value}'

    
    #     # Step 1: Load the polarangle data
    #     # Replace 'polarangle_file' with your actual file path
    #     polarangle_data = nib.load(prefix_model)
    #     polarangle = polarangle_data.darrays[0].data  # polarangle values for all vertices

    #     # Ensure polarangle values align with the surface vertices
    #     assert len(polarangle) == coords.shape[0], "polarangle data doesn't match vertex count."
        

    #     # Correct pairing of distances and polarangles
    #     v2_results = list(zip(min_distances[v2_indices], polarangle[v2_vertices]))
    #     v3_results = list(zip(min_distances[v3_indices], polarangle[v3_vertices]))

    #     # Print or save results for V2 and V3 separately
    #     for i, (dist, ecc) in enumerate(v2_results):
    #         print(f"V2 Vertex {v2_vertices[i]}: Distance to border = {dist:.3f}, polarangle = {ecc:.3f}")
    #     for i, (dist, ecc) in enumerate(v3_results):
    #         print(f"V3 Vertex {v3_vertices[i]}: Distance to border = {dist:.3f}, polarangle = {ecc:.3f}")

    #     # DataFrame for output
    #     output = pd.DataFrame({
    #         "Vertex": np.concatenate([v2_vertices, v3_vertices]),
    #         "Distance_to_Border": np.concatenate([min_distances[v2_indices], min_distances[v3_indices]]),
    #         "polarangle": np.concatenate([polarangle[v2_vertices], polarangle[v3_vertices]])
    #     })
    #     output.to_csv(f"{output_path}_distance.csv", index=False)
    #     print(output_path)
