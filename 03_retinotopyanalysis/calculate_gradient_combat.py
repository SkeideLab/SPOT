"""
Caculated gradients of eccentricity and polarangle
"""
import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
import pyvista as pv
from neuromaps.datasets import fetch_fsaverage


def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        np.logical_or(
            visparc.agg_data() == labels_area[0], visparc.agg_data() == labels_area[1]
        )
    )[0]
    return indices_area


def flatten(arr):
    return arr.flatten()


VISPARC_PATH = ("/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii")
LABELS_V2 = (2, 3)
 # Read the inflated surface
frsurf = '/data/p_02915/templates/template_fsaverage/fsaverage/surf/'
file_inf_L = frsurf + 'lh.inflated'
file_sph_L = frsurf + 'lh.sphere'
file_inf_R = frsurf + 'rh.inflated'
file_sph_R = frsurf + 'rh.sphere'

# nib.freesurfer.io.read_geometry(file) returns to a tuple of two numpy arrays: vertices and faces
inflated_L = nib.freesurfer.io.read_geometry(file_inf_L)
sphere_L = nib.freesurfer.io.read_geometry(file_sph_L) 
inflated_R = nib.freesurfer.io.read_geometry(file_inf_R)
sphere_R = nib.freesurfer.io.read_geometry(file_sph_R) 

# Order of vertices are preserved in the inflated and sphere surfaces.
# RAS based coordinate system is used in the inflated and sphere surfaces.
# The RAS coordinate of a vertex differs from the inflated and sphere surfaces. 
# Caution: For left hemisphere, increasing R heads to the medial part, 
#          and for right hemisphere, increasing R heads to the lateral part.
# Map the preferred numerosity map to the fsaverage sphere surface for our study,
# since the sphere surface has a uniform distribution (regular grid) of vertices.


# Calculate the gradient of scalar field.
# Load geometry data fsaverage

geo_lst = {'sl': sphere_L, 'sr':sphere_R, 'il':inflated_L, 'ir':inflated_R}
sp_list = ['sl', 'sr']

# Prepare dictionary to save the calculated gradients of each cluster for each subject
grad_dict = {}

# Prepare an array containing 
prog_total = {}
df = pd.DataFrame()
groups = ["2nd", "3rd", "preterm",
          "fullterm", 'adolescent', 'adult']
surfaces = fetch_fsaverage(density="164k")
lh, rh = surfaces["sphere"]
for param in ["eccentricity", "polarangle"]: #"_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric"
    print(param)
    for sph in sp_list:
        if sph == "sl":
            hemi ="L"
            template = surface.load_surf_data(lh)
        elif sph == "sr":
            hemi ="R"
            template = surface.load_surf_data(rh)
        file_path = f"/data/p_02915/SPOT/combat_hemi-{hemi}_area-V2_V3_{param}.csv"
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Parse CSV data manually
        data = []
        for line in lines:
            # Convert each line to a list of floats
            row = [float(val) for val in line.strip().split(',')]
            data.append(row)

        # Convert the list of lists to a NumPy array
        combat_harmonized = np.array(data)
        covars = pd.read_csv(f"/data/p_02915/SPOT/covars_hemi-{hemi}.csv")

        # Group data into categories
        groups_sigma = {               
            '2nd': combat_harmonized[:, covars[covars["group"] == "2nd"].index],
            '3rd': combat_harmonized[:, covars[covars["group"] == "3rd"].index],
            'preterm': combat_harmonized[:, covars[covars["group"] == "preterm"].index],
            'fullterm': combat_harmonized[:, covars[covars["group"] == "fullterm"].index],
            'adolescent': combat_harmonized[:, covars[covars["group"] == "adolescent"].index],
            'adult': combat_harmonized[:, covars[covars["group"] == "adult"].index]
        } 
        for i, group in enumerate(groups):
            sub_num = groups_sigma[group].shape[1]
            print(sub_num)
        
            face_sample = geo_lst[sph][1]  # 3 vertices composing a triangle
            vertice_num = np.full((face_sample.shape[0], 1), 3) 
            face = np.hstack((vertice_num, face_sample)) # add the number of vertices, here is 3
            face = face.flatten() # make a face data array such that [num.vertices, vertex1, vertex2, vertex3, num.vertices, vertex4, ..]

            # Create a pyvista mesh object
            mesh = pv.PolyData()
            mesh = pv.PolyData(geo_lst[sph][0], face)
            visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
            indices_v2 = get_indices_roi(LABELS_V2, visparc)
            indices_v2c = np.nonzero(visparc.agg_data() == 3)[0]
            indices_v2d = np.nonzero(visparc.agg_data() == 4)[0]
            group_name=f"{group}_{hemi}"
            parameters=[]
            ccf = np.zeros(template[0].shape[0])
    # FORMAT PATHS FOR INPUT AND OUTPUT       
            for index in range(0, sub_num):
                if group == "adult":
                    subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
                    sub_id = subject_info.at[index, "sub_id"]
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub_id}_{hemi}_{param}_gradient_combat'
                elif group == "adolescent":
                    subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
                    sub_id = subject_info.at[index, "sub_id"]
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub_id}_{hemi}_{param}_gradient_combat'
                elif group == "fullterm":
                    subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_v2.csv')
                    sub_id = subject_info.at[index, "sub_id"]
                    sess_id = subject_info.at[index, "sess_id"]
                    sub = sub_id.replace('sub-','')
                    ses = sess_id.replace('ses-','')
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_{param}_gradient_combat.gii'
                elif group == "preterm":
                    subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_v2.csv')
                    sub_id = subject_info.at[index, "sub_id"]
                    sess_id = subject_info.at[index, "sess_id"]
                    sub = sub_id.replace('sub-','')
                    ses = sess_id.replace('ses-','')
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_{param}_gradient_combat.gii'
                elif group == "3rd":
                    subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
                    sub_id = subject_info.at[index, "sub_id"]
                    sess_id = subject_info.at[index, "sess_id"]
                    sub = sub_id.replace('sub-','')
                    ses = sess_id.replace('ses-','')
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_{param}_gradient_combat.gii'
                elif group == "2nd":
                    subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
                    sub_id = subject_info.at[index, "sub_id"]
                    sess_id = subject_info.at[index, "sess_id"]
                    sub = sub_id.replace('sub-','')
                    ses = sess_id.replace('ses-','')
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_{param}_gradient_combat.gii'
            
            
            # assgine a point data set of preffered numbers, mu                
                ccf[indices_v2] = groups_sigma[group][:,index]
                mu = ccf
                for i in range(mu.shape[0]):
                    if mu[i]==0:
                        mu[i]=np.nan
                    else:
                        pass
                # Add scalar point set to the mesh
                mesh.point_data['mu'] = mu
                # Compute the gradient of the scalar point set
                mu_deri = np.zeros((mesh.points.shape[0],3))
                # calculate the gradient at a vertex - no parametric space of simplex and interpolation based on the shape of PolyData
                # - but compute it with respect to the direct neighbor vertices (connected with an edge)
                for vertex_index in range(mesh.points.shape[0]):
                    if np.isnan(mesh['mu'][vertex_index]):
                        mu_deri[vertex_index,:] = np.nan
                    else:
                        # Coordinates of the i-th vertex
                        v0 = mesh.points[vertex_index]
                        # Indices of neighboring vertices (connected vertices)
                        neighbors_indices = mesh.point_neighbors(vertex_index)
                        # Filter out neighbors with NaN scalar values
                        valid_neighbors_indices = [i for i in neighbors_indices if not np.isnan(mesh['mu'][i])]
                        # Coordinates of neighboring vertices
                        neighbors = mesh.points[valid_neighbors_indices]
                        # Scalar values of neighboring vertices
                        neighbor_scalars = mesh['mu'][valid_neighbors_indices]
                        # Difference vectors from the i-th vertex to its neighbors
                        diff_vectors = neighbors - v0
                        # Scalar differences
                        scalar_diffs = neighbor_scalars - mesh['mu'][vertex_index]
                        # Weighted sum of differences
                        weighted_sum = np.sum(scalar_diffs[:, np.newaxis] * diff_vectors, axis=0)
                        # Sum of squared lengths of difference vectors
                        squared_lengths = np.sum(np.linalg.norm(diff_vectors, axis=1) ** 2)
                        # Gradient calculation
                        gradient = weighted_sum/squared_lengths
                        #print(f"Gradient at the vertex {vertex_index}:", gradient)
                        mu_deri[vertex_index,:] = gradient               
                
                darray = nib.gifti.GiftiDataArray(np.float32(mu_deri))
                params_img = nib.gifti.GiftiImage(darrays=[darray])
                nib.save(params_img, output_path)
                print(output_path)