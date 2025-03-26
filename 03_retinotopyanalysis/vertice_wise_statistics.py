"""
Statistical comparison of ccf results for six groups
"""
import numpy as np
import pandas as pd 
from scipy.stats import kruskal, mannwhitneyu
from multiprocessing import Pool, cpu_count
import math
import nibabel as nib
from nilearn import surface
import pyvista as pv


def mann_whitney_u_to_z(U, n1, n2):
    # Calculate mean (μ_U)
    mu_U = (n1 * n2) / 2

    # Calculate standard deviation (σ_U)
    sigma_U = math.sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)

    # Calculate z-score
    z = (U - mu_U) / sigma_U

    return z


def bootstrap_kruskal_resample(n_resamples, group_sizes, combined_data):
    bootstrap_H_statistics = []
    for _ in range(n_resamples):
        # Resample the combined data
        resampled_indices = np.random.choice(
            len(combined_data), size=len(combined_data), replace=True)
        resampled_data = combined_data[resampled_indices]

        # Split the resampled data back into groups of the original sizes
        resampled_groups = []
        start_idx = 0
        for size in group_sizes:
            resampled_groups.append(resampled_data[start_idx:start_idx + size])
            start_idx += size

        # Compute the Kruskal-Wallis statistic for the resampled groups
        bootstrap_H_statistic, _ = kruskal(*resampled_groups)
        bootstrap_H_statistics.append(bootstrap_H_statistic)
    return bootstrap_H_statistics


def bootstrap_resample_mannwhitney(n_resamples, combined_data, n_groupA, n_groupB):
    bootstrap_mannwhitney_statistics = []
    for _ in range(n_resamples):
        # Resample the combined data
        resample_indices = np.random.choice(
            len(combined_data), size=len(combined_data), replace=True)
        resampled_combined_data = combined_data[resample_indices]

        # Split the resampled data into two groups
        bootstrap_groupA = flatten(resampled_combined_data[:n_groupA])
        bootstrap_groupB = flatten(resampled_combined_data[n_groupA:])

        # Perform Mann-Whitney U test on the bootstrap sample
        bootstrap_mw_statistic, _ = mannwhitneyu(
            bootstrap_groupA, bootstrap_groupB, alternative="two-sided")
        bootstrap_mannwhitney_statistics.append(bootstrap_mw_statistic)
    return bootstrap_mannwhitney_statistics

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

PREFIX_MODEL = (
    "{root_dir}/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k"
)
PREFIX_MODEL_2 = (
    "{root_dir}/ccfmodel/{sub}/"
    "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k"
)
PREFIX_SUB_TEMPLATE = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/"
    "hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii"
)


PATH_v0 = "{prefix_model}_{param}.gii"

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

geo_lst = {'sl': sphere_L, 'sr':sphere_R, 'il':inflated_L, 'ir':inflated_R}
sp_list = ['sl', 'sr']


for param in ["label-eccentricity_desc-real_roi-v2th00_metric", "label-polarangle_desc-real_roi-v2th00_metric"]: #,"label-eccentricity_desc-real_roi-v2th00_metric","label-polarangle_desc-real_roi-v2th00_metric"
    if param =="desc-real_r":
        test_value = "r"
    elif param == "desc-real_sigma":
        test_value = "sigma"
    elif param == "label-eccentricity_desc-real_roi-v2th00_metric":
        test_value = "eccentricity"
    elif param =="label-polarangle_desc-real_roi-v2th00_metric":
        test_value = "polarangle"
    for hemi in ["L", "R"]:
        visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
        indices_v2 = get_indices_roi(LABELS_V2, visparc)
        groups = {}
        for group in ["2nd", "3rd", "preterm", "fullterm", "adolescent", "adult"]:
            parameters = []
            if group == "preterm":
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_v2.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "fullterm":
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_v2.csv')
                sub_num = len(subject_info["sub_id"])           
            elif group == "2nd":
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "3rd":
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "adolescent":
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "adult":
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
                sub_num = len(subject_info["sub_id"])

        # FORMAT PATHS FOR INPUT AND OUTPUT
            for index, row in subject_info.iterrows():
                if group == "adolescent" or group == "adult":
                    sub_id = subject_info.at[index, "sub_id"]                
                    prefix_model = PREFIX_MODEL_2.format(
                    root_dir='/data/p_02915/dhcp_derivatives_SPOT/HCP-D',
                    sub=sub_id,
                    hemi=hemi,
                    )                
                    ccf = surface.load_surf_data(
                        PATH_v0.format(prefix_model=prefix_model, param = param)
                    )
                    ccf_v0 = ccf.astype(np.float64)
                    parameters.append(ccf_v0)

                elif group == "preterm" or group == "fullterm":
                    sub_id = subject_info.at[index, "sub_id"]
                    sess_id = subject_info.at[index, "sess_id"]                
                    sub = sub_id.replace('sub-', '')
                    ses = sess_id.replace('ses-', '')
                    prefix_model = PREFIX_MODEL.format(
                        root_dir="/data/p_02915/dhcp_derivatives_SPOT/Neonates",
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                    )
                    ccf = surface.load_surf_data(
                        PATH_v0.format(prefix_model=prefix_model, param = param)
                    )
                    ccf_v0 = ccf.astype(np.float64)
                    parameters.append(ccf_v0)
                
                elif group == "2nd" or group == "3rd":
                    sub_id = subject_info.at[index, "sub_id"]
                    sess_id = subject_info.at[index, "sess_id"]
                    sub = sub_id.replace('sub-', '')
                    ses = sess_id.replace('ses-', '')
                    prefix_model = PREFIX_MODEL.format(
                        root_dir="/data/p_02915/dhcp_derivatives_SPOT/fetal",
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                    )
                    
                    ccf = surface.load_surf_data(
                        PATH_v0.format(prefix_model=prefix_model, param = param)
                    )
                    ccf_v0 = ccf.astype(np.float64)
                    parameters.append(ccf_v0)
         
            groups[group] = np.array(parameters)
        
        # Generate all pairs of groups
        from itertools import combinations
        group_pairs = list(combinations(groups.keys(), 2))
        for (groupA_name, groupB_name) in group_pairs:
            if hemi == "L":
                sph ="sl"
            elif hemi == "R":
                sph ="sr"

            face_sample = geo_lst[sph][1]  # 3 vertices composing a triangle
            vertice_num = np.full((face_sample.shape[0], 1), 3) 
            face = np.hstack((vertice_num, face_sample)) # add the number of vertices, here is 3
            face = face.flatten() # make a face data array such that [num.vertices, vertex1, vertex2, vertex3, num.vertices, vertex4, ..]

            # Create a pyvista mesh object
            mesh = pv.PolyData()
            mesh = pv.PolyData(geo_lst[sph][0], face)
            visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
            mu_deri = np.zeros((mesh.points.shape[0],1))
            indices_v2 = get_indices_roi(LABELS_V2, visparc)
            groupA = groups[groupA_name]
            groupB = groups[groupB_name]

            for vertex_index in indices_v2:
                groupA_flatt = groupA[:, vertex_index]
                groupB_flatt = groupB[:, vertex_index]
                # Combine data for Kruskal-Wallis test
                combined_data = np.concatenate([groupA, groupB])
                n_groupA = len(groupA)
                n_groupB = len(groupB)

                mw_statistic, mw_p_value = mannwhitneyu(
                    groupA_flatt, groupB_flatt, alternative='two-sided')
                
                if mw_p_value < 0.001:
                    mu_deri[vertex_index] = 1

            output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{groupA_name}_{groupB_name}_{hemi}_{test_value}_vertexwise.gii'
            darray = nib.gifti.GiftiDataArray(np.float32(mu_deri))
            params_img = nib.gifti.GiftiImage(darrays=[darray])
            nib.save(params_img, output_path)
            print(output_path)
