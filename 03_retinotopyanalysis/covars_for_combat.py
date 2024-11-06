import numpy as np
import os
import pandas as pd
import nibabel as nib
from nilearn import surface

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


def flatten(arr):
    return arr.flatten()

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


PATH_v0 = "{prefix_model}_desc-{model}_v0i.gii"

LABELS_V1 = [1]
LABELS_V2 = [2]
LABELS_V3 = [3]


for hemi in ["L", "R"]:
    subject_data_list = []
    tsnr_v2_list = []
    tsnr_v1_list = []
    tsnr_v3_list = []
    tsnr_v2_v1_list = []
    tsnr_roi_list = []
    visparc = nib.load(PREFIX_SUB_TEMPLATE.format(
        hemi=hemi,
        ))
    indices_v1 = get_indices_roi(LABELS_V1, visparc)
    indices_v2 = get_indices_roi(LABELS_V2, visparc)
    indices_v3 = get_indices_roi(LABELS_V3, visparc)
    
    for group in ["2nd", "3rd", "preterm", "fullterm", "adolescent", "adult"]:
        if group == "preterm":
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_v2.csv')
            sub_num = len(subject_info["sub_id"])
            subject_info["sex"] = subject_info['sex'].map({"male": 0, "female": 1})
        elif group == "fullterm":
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_v2.csv')
            sub_num = len(subject_info["sub_id"])            
            subject_info["sex"] = subject_info['sex'].map({"male": 0, "female": 1})
        elif group == "2nd":
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
            sub_num = len(subject_info["sub_id"])
            subject_info["sex"] = subject_info['sex'].map({"M": 0, "F": 1})
        elif group == "3rd":
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
            sub_num = len(subject_info["sub_id"])
            subject_info["sex"] = subject_info['sex'].map({"M": 0, "F": 1})
        elif group == "adolescent":
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
            sub_num = len(subject_info["sub_id"])
            subject_info["sex"] = subject_info['sex'].map({"M": 0, "F": 1})
        elif group == "adult":
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
            sub_num = len(subject_info["sub_id"])
            subject_info["sex"] = subject_info['sex'].map({"M": 0, "F": 1})

        
    # FORMAT PATHS FOR INPUT AND OUTPUT
        for index, row in subject_info.iterrows():
            if group == "adolescent" or group == "adult":
                sub_id = subject_info.at[index, "sub_id"]                
                prefix_model = PREFIX_MODEL_2.format(
                root_dir='/data/p_02915/dhcp_derivatives_SPOT/HCP-D',
                sub=sub_id,
                hemi=hemi,
                )                
                ccf_v0 = surface.load_surf_data(
                    PATH_v0.format(prefix_model=prefix_model, model="real")
                ).astype(int, copy=False)
                tSNR = surface.load_surf_data(f"/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub_id}_{hemi}_tSNR_fsaverage.gii")
                sex = subject_info.at[index, "sex"]
                age = subject_info.at[index, "scan_age"] * 52 / 12 # convert month to week
                site = 1

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
                ccf_v0 = surface.load_surf_data(
                    PATH_v0.format(prefix_model=prefix_model, model="real")
                ).astype(int, copy=False)
                tSNR = surface.load_surf_data(f"/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_tSNR_fsaverage.gii")
                
                sex = subject_info.at[index, "sex"]
                age = subject_info.at[index, "scan_age"]
                site = 2

            
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
                
                ccf_v0 = surface.load_surf_data(
                    PATH_v0.format(prefix_model=prefix_model, model="real")
                ).astype(int, copy=False)
                tSNR = surface.load_surf_data(f"/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_tSNR_fsaverage.gii")
                sex = subject_info.at[index, "sex"] 
                age = subject_info.at[index, "scan_age"]
                site = 2

            centers_v2 = ccf_v0[indices_v2]
            # restrict to V1
            tsnr_v1 = tSNR[indices_v1]
            tsnr_v2 = tSNR[indices_v2]
            tsnr_v3 = tSNR[indices_v3]
            tsnr_v2_v1 = np.full_like(tsnr_v2, np.mean(tsnr_v1))
            #tsnr_roi = np.hstack([tsnr_v1, tsnr_v2])            

            subject_data = [group, site, age, sex, np.mean(tsnr_v1), np.mean(tsnr_v2)]
            subject_data_list.append(subject_data)
            #tsnr_v2_v1_list.append(tsnr_v2_v1)
            tsnr_v2_list.append(tsnr_v2)
            #tsnr_roi_list.append(tsnr_roi)
            tsnr_v1_list.append(tsnr_v1)
            tsnr_v3_list.append(tsnr_v3)

    # Optionally, convert lists to DataFrame (if necessary)
    subject_data_list = np.vstack(subject_data_list)
    # Convert any numpy array to pandas DataFrame
    #arr_tsnr_v2_v1 = np.vstack(tsnr_v2_v1_list)
    arr_tsnr_v2 = np.vstack(tsnr_v2_list)
    #arr_tsnr_roi = np.vstack(tsnr_roi_list)
    arr_tsnr_v3 = np.vstack(tsnr_v3_list)
    arr_tsnr_v1 = np.vstack(tsnr_v1_list)


    #combined_subject_data = pd.DataFrame(subject_data_list, columns=["group", "site", "age", "sex", "tsnr1", "tsnr2"])
    #combined_subject_data.to_csv(f"/data/p_02915/SPOT/covars_hemi-{hemi}.csv")
    #df_tsnr_v2_v1 = pd.DataFrame(arr_tsnr_v2_v1)
    #df_tsnr_v2_v1.to_csv(f"/data/p_02915/SPOT/covars_hemi-{hemi}_tsnr_v2_v1.csv")
    df_tsnr_v2 = pd.DataFrame(arr_tsnr_v2)
    df_tsnr_v2.to_csv(f"/data/p_02915/SPOT/covars_hemi-{hemi}_tsnr_v2.csv")
    #df_tsnr_roi = pd.DataFrame(arr_tsnr_roi)
    #df_tsnr_roi.to_csv(f"/data/p_02915/SPOT/covars_hemi-{hemi}_tsnr_roi.csv")
    df_tsnr_v3 = pd.DataFrame(arr_tsnr_v3)
    df_tsnr_v3.to_csv(f"/data/p_02915/SPOT/covars_hemi-{hemi}_tsnr_v3.csv")
    df_tsnr_v1 = pd.DataFrame(arr_tsnr_v1)
    df_tsnr_v1.to_csv(f"/data/p_02915/SPOT/covars_hemi-{hemi}_tsnr_v1.csv")
    

