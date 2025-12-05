import nibabel as nib
import numpy as np
import sys
import pandas as pd
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

def calc_msd(param_area1, param_area2):
    # mean squared difference between two different model results in an area
    n1 = param_area1.size
    n2 = param_area2.size
    assert (
        n1 == n2
    ), "areas don't have the same number of vertices and cannot be compared!"

    return np.sum((param_area1 - param_area2) ** 2) / n1


groups = ["neonates_young","neonates_old", "fetal_young", "fetal_old", "adolescent", "adult"]
for group in groups:
    print(group)
    if group == "neonates_young":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel{cross}/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
        )
        PREFIX_SUB_TEMPLATE = (
        "/data/p_02915/dhcp_derivatives_SPOT/Neonates/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
        "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native")
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_v2.csv')
        sub_num = len(subject_info["sub_id"])   
    elif group == "neonates_old":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel{cross}/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
        )
        PREFIX_SUB_TEMPLATE = (
        "/data/p_02915/dhcp_derivatives_SPOT/Neonates/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
        "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native")
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_v2.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal_young":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel{cross}/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
        )
        PREFIX_SUB_TEMPLATE = (
        "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
        "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native")
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
        sub_num = len(subject_info["sub_id"])        
    elif group == "fetal_old":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel{cross}/sub-{sub}/ses-{ses}/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
        )
        PREFIX_SUB_TEMPLATE = (
        "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
        "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native")
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "adolescent":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel{cross}/{sub}/"
            "{sub}_hemi-{hemi}_mesh-native_dens-native"
        )
        PREFIX_SUB_TEMPLATE = (
        "/data/p_02915/SPOT/00_HCP/template/"
        "hemi-{hemi}_mesh-native_dens-native")
        subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "adult":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel{cross}/{sub}/"
            "{sub}_hemi-{hemi}_mesh-native_dens-native"
        )
        PREFIX_SUB_TEMPLATE = (
        "/data/p_02915/SPOT/00_HCP/template/"
        "hemi-{hemi}_mesh-native_dens-native"
        )
    
        subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
        sub_num = len(subject_info["sub_id"])


    PATH_ECCENTRICITY = "{prefix_sub_template}_desc-eccentretinotbenson2014_seg.shape.gii"
    PATH_ANGLE = "{prefix_sub_template}_desc-angleretinotbenson2014_seg.shape.gii"
    PATH_VISPARC = (
        "{prefix_sub_template}_desc-retinotbenson2014_label-visarea_dparc.label.gii"
    )
    OUTPUT_ECC = (
        "{prefix_model}_label-eccentricity_desc-{model}_roi-v2th00_metric.gii"
    )
    OUTPUT_PANGLE = (
        "{prefix_model}_label-polarangle_desc-{model}_roi-v2th00_metric.gii"
    )
    LABELS_V2 = [2]
    LABELS_V3 = [3]

    columns = ["Subject_ID",
            "L_eccentricity_real_benson",
            "R_eccentricity_real_benson",
            "L_polarangle_real_benson",
            "R_polarangle_real_benson",
            "L_eccentricity_real_simulated",
            "R_eccentricity_real_simulated",
            "L_polarangle_real_simulated",
            "R_polarangle_real_simulated",
            "L_eccentricity_cross1_benson",
            "R_eccentricity_cross1_benson",
            "L_polarangle_cross1_benson",
            "R_polarangle_cross1_benson",
            "L_eccentricity_cross2_benson",
            "R_eccentricity_cross2_benson",
            "L_polarangle_cross2_benson",
            "R_polarangle_cross2_benson",
            "L_eccentricity_real_cross1",
            "R_eccentricity_real_cross1",
            "L_polarangle_real_cross1",
            "R_polarangle_real_cross1",
            "L_eccentricity_real_cross2",
            "R_eccentricity_real_cross2",
            "L_polarangle_real_cross2",
            "R_polarangle_real_cross2",
            ]
    msd_value = pd.DataFrame(columns=columns)
    for index, row in subject_info.iterrows():
        retinotopy_dict = {}
        in_out = {
            "eccentricity": OUTPUT_ECC,
            "polarangle": OUTPUT_PANGLE,
        }

        for hemi in ["L", "R"]:
            hemi_results = {}  # Results for the current hemisphere
            if group == "adolescent" or group == "adult":
                sub_id = subject_info.at[index, "sub_id"]
                sub = sub_id.replace('sub-', '')
                prefix_sub_template = PREFIX_SUB_TEMPLATE.format(
                    hemi=hemi,
                )
            else:                
                sub_id = subject_info.at[index, "sub_id"]
                sub = sub_id.replace('sub-', '')
                sess_id = subject_info.at[index, "sess_id"]
                ses = sess_id.replace('ses-', '')
                # FORMAT PATHS FOR INPUT AND OUTPUT
                prefix_sub_template = PREFIX_SUB_TEMPLATE.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )

            # LOAD DATA
            visparc = nib.load(
                PATH_VISPARC.format(prefix_sub_template=prefix_sub_template)
            )
            indices_v2 = get_indices_roi(LABELS_V2, visparc)
            indices_v3 = get_indices_roi(LABELS_V3, visparc)

            retinotopy_dict[(hemi, "benson", "eccentricity")] = surface.load_surf_data(
                PATH_ECCENTRICITY.format(prefix_sub_template=prefix_sub_template)
            )

            retinotopy_dict[(hemi, "benson", "polarangle")] = surface.load_surf_data(
                PATH_ANGLE.format(prefix_sub_template=prefix_sub_template)
            )

            for model in ["cross1", "cross2"]:
                if model == "cross1":
                    cross = 1
                elif model == "cross2":
                    cross =2
                if group == "adolescent" or group == "adult":
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub,
                        hemi=hemi,
                        cross=cross,                    )
                else:    
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                        cross=cross,
                    )
                for param in ["eccentricity", "polarangle"]:

                    retinotopy_out = surface.load_surf_data(in_out[param].format(
                        prefix_model=prefix_model,
                        model="real",
                    ))

                    # keep data for model comparisons later
                    retinotopy_dict[(hemi, model, param)] = retinotopy_out #ret_wholeb

            for model in ["real", "simulated"]:
                if group == "adolescent" or group == "adult":
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub,
                        hemi=hemi,
                        cross="",
                    )
                else:    
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                        cross="",
                    )

                for param in ["eccentricity", "polarangle"]:
                    retinotopy_out = surface.load_surf_data(in_out[param].format(
                        prefix_model=prefix_model,
                        model=model,
                    ))
                    # keep data for model comparisons later
                    retinotopy_dict[(hemi, model, param)] = retinotopy_out #ret_wholeb

            # CALCULATE DIFFERENCES BETWEEN MODELS PER HEMI
            # calculate mean square differences between models for eccentricity and angle
            try:
                for model_1, model_2 in [("real", "benson"), ("real", "simulated"), ("cross1","benson"), ("cross2","benson"), ("real", "cross1"), ("real", "cross2")]:
                    for param in ["eccentricity", "polarangle"]:
                        msd = calc_msd(
                            retinotopy_dict[(hemi, model_1, param)][indices_v3],
                            retinotopy_dict[(hemi, model_2, param)][indices_v3],
                        )
                        # Store the result for the current param
                        msd_value.at[index,
                                    f"{hemi}_{param}_{model_1}_{model_2}"] = msd

            except KeyError:
                print(
                    f"Data missing for hemisphere {hemi}, not calculating differences for this..."
                )

    msd_value.to_csv(f'/data/p_02915/SPOT/MSD_{group}_V3.csv')
