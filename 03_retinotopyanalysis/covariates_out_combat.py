import os
import pandas as pd
import numpy as np
import nibabel as nib
from nilearn import surface
from neuroCombat import neuroCombat


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



PREFIX_MODEL = (
    "{root_dir}/ccfmodel_var/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k"
)
PREFIX_MODEL_2 = (
    "{root_dir}/ccfmodel_var/{sub}/"
    "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k"
)
PREFIX_SUB_TEMPLATE = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/"
    "hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii"
)


PATH_v0 = "{prefix_model}_{param}.gii"


for param in ["desc-real_r"]: #  
    if param =="desc-real_r":
        test_value = "r"
    elif param == "desc-real_sigma":
        test_value = "sigma"
    elif param == "label-eccentricity_desc-real_roi-v2th00_metric":
        test_value = "eccentricity"
    elif param =="label-polarangle_desc-real_roi-v2th00_metric":
        test_value = "polarangle"
    for area in ["V2_V3"]:
        if area == "V2":
            LABELS_V2 = [2]
        elif area == "V3":
            LABELS_V2 = [3]
        elif area == "V2_V3":
            LABELS_V2 = [2, 3]

        for hemi in ["L", "R"]:
            print(hemi)
            parameters = []
            combined_subject_data = []
            combat_harmonized=[]
            sub_total = 0
            tsnr_v2_v1_list = []
            visparc = nib.load(PREFIX_SUB_TEMPLATE.format(
                hemi=hemi,
                ))
            indices_v2 = get_indices_roi(LABELS_V2, visparc)
            
            for group in ["2nd", "3rd", "preterm", "fullterm", "adolescent", "adult"]:
                if group == "preterm":
                    subject_info = pd.read_csv(
                        '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_no_drop_v2.csv')
                    sub_num = len(subject_info["sub_id"])
                elif group == "fullterm":
                    subject_info = pd.read_csv(
                        '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_no_drop_v2.csv')
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

                sub_total = sub_total+sub_num
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
                        ccf_v0 = ccf[indices_v2].astype(np.float64)
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
                        ccf_v0 = ccf[indices_v2].astype(np.float64)
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
                        ccf_v0 = ccf[indices_v2].astype(np.float64)
                        parameters.append(ccf_v0)
                print(group)
                print(np.nanmean(np.vstack(parameters)>0))

            processed_parameters = []

            for par in parameters:
                # If the element is empty, replace it with the default value
                if par.size == 0:
                    processed_parameters.append(0)
                else:
                    # Replace None or np.nan values with a default value (e.g., 0)
                    cleaned_param = np.nan_to_num(par, nan=0.0)  # Replaces NaN with 0.0
                    cleaned_param = np.where(cleaned_param == None, 0, cleaned_param)  # Replaces None with 0
                    processed_parameters.append(cleaned_param)

            combined_subject_data = np.vstack(processed_parameters)
            if test_value == "r":
                group_avg_r = np.mean(combined_subject_data, axis=0)  # average for each vertex
                r_min = np.min(group_avg_r)
                r_max = np.max(group_avg_r)
                print(f"Group-averaged r-values for hemi-{hemi}, area-{area}:")
                print(f"  → r_min = {r_min:.6f}, r_max = {r_max:.6f}")

                # Clip the data to avoid r values exactly equal to -1 or 1 (which would cause division by zero)
                combined_subject_data_clipped = np.clip(combined_subject_data, -0.9999, 0.9999)

                # Apply Fisher Z-transform: Z = 0.5 * np.log((1 + r) / (1 - r))
                fisher_transformed_data = 0.5 * np.log((1 + combined_subject_data_clipped) / (1 - combined_subject_data_clipped))

                # Save the Fisher-transformed data
                save_fisher_data = pd.DataFrame(fisher_transformed_data.T)
                #save_fisher_data.to_csv(f"/data/p_02915/SPOT/Result/raw_hemi-{hemi}_area-{area}_{test_value}.csv", index=False, header=False)
            else:
                fisher_transformed_data  = combined_subject_data
                save_fisher_data = pd.DataFrame(combined_subject_data.T)
                save_fisher_data.to_csv(f"/data/p_02915/SPOT/Result/raw_hemi-{hemi}_area-{area}_{test_value}.csv", index=False, header=False)
            #epsilon = 1e-6  # Small constant to avoid zero variance
            #combined_subject_data = combined_subject_data + epsilon
            #save_raw_data = pd.DataFrame(combined_subject_data.T)
            #save_raw_data.to_csv(f"/data/p_02915/SPOT/raw_hemi-{hemi}_area-{area}_{test_value}_raw.csv", index=False, header=False)
           
            print(combined_subject_data.shape)
            # Create a covariates DataFrame
            covars = pd.read_csv(f"/data/p_02915/SPOT/Result/covars_hemi-{hemi}.csv")
            
            # Check unique values in covariates
            print(covars[['site', 'age', 'sex', 'tsnr1', 'tsnr2']].nunique())

            # Ensure all necessary columns are of appropriate type
            covars['site'] = covars['site'].astype('category')
            covars['sex'] = covars['sex'].astype('category')

            # Step 4: Apply ComBat harmonization to the correlation data
            combat_harmonized = neuroCombat(dat=fisher_transformed_data.T,          # Transpose the data
                                            covars=covars[['site', 'age', 'sex', 'tsnr1', 'tsnr2']],               # Covariates DataFrame
                                            batch_col='site',            # Adjust for site
                                            continuous_cols=['age'],     # Preserve 'age'
                                            categorical_cols=['sex'])["data"] # Preserve 'gender'
            save_combat_data = pd.DataFrame(combat_harmonized)
            #save_combat_data.to_csv(f"/data/p_02915/SPOT/Result/combat_hemi-{hemi}_area-{area}_{test_value}.csv", index=False, header=False)
        

