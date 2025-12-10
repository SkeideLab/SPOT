"""
Cacluate tSNR
"""
import pandas as pd
import numpy as np
from nilearn import surface
import nibabel as nib


# Function to convert elements to float, setting non-convertible elements to np.nan
def to_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.nan
    
def flatten(arr):
    return arr.flatten()
  
original_H_value=[]
original_z_score=[]
original_p_value=[]
bootstrap_p_value=[]
bootstrap_ci_lows=[]
bootstrap_ci_highs=[]
index_label = []


for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]: 
    results_list = []    
    if group == "neonates<37":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/Neonates/dhcp_surface/sub-{sub}/ses-{ses}/func/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        )
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_no_drop_v2.csv')
        sub_num = len(subject_info["sub_id"])   
    elif group == "neonates>37":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/Neonates/dhcp_surface/sub-{sub}/ses-{ses}/func/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        )
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_no_drop_v2.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "fetal<29":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/func/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        )
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
        sub_num = len(subject_info["sub_id"])        
    elif group == "fetal>29":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/func/"
            "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        )
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "12-16y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/"
            "{sub}_hemi-{hemi}-Atlas-MSMAll_hp0_clean_bold.func.gii"
        )
        PREFIX_MODEL_2 = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/"
            "{sub}_hemi-{hemi}_mean.func.gii"
        )
        subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "18-21y":
        PREFIX_MODEL = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/"
            "{sub}_hemi-{hemi}-Atlas-MSMAll_hp0_clean_bold.func.gii"
        )
        PREFIX_MODEL_2 = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/"
            "{sub}_hemi-{hemi}_mean.func.gii"
        )
        subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
        sub_num = len(subject_info["sub_id"])
            
    parameters_L=[]
    parameters_R=[]
    # FORMAT PATHS FOR INPUT AND OUTPUT       
    for index, row in subject_info.iterrows():                
        if group == "12-16y" or group == "18-21y":
            for hemi in ["L", "R"]:        
                sub_id = subject_info.at[index, "sub_id"]
                prefix_model = PREFIX_MODEL.format(
                    sub=sub_id,
                    hemi=hemi,
                )
                input_path = prefix_model
                mean_model = PREFIX_MODEL_2.format(
                    sub=sub_id,
                    hemi=hemi,
                )

                if hemi == "L":
                    left = surface.load_surf_data(input_path).astype(np.float64)
                    averaged_brain = surface.load_surf_data(mean_model).astype(np.float64)
                    left = left + averaged_brain[:, np.newaxis]
                    cleaned_data = np.array([[to_float(x) for x in row] for row in left], dtype=np.float64)
                    mean_signal = np.mean(cleaned_data, axis=1)                    
                    std_signal = np.std(cleaned_data, axis=1)
                    std_signal[std_signal == 0] =1e-6
                    tSNR = mean_signal / std_signal
                    print(tSNR)
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub_id}_{hemi}_tSNR.gii'
                    darray = nib.gifti.GiftiDataArray(np.float32(tSNR))
                    params_img = nib.gifti.GiftiImage(darrays=[darray])
                    nib.save(params_img, output_path)
                    parameters_L.append(np.mean(tSNR))

                elif hemi == "R":
                    right = surface.load_surf_data(input_path).astype(np.float64)
                    averaged_brain = surface.load_surf_data(mean_model).astype(np.float64)
                    right = right + averaged_brain[:, np.newaxis]
                    cleaned_data = np.array([[to_float(x) for x in row] for row in right], dtype=np.float64)
                    mean_signal = np.mean(cleaned_data, axis=1)
                    std_signal = np.std(cleaned_data, axis=1)
                    std_signal[std_signal == 0] =1e-6
                    tSNR = mean_signal / std_signal
                    print(tSNR)
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub_id}_{hemi}_tSNR.gii'
                    darray = nib.gifti.GiftiDataArray(np.float32(tSNR))
                    params_img = nib.gifti.GiftiImage(darrays=[darray])
                    nib.save(params_img, output_path)
                    parameters_R.append(np.mean(tSNR))
        
        else:
            for hemi in ["L", "R"]:        
                sub_id = subject_info.at[index, "sub_id"]
                sess_id = subject_info.at[index, "sess_id"]
                sub = sub_id.replace('sub-','')
                ses = sess_id.replace('ses-','')
                prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                input_path = prefix_model

                if hemi == "L":
                    left = surface.load_surf_data(input_path).astype(np.float64)
                    cleaned_data = np.array([[to_float(x) for x in row] for row in left], dtype=np.float64)
                    mean_signal = np.mean(cleaned_data, axis=1)
                    std_signal = np.std(cleaned_data, axis=1)
                    std_signal[std_signal == 0] =1e-6
                    tSNR = mean_signal / std_signal
                    print(tSNR)
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_tSNR.gii'
                    darray = nib.gifti.GiftiDataArray(np.float32(tSNR))
                    params_img = nib.gifti.GiftiImage(darrays=[darray])
                    nib.save(params_img, output_path)
                    parameters_L.append(np.mean(tSNR))

                elif hemi == "R":
                    right = surface.load_surf_data(input_path).astype(np.float64)
                    cleaned_data = np.array([[to_float(x) for x in row] for row in right], dtype=np.float64)
                    mean_signal = np.mean(cleaned_data, axis=1)
                    std_signal = np.std(cleaned_data, axis=1)
                    std_signal[std_signal == 0] =1e-6
                    tSNR = mean_signal / std_signal
                    print(tSNR)
                    output_path=f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_tSNR.gii'
                    darray = nib.gifti.GiftiDataArray(np.float32(tSNR))
                    params_img = nib.gifti.GiftiImage(darrays=[darray])
                    nib.save(params_img, output_path)
                    parameters_R.append(np.mean(tSNR))

    df = {f"{group}_L":parameters_L, f"{group}_R": parameters_R}
results_df = pd.DataFrame(df)
results_df.to_csv(f'/data/p_02915/SPOT/02_averaged_tSNR.csv')