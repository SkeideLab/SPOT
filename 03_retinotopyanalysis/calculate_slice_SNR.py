"""
Calculate slice SNR
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

# Lists to store results
parameters_L = []
parameters_R = []

for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]: #, "neonates<37", "neonates>37", "fetal<29", "fetal>29"]
    results_list = []    
    if group == "neonates<37":
        PREFIX_MODEL = "/data/p_02915/dhcp_derivatives_SPOT/Neonates/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_no_drop_v2.csv')
    elif group == "neonates>37":
        PREFIX_MODEL = "/data/p_02915/dhcp_derivatives_SPOT/Neonates/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_no_drop_v2.csv')
    elif group == "fetal<29":
        PREFIX_MODEL = "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
    elif group == "fetal>29":
        PREFIX_MODEL = "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_bold.func.gii"
        subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
    elif group == "12-16y":
        PREFIX_MODEL = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/{sub}_hemi-{hemi}-Atlas-MSMAll_hp0_clean_bold.func.gii"
        subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
        PREFIX_MODEL_2 = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/"
            "{sub}_hemi-{hemi}_mean.func.gii"
        )
    elif group == "18-21y":
        PREFIX_MODEL = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/{sub}_hemi-{hemi}-Atlas-MSMAll_hp0_clean_bold.func.gii"
        subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
        PREFIX_MODEL_2 = (
            "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/"
            "{sub}_hemi-{hemi}_mean.func.gii"
        )

    # Calculate SNR for each subject
    for index, row in subject_info.iterrows():                
        for hemi in ["L", "R"]:        
            sub_id = subject_info.at[index, "sub_id"]
            sess_id = subject_info.at[index, "sess_id"] if "sess_id" in subject_info.columns else None
            
            if sess_id:
                sub = sub_id.replace('sub-', '')
                ses = sess_id.replace('ses-', '')
                input_path = PREFIX_MODEL.format(sub=sub, ses=ses, hemi=hemi)
            else:
                input_path = PREFIX_MODEL.format(sub=sub_id, hemi=hemi)

            # Load image data
            image_data = surface.load_surf_data(input_path).astype(np.float64)

            if group in ["12-16y", "18-21y"]:
                mean_model = PREFIX_MODEL_2.format(
                    sub=sub_id,
                    hemi=hemi,
                )
                averaged_brain = surface.load_surf_data(mean_model).astype(np.float64)
                image_data = image_data + averaged_brain[:, np.newaxis]
            # Convert to float and clean data
            cleaned_data = np.array([[to_float(x) for x in row] for row in image_data], dtype=np.float64)

            # Compute mean signal intensity in the region of interest
            mean_signal = np.mean(cleaned_data, axis=1)

            # Estimate noise from background (assuming lowest 10% of intensity values represent noise)
            noise_region = np.percentile(cleaned_data, 10, axis=1)
            std_noise = np.std(noise_region)

            # Avoid division by zero
            std_noise = std_noise if std_noise > 1e-6 else 1e-6

            # Compute SNR
            slice_SNR = mean_signal / std_noise
            print(slice_SNR)

            # Save results
            # Determine output file name based on the group type
            if group in ["12-16y", "18-21y"]:
                output_path = f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub_id}_{hemi}_sliceSNR.gii'
            else:
                output_path = f'/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_sliceSNR.gii'

            # Save the slice SNR values
            darray = nib.gifti.GiftiDataArray(np.float32(slice_SNR))
            params_img = nib.gifti.GiftiImage(darrays=[darray])
            nib.save(params_img, output_path)

            darray = nib.gifti.GiftiDataArray(np.float32(slice_SNR))
            params_img = nib.gifti.GiftiImage(darrays=[darray])
            nib.save(params_img, output_path)

            if hemi == "L":
                parameters_L.append(np.mean(slice_SNR))
            elif hemi == "R":
                parameters_R.append(np.mean(slice_SNR))

    # Save to CSV
    df = {f"{group}_L": parameters_L, f"{group}_R": parameters_R}
    results_df = pd.DataFrame(df)
    results_df.to_csv('/data/p_02915/SPOT/Result/02_averaged_slice_SNR.csv')
