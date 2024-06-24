import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np

PREFIX_MODEL = (
    "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
    "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k"
)
#sub-CC00056XX07_ses-10700_hemi-L_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric
INPUT = (
    "{prefix_model}_label-{parm}_desc-{model}_roi-v2th00_metric.gii"
)
OUTPUT = (
    "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/Averaged_{hemi}_label-{parm}_desc-{model}_roi-v2th00_metric_old.gii"
)
subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
sub_num = len(subject_info["sub_id"])

sum_data = []  # Initialize outside the loop
for param in ["eccentricity", "polarangle"]:          
    for hemi in ["L", "R"]:
        # FORMAT PATHS FOR INPUT AND OUTPUT       
           
        for model in ["real"]:
            sum_data=[]    
            for index, row in subject_info.iterrows(): 
                sub_id = subject_info.at[index, "sub_id"]
                prefix_model = PREFIX_MODEL.format(
                    sub=sub_id,
                    hemi=hemi,
                )
                output_path = OUTPUT.format(
                    prefix_model=prefix_model, hemi=hemi, model=model, parm=param
                )
                input_path =  INPUT.format(
                    prefix_model=prefix_model, model=model, parm=param
                )
                try:
                    ccf_v0 = surface.load_surf_data(input_path).astype(np.float64)
                    sum_data.append(ccf_v0)

                except ValueError:
                    print(
                        f"No model results found under {INPUT.format(prefix_model=prefix_model, model=model, parm = param)}."
                    )
                    print(" Skipping this...", flush=True)
                    continue
            sum_data_array = np.array(sum_data)
            averaged_data = np.mean(sum_data_array, axis=0)
            print(sum_data_array)
            print(averaged_data)
            # save as gifti
            img_gifti = nib.gifti.GiftiImage(
                darrays=[
                    nib.gifti.GiftiDataArray(np.float32(np.squeeze(averaged_data)))
                ]
            )

            nib.save(img_gifti, output_path)


                

                
                

