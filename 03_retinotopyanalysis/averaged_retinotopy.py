"""
Average ccf results
"""

import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np

PREFIX_MODEL = (
    "/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel_var/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k"
)
INPUT = (
    "{prefix_model}_label-{param}_desc-{model}_roi-v2th00_metric.gii"
)
#INPUT = (
#    "{prefix_model}_desc-{model}_{param}.gii"
# )
OUTPUT = (
    "/data/p_02915/dhcp_derivatives_SPOT/Neonates/ccfmodel_var/Averaged_older_n_{hemi}_label-{param}_desc-{model}_roi-v2th00_metric.gii"
)
# OUTPUT = (
#    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_younger_fetal_{hemi}_desc-{model}_{param}.gii"
# )
# INPUT = (
#    "{prefix_model}_desc-real_{parm}.gii"
# )
# OUTPUT = (
#    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_{hemi}_desc-real_{parm}.gii"
# )
subject_info = pd.read_csv(
    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_no_drop_v2.csv')
sub_num = len(subject_info["sub_id"])

sum_data = []  # Initialize outside the loop
for param in ["eccentricity", "polarangle"]:
#for param in ["sigma"]:
    # for param in ["rss","sigma"]:
    for hemi in ["L", "R"]:
        # FORMAT PATHS FOR INPUT AND OUTPUT

        for model in ["real"]:
            sum_data = []
            for index, row in subject_info.iterrows():
                sub_id = subject_info.at[index, "sub_id"]
                sess_id = subject_info.at[index, "sess_id"]
                sub = sub_id.replace('sub-', '')
                ses = sess_id.replace('ses-', '')
                prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                output_path = OUTPUT.format(
                    prefix_model=prefix_model, hemi=hemi, model=model, param=param
                )
                input_path = INPUT.format(
                    prefix_model=prefix_model, model=model, param=param
                )
                try:
                    ccf = surface.load_surf_data(
                        input_path).astype(np.float64)                  
                
                    sum_data.append(ccf)

                except ValueError as e:
                    print(f"ValueError: {e}")
                    print(
                        f"No model results found under {INPUT.format(prefix_model=prefix_model, model=model, param = param)}."
                    )
                    print(" Skipping this...", flush=True)
                    continue
            sum_data_array = np.array(sum_data)            
            #sum_data_array = np.where(sum_data_array == 0, np.nan, sum_data_array)
            averaged_data = np.mean(sum_data_array, axis=0)
            print(sum_data_array)
            print(averaged_data)
            # save as gifti
            img_gifti = nib.gifti.GiftiImage(
                darrays=[
                    nib.gifti.GiftiDataArray(
                        np.float32(np.squeeze(averaged_data)))
                ]
            )

            nib.save(img_gifti, output_path)
