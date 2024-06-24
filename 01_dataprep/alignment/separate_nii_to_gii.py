"""
Convert Cifti file to gifti (for prenatal dataset)
"""

import pandas as pd
import subprocess

path_wbcommand = "/bin/wb_command"
subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')

for index, row in subject_info.iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    sess_id = subject_info.at[index, "sess_id"]
    sub = sub_id.replace('sub-','')
    ses = sess_id.replace('ses-','')
    for model in ["sulc", "curv"]:
        Input_path = (f"/data/pt_02880/Package_1225541/fmriresults01/dhcp_anat_pipeline/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_{model}.dscalar.nii")
        output_left = (f"/data/pt_02880/Package_1225541/fmriresults01/dhcp_anat_pipeline/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-left_{model}.shape.gii")
        output_right = (f"/data/pt_02880/Package_1225541/fmriresults01/dhcp_anat_pipeline/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-right_{model}.shape.gii")
          
   
        cmd_project_results = [
                path_wbcommand,
                "-cifti-separate",
                Input_path,
                "COLUMN", "-metric", "CORTEX_LEFT",
                output_left,
                "-metric", "CORTEX_RIGHT",
                output_right,
                ]
        subprocess.run(cmd_project_results, check=True)
