"""
Cacluate tSNR
"""
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path


def flatten(arr):
    return arr.flatten()
 
path_HCPtemplates_standardmeshatlases = (
    "/data/u_yoos_software/HCPpipelines/"
    "global/templates/standard_mesh_atlases"
)
path_fsaverage = ("/data/p_02915/templates/template_fsaverage/fsaverage")
path_wbcommand = "/bin/wb_command"

for group in ["12-16y", "18-21y"]: 
    results_list = []    
    if group == "12-16y":
        project_results_script = str(
            "/data/p_02915/SPOT/00_HCP/HCP_concatenate.sh"
        )
        path_raw_data = Path("/data/pt_02880/HCP_D/fmriresults01")
        path_derivatives = Path("/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/")
        subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
        sub_num = len(subject_info["sub_id"])
    elif group == "18-21y":
        project_results_script = str(
            "/data/p_02915/SPOT/00_HCP/HCP_concatenate.sh"
        )
        path_raw_data = Path("/data/pt_02880/HCP_D/fmriresults01")
        path_derivatives = Path("/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/")
        subject_info = pd.read_csv('/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
        sub_num = len(subject_info["sub_id"])

    if group == "12-16y":
        # FORMAT PATHS FOR INPUT AND OUTPUT       
        for index in range(57,60):          
            sub_id = subject_info.at[index, "sub_id"]
            path_anat_data = str(path_raw_data / sub_id /
                        "MNINonLinear" / "fsaverage_LR32k")
            path_func_data = str(path_raw_data / sub_id /
                        "MNINonLinear" / "Results")
            cmd_project_results = [
                project_results_script,
                str(sub_id),
                str(path_anat_data),
                str(path_func_data),
                str(path_derivatives),
                str(path_wbcommand),
            ]
            subprocess.run(cmd_project_results, check=True)
    elif group == "18-21y":
        for index, row in subject_info.iterrows():   
            sub_id = subject_info.at[index, "sub_id"]
            path_anat_data = str(path_raw_data / sub_id /
                         "MNINonLinear" / "fsaverage_LR32k")
            path_func_data = str(path_raw_data / sub_id /
                         "MNINonLinear" / "Results")
            cmd_project_results = [
                project_results_script,
                str(sub_id),
                str(path_anat_data),
                str(path_func_data),
                str(path_derivatives),
                str(path_wbcommand),
            ]
            subprocess.run(cmd_project_results, check=True)
