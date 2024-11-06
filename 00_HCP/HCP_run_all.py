import re
import subprocess
from pathlib import Path
import sys
import pandas as pd


path_derivatives = Path("/data/p_02915/dhcp_derivatives_SPOT/HCP-D")
# path_derivatives / "dhcp_anat_pipeline"
path_raw_data = Path("/data/pt_02880/HCP_D/fmriresults01")
path_output_data = path_derivatives / "hcp_surface"
path_templates = Path("/data/p_02915/templates")

# templates
path_fsaverage = path_templates / "template_fsaverage/fsaverage"
path_HCPtemplates_standardmeshatlases = (
    "/data/u_yoos_software/HCPpipelines/"
    "global/templates/standard_mesh_atlases"
)
# --------------PATHS TO SOFTWARE--------------------------------------------------
path_newmsm = "/data/u_yoos_software/fsl/bin/msm-env/bin/newmsm"
path_wbcommand = "/bin/wb_command"
path_mirtk = "/afs/cbs.mpg.de/software/mirtk/0.20231123/debian-bullseye-amd64/bin/mirtk"

subject_info = pd.read_csv(
    '/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
sub_num = int(sys.argv[1])

script_dir = Path(__file__).resolve().parent
file_prep_script = str(script_dir / "HCP_file_prep.sh")
resample_template_script = str(script_dir / "HCP_warp_templates_to_native.sh")
ccfmodel_script = str(script_dir / "HCP_ccf_model.py")
fillretinotopy_script = str(
    script_dir / "HCP_retinotopy_analysis.py"
)
project_results_script = str(
    script_dir / "HCP_native_to_fsaverage.sh"
)
project_results_script2 = str(
    script_dir / "HCP_parameter_to_fsaverage.sh"
)
for index, row in subject_info.iloc[sub_num:sub_num + 1].iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    path_anat_data = str(path_raw_data / sub_id /
                         "MNINonLinear" / "fsaverage_LR32k")
    path_func_data = str(path_raw_data / sub_id /
                         "MNINonLinear" / "Results" / "rfMRI_REST")

    ############################### Resample CIFTI to GIFTI#############################
    cmd_fileprep = [
        str(option)
        for option in [
            file_prep_script,
            sub_id,
            path_anat_data,
            path_func_data,
            path_output_data,
            path_wbcommand,
        ]
    ]
    #subprocess.run(cmd_fileprep, check=True)

    ############################### RESAMPLING RETINOTOPY TEMPLATE#############################
    cmd_warp_template = [
        str(resample_template_script),
        str(sub_id),
        str(path_anat_data),
        str(path_output_data),
        str(path_HCPtemplates_standardmeshatlases),
        str(path_fsaverage),
    ]
    #subprocess.run(cmd_warp_template, check=True)

    ######################## SIMULATING ALTERNATIVE DATA###################################
    cmd_simulation_model = [
        "python",
        "/data/p_02915/SPOT/00_HCP/HCP_simulation_model.py",
        "-da",
        str(path_anat_data),
        "-d",
        str(path_derivatives),
        "-sub",
        str(sub_id)
    ]
    #subprocess.run(cmd_simulation_model, check=True)

    cmd_ccfmodel = [
        "python",
        ccfmodel_script,
        "-d",
        str(path_derivatives),
        "-sub",
        sub_id,
        "-th",
        '0'
    ]

    subprocess.run(cmd_ccfmodel, check=True)

    cmd_fillretinotopy = [
        "python",
        fillretinotopy_script,
        "-d",
        str(path_derivatives),
        "-sub",
        str(sub_id),
        "-th",
        '00',
    ]
    subprocess.run(cmd_fillretinotopy, check=True)

    cmd_project_results = [
        project_results_script,
        str(sub_id),
        str(path_anat_data),
        str(path_derivatives),
        str(path_HCPtemplates_standardmeshatlases),
        str(path_fsaverage),
        str(path_wbcommand),
    ]
    subprocess.run(cmd_project_results, check=True)

    cmd_project_results2 = [
        project_results_script2,
        str(sub_id),
        str(path_anat_data),
        str(path_derivatives),
        str(path_HCPtemplates_standardmeshatlases),
        str(path_fsaverage),
        str(path_wbcommand),
    ]
    subprocess.run(cmd_project_results2, check=True)

