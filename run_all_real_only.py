import re
import subprocess
from pathlib import Path
import sys
import pandas as pd


# /data/u_kieslinger_software/miniconda3/bin/conda run -n fmritools --live-stream python /data/u_kieslinger_software/code/spot/run_prep_and_ccf.py
def get_age(index, subject_info):
    age = subject_info.at[index, "scan_age"]
    if age > 44:
        age = 44
    return f"{age:.0f}"  # round to whole week


path_derivatives = Path("/data/p_02915/dhcp_derivatives_SPOT")
path_anat_data = Path("/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_anat_pipeline")  # path_derivatives / "dhcp_anat_pipeline"
path_func_data = Path("/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_fmri_pipeline") # path_derivatives / "dhcp_fmri_pipeline"
path_output_data = path_derivatives / "dhcp_surface"
path_templates = Path("/data/p_02915/templates")

# templates
name_volume_template_40wks = "dhcp40wk"
file_volume_template_40wks = (
    path_templates
    / "template_augmentedvolumetricatlas_dhcp/atlas/T2/template-40.nii.gz"
)
name_surface_template = "dhcpSym"
path_surface_template = (
    path_templates
    / "template_corticalsurfaceneonatessym_williams2023_dhcp/dhcpSym_template"
)
path_fsaverage = path_templates / "template_fsaverage/fsaverage"
path_HCPtemplates_standardmeshatlases = (
    "/data/u_yoos_software/HCPpipelines/"
    "global/templates/standard_mesh_atlases"
)
# --------------PATHS TO SOFTWARE--------------------------------------------------
path_newmsm = "/data/u_yoos_software/fsl/bin/msm-env/bin/newmsm"
path_wbcommand = "/bin/wb_command"
path_mirtk = "/afs/cbs.mpg.de/software/mirtk/0.20231123/debian-bullseye-amd64/bin/mirtk"

subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT.csv')
sub_num = int(sys.argv[1])

script_dir = Path(__file__).resolve().parent
surfaceprep_script = str(script_dir / "01_dataprep" / "run_surfaceprep_real_only.sh")
ccfmodel_script = str(script_dir / "02_ccfanalysis" / "run_model_real_only.py")
fillretinotopy_script = str(
    script_dir / "03_retinotopyanalysis" / "analyse_retinotopy.py"
)
project_results_script = str(
    script_dir / "03_retinotopyanalysis" / "project_ccf_to_fsaverage.sh"
)
for index, row in subject_info.iloc[sub_num:sub_num + 1].iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    sess_id = subject_info.at[index, "sess_id"]
    sub = sub_id.replace('sub-','')
    ses = sess_id.replace('ses-','')
# for sub, ses in sub_ses_todo:

    path_scaninfo = path_anat_data / f"sub-{sub}" / f"sub-{sub}_sessions.tsv"
    age = get_age(index, subject_info)
    


    cmd_surfaceprep = [
        str(option)
        for option in [
            surfaceprep_script,
            sub,
            ses,
            age,
            path_derivatives,
            path_anat_data,
            path_func_data,
            path_output_data,
            file_volume_template_40wks,
            name_volume_template_40wks,
            path_surface_template,
            name_surface_template,
            path_HCPtemplates_standardmeshatlases,
            path_fsaverage,
            path_newmsm,
            path_wbcommand,
            path_mirtk,
        ]
    ]
    subprocess.run(cmd_surfaceprep, check=True)

    cmd_ccfmodel = [
        "python",
        ccfmodel_script,
        "-d",
        str(path_derivatives),
        "-sub",
        sub,
        "-ses",
        ses,
    ]
    subprocess.run(cmd_ccfmodel, check=True)
