import subprocess
from pathlib import Path
import pandas as pd
import sys

script_dir = Path(__file__).resolve().parent

project_results_script = str(
    script_dir / "03_retinotopyanalysis" / "project_correlation_to_fsaverage.sh"
)
path_derivatives = Path("/data/p_02915/dhcp_derivatives_SPOT/fetal")
path_anat_data = Path("/data/pt_02880/Package_1225541/fmriresults01/dhcp_anat_pipeline/")  # path_derivatives / "dhcp_anat_pipeline"
path_func_data = Path("/data/pt_02880/Package_1225541/fmriresults01/dhcp_fmri_pipeline/") # path_derivatives / "dhcp_fmri_pipeline"
path_output_data = path_derivatives / "dhcp_surface"
path_templates = Path("/data/p_02915/templates")
path_fsaverage = path_templates / "template_fsaverage/fsaverage"
path_HCPtemplates_standardmeshatlases = (
    "/data/u_yoos_software/HCPpipelines/"
    "global/templates/standard_mesh_atlases"
)
# --------------PATHS TO SOFTWARE--------------------------------------------------
path_newmsm = "/data/u_yoos_software/fsl/bin/msm-env/bin/newmsm"
path_wbcommand = "/bin/wb_command"
path_mirtk = "/afs/cbs.mpg.de/software/mirtk/0.20231123/debian-bullseye-amd64/bin/mirtk"

# Iterate over the index array using a for loop
subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
sub_num = int(sys.argv[1])

for index, row in subject_info.iloc[sub_num:sub_num + 1].iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    sess_id = subject_info.at[index, "sess_id"]
    sub = sub_id.replace('sub-','')
    ses = sess_id.replace('ses-','')
    path_input_file=(f"/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_fmri_pipeline/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_task-rest_desc-preproc_bold.nii.gz")
    path_output_file=(f"/data/p_02915/dhcp_derivatives_SPOT/dhcp_surface/sub-{sub}/ses-{ses}/func/sub-{sub}_ses-{ses}_mcf")
# for sub, ses in sub_ses_todo:

    path_scaninfo = path_anat_data / f"sub-{sub}" / f"sub-{sub}_sessions.tsv"

    cmd_project_results = [
    "FSL",
    "mcflirt",
    "-in",
    path_input_file,
    "-out",
    path_output_file,
    "-plots",
    ]
    subprocess.run(cmd_project_results, check=True)