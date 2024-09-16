import subprocess
from pathlib import Path
import pandas as pd
import sys


script_dir = Path(__file__).resolve().parent

project_results_script = str(
    script_dir / "HCP_parameter_to_fsaverage.sh"
)
path_derivatives = Path("/data/p_02915/dhcp_derivatives_SPOT/HCP-D")
path_raw_data = Path("/data/pt_02880/HCP_D/fmriresults01")  # path_derivatives / "dhcp_anat_pipeline"
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

subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
#sub_num = int(sys.argv[1])

for index, row in subject_info.iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    sub = sub_id.replace('sub-','')
# for sub, ses in sub_ses_todo:

    sub_id = subject_info.at[index, "sub_id"]
    path_anat_data = str(path_raw_data / sub_id / "MNINonLinear" / "fsaverage_LR32k")

    cmd_project_results = [
    project_results_script,
    sub,
    str(path_anat_data),
    str(path_derivatives),
    str(path_HCPtemplates_standardmeshatlases),
    str(path_fsaverage),
    str(path_wbcommand),
]
    subprocess.run(cmd_project_results, check=True)

