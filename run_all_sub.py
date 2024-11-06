import subprocess
from pathlib import Path
import pandas as pd
import sys

script_dir = Path(__file__).resolve().parent

project_results_script = str(
    script_dir / "03_retinotopyanalysis" / "project_correlation_to_fsaverage.sh"
)
path_derivatives = Path("/data/p_02915/dhcp_derivatives_SPOT/Neonates")
path_anat_data = Path("/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_anat_pipeline")  # path_derivatives / "dhcp_anat_pipeline"
path_func_data = Path("/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_fmri_pipeline/") # path_derivatives / "dhcp_fmri_pipeline"
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

project_results_script_2 = str(
    script_dir / "03_retinotopyanalysis" / "project_tSNR.sh"
)
# Iterate over the index array using a for loop
subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_v2.csv')
sub_num = int(sys.argv[1])

for index, row in subject_info.iloc[sub_num:sub_num + 1].iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    sess_id = subject_info.at[index, "sess_id"]
    sub = sub_id.replace('sub-','')
    ses = sess_id.replace('ses-','')
    cmd_project_results2 = [
        project_results_script_2,
        sub,
        ses,
        str(path_anat_data),
        str(path_derivatives),
        str(path_HCPtemplates_standardmeshatlases),
        str(path_fsaverage),
        str(path_wbcommand),
    ]
    subprocess.run(cmd_project_results2, check=True)

