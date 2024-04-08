import re
import subprocess
from pathlib import Path

import pandas as pd


# /data/u_kieslinger_software/miniconda3/bin/conda run -n fmritools --live-stream python /data/u_kieslinger_software/code/spot/run_prep_and_ccf.py
def get_age(ses, path_scaninfo):
    scaninfo = pd.read_csv(str(path_scaninfo), sep="\t")
    age = scaninfo.query(f"session_id=={ses}").at[0, "scan_age"]
    return f"{age:.0f}"  # round to whole week


path_derivatives = Path("/data/p_02495/dhcp_derivatives")
path_anat_data = path_derivatives / "dhcp_anat_pipeline"
path_func_data = path_derivatives / "dhcp_fmri_pipeline"
path_output_data = path_derivatives / "dhcp_surface"
path_templates = Path("/data/p_02495/templates")

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
    "/data/u_kieslinger_software/code/HCPpipelines/"
    "global/templates/standard_mesh_atlases"
)
# --------------PATHS TO SOFTWARE--------------------------------------------------
path_newmsm = "/data/u_kieslinger_software/fsldevdir/bin/newmsm"
path_wbcommand = "/bin/wb_command"
path_mirtk = "/afs/cbs.mpg.de/software/mirtk/0.20231123/debian-bullseye-amd64/bin/mirtk"

sub_ses_anat = {
    (
        re.search(r"(?<=sub-)[^_/]*", str(anat_file)).group(),
        re.search(r"(?<=ses-)[^_/]*", str(anat_file)).group(),
    )
    for anat_file in path_anat_data.glob("*/*/*/sub-*ses-*.gii")
}
sub_ses_func = {
    (
        re.search(r"(?<=sub-)[^_/]*", str(anat_file)).group(),
        re.search(r"(?<=ses-)[^_/]*", str(anat_file)).group(),
    )
    for anat_file in path_func_data.glob("*/*/*/sub-*ses-*.nii.gz")
}
sub_ses_todo = sub_ses_anat.intersection(sub_ses_func)

script_dir = Path(__file__).resolve().parent
surfaceprep_script = str(script_dir / "01_dataprep" / "run_surfaceprep.sh")
ccfmodel_script = str(script_dir / "02_ccfanalysis" / "run_model.py")
fillretinotopy_script = str(
    script_dir / "03_retinotopyanalysis" / "analyse_retinotopy.py"
)
project_results_script = str(
    script_dir / "03_retinotopyanalysis" / "project_ccf_to_fsaverage.sh"
)

for sub, ses in sub_ses_todo:

    path_scaninfo = path_anat_data / f"sub-{sub}" / f"sub-{sub}_sessions.tsv"
    age = get_age(ses, path_scaninfo)

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

    cmd_fillretinotopy = [
        "python",
        fillretinotopy_script,
        "-d",
        str(path_derivatives),
        "-sub",
        sub,
        "-ses",
        ses,
    ]
    subprocess.run(cmd_fillretinotopy, check=True)

    cmd_project_results = [
        project_results_script,
        sub,
        ses,
        str(path_anat_data),
        str(path_derivatives),
        str(path_HCPtemplates_standardmeshatlases),
        str(path_fsaverage),
        str(path_wbcommand),
    ]
    subprocess.run(cmd_project_results, check=True)
