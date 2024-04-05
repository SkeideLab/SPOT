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
    script_dir / "03_retinotopyanalysis" / "fill_ccfmodel_with_retinotopydata.py"
)

for sub, ses in sub_ses_todo:

    path_scaninfo = path_anat_data / f"sub-{sub}" / f"sub-{sub}_sessions.tsv"
    age = get_age(ses, path_scaninfo)

    cmd_surfaceprep = [surfaceprep_script, sub, ses, age, str(path_derivatives)]
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
