import re
import subprocess
from pathlib import Path
import sys
import pandas as pd
import time
from datetime import timedelta

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
    '/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
sub_num = int(sys.argv[1])

script_dir = Path(__file__).resolve().parent
file_prep_script = str(script_dir / "surface_prep.sh")
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
runtime_records = []
for index, row in subject_info.iloc[sub_num:sub_num + 1].iterrows():
    sub_id = str(subject_info.at[index, "sub_id"])
    print(f"\n============================")
    print(f"Processing subject: {sub_id}")
    print(f"============================")

    # Define paths
    path_anat_data = str(path_raw_data / sub_id / "MNINonLinear" / "Native")
    path_func_data = str(path_raw_data / sub_id / "MNINonLinear" / "Results" / "rfMRI_REST")

    # -------------------- SURFACE PREP (optional) --------------------
    cmd_surfaceprep = [
        str(file_prep_script),
        sub_id,
        str(path_derivatives),
        path_anat_data,
        path_func_data,
        str(path_output_data),
        str(path_HCPtemplates_standardmeshatlases),
        str(path_fsaverage),
        str(path_newmsm),
        str(path_wbcommand),
        str(path_mirtk),
    ]
    subprocess.run(cmd_surfaceprep, check=True)

    # -------------------- SIMULATION MODEL (optional) --------------------
    cmd_simulation_model = [
        "/data/u_yoos_software/miniforge3/envs/SPOT/bin/python",
        "/data/p_02915/SPOT/00_HCP/HCP_simulation_model.py",
        "-da", path_anat_data,
        "-d", str(path_derivatives),
        "-sub", sub_id,
    ]
    subprocess.run(cmd_simulation_model, check=True)

    # -------------------- CCF MODEL (TRACK RUNTIME) --------------------
    cmd_ccfmodel = [
        "/data/u_yoos_software/miniforge3/envs/SPOT/bin/python",
        ccfmodel_script,
        "-d", str(path_derivatives),
        "-sub", sub_id,
        "-th", "0"
    ]

    print(f"\n>>> Starting CCF model for {sub_id}...")
    start_time = time.time()

    try:
        subprocess.run(cmd_ccfmodel, check=True)
        success = True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running CCF model for {sub_id}: {e}")
        success = False

    end_time = time.time()
    elapsed = timedelta(seconds=end_time - start_time)
    print(f"✅ Completed CCF model for {sub_id} in {elapsed}")

    runtime_records.append({
        "subject_id": sub_id,
        "success": success,
        "runtime_seconds": round(end_time - start_time, 2),
        "runtime_hms": str(elapsed)
    })

    # -------------------- RETINOTOPY FILL (optional) --------------------
    cmd_fillretinotopy = [
        "/data/u_yoos_software/miniforge3/envs/SPOT/bin/python",
        fillretinotopy_script,
        "-d", str(path_derivatives),
        "-sub", sub_id,
        "-th", "00",
    ]
    subprocess.run(cmd_fillretinotopy, check=True)

    # -------------------- PROJECT RESULTS (optional) --------------------
    cmd_project_results = [
        project_results_script,
        sub_id,
        path_anat_data,
        str(path_derivatives),
        str(path_HCPtemplates_standardmeshatlases),
        str(path_fsaverage),
        str(path_wbcommand),
    ]
    subprocess.run(cmd_project_results, check=True)

    cmd_project_results2 = [
        project_results_script2,
        sub_id,
        path_anat_data,
        str(path_derivatives),
        str(path_HCPtemplates_standardmeshatlases),
        str(path_fsaverage),
        str(path_wbcommand),
    ]
    subprocess.run(cmd_project_results2, check=True)