import numpy as np
from nilearn import surface
import pandas as pd

FUNC_PATH = (
    "{root_dir}/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-{simulated}_rss.gii"
)
subject_info = pd.read_csv("/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv")
median_r = pd.DataFrame(columns=["sub_id", "left_real", "right_real", "left_simulated", "right_simulated"])
path_derivatives = "/data/p_02915/dhcp_derivatives_SPOT"
for index, row in subject_info.iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    sess_id = subject_info.at[index, "sess_id"]
    sub = sub_id.replace('sub-','')
    ses = sess_id.replace('ses-','')
    median_r.at[index, "sub_id"] = sub_id
    for hemi in ["L", "R"]:
            if hemi == "L": 
                hemi_down = "left"
            elif hemi == "R":
                hemi_down = "right"
            ids = {
                "sub": sub,
                "ses": ses,
                "hemi": hemi,
                "root_dir": "/data/p_02915/dhcp_derivatives_SPOT/fetal"
                }

            # load retinotopy data
            for datasource in ["real", "simulated"]:
                print(f"Modeling {datasource} data on hemisphere {hemi}.")
                simulated = "simulated" if datasource == "simulated" else "real"
                try:
                    surface_data = surface.load_surf_data(FUNC_PATH.format(**ids, simulated=simulated))
                    num_vertices = surface_data.shape
                    surface_data = surface_data.astype(float)
                    print(surface_data.shape)
                    surface_data = surface_data[surface_data != 0]
                    print(surface_data.shape)
                    median_correlation = np.median(surface_data)
                    median_r.at[index, f"{hemi_down}_{datasource}"] = median_correlation
                except ValueError:
                    print(
                        f"No model results found under {datasource}."
                    )
                    print(" Skipping this...", flush=True)
                    continue

            # get bold data for V1 and V2
            
            
            
#median_r.to_csv('/data/p_02915/SPOT/median_r_fetal_young.csv')