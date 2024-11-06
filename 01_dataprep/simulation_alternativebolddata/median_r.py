import numpy as np
from nilearn import surface
import pandas as pd
import nibabel as nib

def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
     # Ensure labels_area is a list
    if not isinstance(labels_area, list):
        labels_area = [labels_area]
    
    # Collect indices for all labels in labels_area
    indices_area = np.concatenate([
        np.nonzero(visparc.agg_data() == label)[0]
        for label in labels_area
    ])

    return indices_area
LABELS_V2 = [2, 3]
FUNC_PATH = (
    "{root_dir}/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-{simulated}_r.gii"
)
subject_info = pd.read_csv(
    "/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv")
median_r = pd.DataFrame(
    columns=["sub_id", "left_real", "right_real", "left_simulated", "right_simulated"])
path_derivatives = "/data/p_02915/dhcp_derivatives_SPOT/fetal"
VISPARC_PATH = (
    "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii")
percent= []
for index, row in subject_info.iterrows():

    sub_id = subject_info.at[index, "sub_id"]
    sess_id = subject_info.at[index, "sess_id"]
    sub = sub_id.replace('sub-', '')
    ses = sess_id.replace('ses-', '')
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
            "root_dir": path_derivatives
        }
        visparc = nib.load(
                    VISPARC_PATH.format(**ids))
        indices_v2 = get_indices_roi(LABELS_V2, visparc)

        # load retinotopy data
        for datasource in ["real"]:
            print(f"Modeling {datasource} data on hemisphere {hemi}.")
            simulated = "simulated" if datasource == "simulated" else "real"
            try:
                surface_data = surface.load_surf_data(
                    FUNC_PATH.format(**ids, simulated=simulated))[indices_v2].astype(np.float64)
                num_vertices = surface_data.shape
                surface_data = surface_data.astype(float)
                total=len(surface_data)
                surface_data = surface_data[surface_data > 0.1]
                over=len(surface_data)
                percent.append(over/total*100)
                median_correlation = np.median(surface_data)
                median_r.at[index,
                            f"{hemi_down}_{datasource}"] = median_correlation
            except ValueError:
                print(
                    f"No model results found under {FUNC_PATH.format(**ids, simulated=simulated)}."
                )
                print(" Skipping this...", flush=True)
                continue
print(np.mean(percent))

        # get bold data for V1 and V2


# median_r.to_csv('/data/p_02915/SPOT/median_r_fetal_young.csv')
