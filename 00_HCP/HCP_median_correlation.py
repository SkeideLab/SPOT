import numpy as np
from nilearn import surface
import pandas as pd

FUNC_PATH = (
    "{root_dir}/ccfmodel/{sub}/"
    "{sub}_hemi-{hemi}_mesh-native_dens-native_desc-{simulated}_rss.gii"
)


def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        np.logical_or(
            visparc.agg_data() == labels_area[0], visparc.agg_data(
            ) == labels_area[1]
        )
    )[0]
    return indices_area


VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/{sub}/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")
LABELS_V2 = (3, 4)
subject_info = pd.read_csv(
    "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv")
median_r = pd.DataFrame(
    columns=["sub_id", "left_real", "right_real", "left_simulated", "right_simulated"])
path_derivatives = "/data/p_02915/dhcp_derivatives_SPOT"
for index, row in subject_info.iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    sub = sub_id.replace('sub-', '')
    median_r.at[index, "sub_id"] = sub_id
    for hemi in ["L", "R"]:
        if hemi == "L":
            hemi_down = "left"
        elif hemi == "R":
            hemi_down = "right"
        ids = {
            "sub": sub_id,
            "hemi": hemi,
            "root_dir": "/data/p_02915/dhcp_derivatives_SPOT/HCP-D"
        }

        # load retinotopy data
        for datasource in ["real", "simulated"]:
            print(f"Modeling {datasource} data on hemisphere {hemi}.")
            simulated = "simulated" if datasource == "simulated" else "real"
            try:
                surface_data = surface.load_surf_data(
                    FUNC_PATH.format(**ids, simulated=simulated))
                num_vertices = surface_data.shape
                surface_data = surface_data.astype(float)
                print(surface_data.shape)
                surface_data = surface_data[surface_data != 0]
                print(surface_data.shape)
                median_correlation = np.median(surface_data)
                median_r.at[index,
                            f"{hemi_down}_{datasource}"] = median_correlation
            except ValueError:
                print(
                    f"No model results found under {datasource}."
                )
                print(" Skipping this...", flush=True)
                continue

        # get bold data for V1 and V2


# median_r.to_csv('/data/p_02915/SPOT/median_r_HCP_old.csv')
