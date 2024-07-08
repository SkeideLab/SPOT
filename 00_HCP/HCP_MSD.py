"""
Project V1 retinotopy template values to V2 based on cortico-connective field results. 
Print the mean square differences between 
ccf results from real vs simulated data and between 
ccf results from real data vs benson template data.
"""

from argparse import ArgumentParser
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import surface

# FIXED INPUT AND OUTPUT STRUCTURE
PREFIX_MODEL = (
    "{root_dir}/ccfmodel/{sub}/"
    "{sub}_hemi-{hemi}_mesh-native_dens-native"
)
PREFIX_SUB_TEMPLATE = (
    "/data/p_02915/SPOT/00_HCP/template/"
    "hemi-{hemi}_mesh-native_dens-native"
)

PATH_v0 = "{prefix_model}_desc-{model}_v0i.gii"
PATH_r = "{prefix_model}_desc-{model}_r.gii"
OUTPUT_ECC = (
    "{prefix_model}_label-eccentricity_desc-{model}_roi-v2th{threshold}_metric.gii"
)
OUTPUT_PANGLE = (
    "{prefix_model}_label-polarangle_desc-{model}_roi-v2th{threshold}_metric.gii"
)

PATH_ECCENTRICITY = "{prefix_sub_template}_desc-eccentretinotbenson2014_seg.shape.gii"
PATH_ANGLE = "{prefix_sub_template}_desc-angleretinotbenson2014_seg.shape.gii"
PATH_VISPARC = (
    "{prefix_sub_template}_desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
)
LABELS_V1 = (1, 2)
LABELS_V2 = (3, 4)


def parse_args():
    """Parses arguments from the command line."""

    parser = ArgumentParser()
    parser.add_argument(
        "-d",
        "--derivatives_directory",
        required=True,
        help="Superdirectory for top-level dhcp derivatives datasets ",
    )
    parser.add_argument(
        "-sub",
        required=True,
        help="Subject being processed",
    )
    parser.add_argument(
        "-th",
        "--threshold",
        required=False,
        default=0.1,
        help=(
            "Correlation coefficient (r) threshold from model results. "
            "Vertices that pass the threshold will be filled with retinotopy data."
        ),
    )

    args = parser.parse_args()
    return args


def get_indices_roi(labels_area, visparc_array):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (numpy.array): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        np.logical_or(visparc_array ==
                      labels_area[0], visparc_array == labels_area[1])
    )[0]
    return indices_area


def calc_msd(param_area1, param_area2):
    # mean squared difference between two different model results in an area
    n1 = param_area1.size
    n2 = param_area2.size
    assert (
        n1 == n2
    ), "areas don't have the same number of vertices and cannot be compared!"

    return np.sum((param_area1 - param_area2) ** 2) / n1


# Now `results` dictionary contains the output of `calc_msd` organized by subject, hemisphere, model combination, and parameter
# You can access the results using appropriate keys corresponding to each level of the nested structure
subject_info = pd.read_csv(
    "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv")
path_derivatives = "/data/p_02915/dhcp_derivatives_SPOT/HCP-D"
columns = ["Subject_ID",
           "L_eccentricity_real_benson",
           "R_eccentricity_real_benson",
           "L_polarangle_real_benson",
           "R_polarangle_real_benson",
           "L_eccentricity_real_simulated",
           "R_eccentricity_real_simulated",
           "L_polarangle_real_simulated",
           "R_polarangle_real_simulated",]
msd_value = pd.DataFrame(columns=columns)
for index, row in subject_info.iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    sub = sub_id.replace('sub-', '')
    msd_value.at[index, "Subject_ID"] = sub

    retinotopy_dict = {}
    in_out = {
        "eccentricity": OUTPUT_ECC,
        "polarangle": OUTPUT_PANGLE,
    }

    for hemi in ["L", "R"]:
        hemi_results = {}  # Results for the current hemisphere
        # FORMAT PATHS FOR INPUT AND OUTPUT
        prefix_model = PREFIX_MODEL.format(
            root_dir=path_derivatives,
            sub=sub,
            hemi=hemi,
        )
        prefix_sub_template = PREFIX_SUB_TEMPLATE.format(
            root_dir=path_derivatives,
            hemi=hemi,
        )

        # LOAD DATA
        visparc = surface.load_surf_data(
            PATH_VISPARC.format(prefix_sub_template=prefix_sub_template)
        )
        indices_v1 = get_indices_roi(LABELS_V1, visparc)
        indices_v2 = get_indices_roi(LABELS_V2, visparc)

        retinotopy_dict[(hemi, "benson", "eccentricity")] = surface.load_surf_data(
            PATH_ECCENTRICITY.format(prefix_sub_template=prefix_sub_template)
        )

        retinotopy_dict[(hemi, "benson", "polarangle")] = surface.load_surf_data(
            PATH_ANGLE.format(prefix_sub_template=prefix_sub_template)
        )

        for model in ["real", "simulated"]:

            try:
                ccf_v0 = surface.load_surf_data(
                    PATH_v0.format(prefix_model=prefix_model, model=model)
                ).astype(int, copy=False)

                ccf_r = surface.load_surf_data(
                    PATH_r.format(prefix_model=prefix_model, model=model)
                )
            except ValueError:
                print(
                    f"No model results found under {PATH_v0.format(prefix_model=prefix_model, model=model)}."
                )
                print(" Skipping this...", flush=True)
                continue

            # get indices into V1 for each area V2 vertex
            # indices refer to V1 vertex index inside V1, not whole brain!
            centers_v2 = ccf_v0[indices_v2]

            # indices to voxels in v1 array that are the centers for each v2 voxel
            # indexv1_centers_v2 = indices_v1[centers_v2]
            for param in ["eccentricity", "polarangle"]:

                retinotopy_out = in_out[param].format(
                    prefix_model=prefix_model,
                    model=model,
                    threshold=str("0.0").replace(".", ""),
                )

                # restrict to V1
                ret_v1 = retinotopy_dict[(hemi, "benson", param)][indices_v1]

                # project data from V1 vertices to V2
                ret_v2 = ret_v1[centers_v2]

                # make empty wholebrain
                ret_wholeb = np.zeros(
                    retinotopy_dict[(hemi, "benson", param)].shape)

                # fill retinotopy values into v2
                ret_wholeb[indices_v2] = ret_v2

                # keep data for model comparisons later
                retinotopy_dict[(hemi, model, param)] = ret_wholeb

        # CALCULATE DIFFERENCES BETWEEN MODELS PER HEMI
        # calculate mean square differences between models for eccentricity and angle
        try:
            for model_1, model_2 in [("real", "benson"), ("real", "simulated")]:
                for param in ["eccentricity", "polarangle"]:
                    msd = calc_msd(
                        retinotopy_dict[(hemi, model_1, param)][indices_v2],
                        retinotopy_dict[(hemi, model_2, param)][indices_v2],
                    )
                    # Store the result for the current param
                    msd_value.at[index,
                                 f"{hemi}_{param}_{model_1}_{model_2}"] = msd

        except KeyError:
            print(
                f"Data missing for hemisphere {hemi}, not calculating differences for this..."
            )

msd_value.to_csv('/data/p_02915/SPOT/MSD_HCP_old.csv')
