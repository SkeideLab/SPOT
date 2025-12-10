"""
Cacluate MSD values
"""

from argparse import ArgumentParser
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import surface

# FIXED INPUT AND OUTPUT STRUCTURE
PREFIX_MODEL = (
    "{root_dir}/ccfmodel_var/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
) 
PREFIX_SUB_TEMPLATE = (
    "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
)

PATH_v0 = "{prefix_model}_desc-{model}_v0i.gii"
PATH_r = "{prefix_model}_desc-{model}_r.gii"
OUTPUT_ECC = (
    "{prefix_model}_label-eccentricity_desc-{model}_roi-v2th00_metric.gii"
)
OUTPUT_PANGLE = (
    "{prefix_model}_label-polarangle_desc-{model}_roi-v2th00_metric.gii"
)
OUTPUT_SIGMA = (
    "{prefix_model}_desc-{model}_sigma.gii"
)

PATH_ECCENTRICITY = "{prefix_sub_template}_desc-eccentretinotbenson2014_seg.shape.gii"
PATH_ANGLE = "{prefix_sub_template}_desc-angleretinotbenson2014_seg.shape.gii"
PATH_SIGMA = "{prefix_sub_template}_desc-sigmaretinotbenson2014_seg.shape.gii"
PATH_VISPARC = (
    "{prefix_sub_template}_desc-retinotbenson2014_label-visarea_dparc.label.gii"
)
LABELS_V1 = [1]
LABELS_V2 = [2, 3]


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
        "-ses",
        required=True,
        help="Session being processed",
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
subject_info = pd.read_csv("/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv")
path_derivatives = "/data/p_02915/dhcp_derivatives_SPOT/fetal"
columns = ["Subject_ID",
           #"L_sigma_real_benson",
           #"R_sigma_real_benson",
           #"L_sigma_real_simulated",
           #"R_sigma_real_simulated",
           "L_eccentricity_real_benson", 
           "R_eccentricity_real_benson", 
           "L_polarangle_real_benson", 
           "R_polarangle_real_benson",
           "L_eccentricity_real_simulated", 
           "R_eccentricity_real_simulated", 
           "L_polarangle_real_simulated", 
           "R_polarangle_real_simulated",
           ]
msd_value = pd.DataFrame(columns=columns)
for index, row in subject_info.iterrows():
    sub_id = subject_info.at[index, "sub_id"]
    sess_id = subject_info.at[index, "sess_id"]
    sub = sub_id.replace('sub-','')
    ses = sess_id.replace('ses-','')
    msd_value.at[index,"Subject_ID"] = sub

    retinotopy_dict = {}
    in_out = {
        "eccentricity": OUTPUT_ECC,
        "polarangle": OUTPUT_PANGLE,
        #"sigma": OUTPUT_SIGMA,
    }    

    for hemi in ["L", "R"]:
        hemi_results = {}  # Results for the current hemisphere
        # FORMAT PATHS FOR INPUT AND OUTPUT
        prefix_model = PREFIX_MODEL.format(
            root_dir=path_derivatives,
            sub=sub,
            ses=ses,
            hemi=hemi,
        )
        prefix_sub_template = PREFIX_SUB_TEMPLATE.format(
            root_dir=path_derivatives,
            sub=sub,
            ses=ses,
            hemi=hemi,
        )

        # LOAD DATA
        visparc = nib.load(
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

            retinotopy_dict[(hemi, model, 'eccentricity')] = surface.load_surf_data(
                OUTPUT_ECC.format(prefix_model=prefix_model, model = model)
            )
            retinotopy_dict[(hemi, model, "polarangle")] = surface.load_surf_data(
               OUTPUT_PANGLE.format(prefix_model=prefix_model, model=model)
            )

        try:
            for model_1, model_2 in [("real", "benson"), ("real", "simulated")]:
                for param in ["eccentricity","polarangle"]:
                    msd = calc_msd(
                        retinotopy_dict[(hemi, model_1, param)][indices_v2],
                        retinotopy_dict[(hemi, model_2, param)][indices_v2],
                    )
                    msd_value.at[index,f"{hemi}_{param}_{model_1}_{model_2}"] = msd  # Store the result for the current param  
                
        except KeyError:
            print(
                f"Data missing for hemisphere {hemi}, not calculating differences for this..."
            )        
       
msd_value.to_csv('/data/p_02915/SPOT/MSD_fetal_young.csv')