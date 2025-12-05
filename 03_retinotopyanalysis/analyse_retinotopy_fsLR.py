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

from nilearn import surface

# FIXED INPUT AND OUTPUT STRUCTURE
PREFIX_MODEL = (
    "{root_dir}/ccfmodel_fsLR/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
)
PREFIX_SUB_TEMPLATE = (
    "/data/p_02915/SPOT/00_HCP/template/"
    "hemi-{hemi}_mesh-native_dens-native"
)

PATH_v0 = "{prefix_model}_desc-{model}_v0i.gii"
PATH_r = "{prefix_model}_desc-{model}_r.gii"
PATH_sigma = "{prefix_model}_desc-{model}_sigma.gii"
OUTPUT_ECC = (
    "{prefix_model}_label-eccentricity_desc-{model}_roi-v2th{threshold}_metric.gii"
)
OUTPUT_PANGLE = (
    "{prefix_model}_label-polarangle_desc-{model}_roi-v2th{threshold}_metric.gii"
)


PATH_ECCENTRICITY = "{prefix_sub_template}_desc-eccentretinotbenson2014_seg.shape.gii"
PATH_ANGLE = "{prefix_sub_template}_desc-angleretinotbenson2014_seg.shape.gii"
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


def print_differences(retinotopy_dict, hemi, indices_v2):
    for model_1, model_2 in [("real", "benson"), ("real", "simulated")]:
        for param in ["eccentricity", "polarangle"]:
            msd = calc_msd(
                retinotopy_dict[(hemi, model_1, param)][indices_v2],
                retinotopy_dict[(hemi, model_2, param)][indices_v2],
            )
            print(f"MSD for {param}: comparison {model_1} vs {model_2}: {msd}")


def main():
    args = parse_args()
    # -d /data/p_02495/dhcp_derivatives -sub CC00065XX08 -ses 18600

    retinotopy_dict = {}
    in_out = {
        "eccentricity": OUTPUT_ECC,
        "polarangle": OUTPUT_PANGLE,
    }
    for hemi in ["L", "R"]:
        # FORMAT PATHS FOR INPUT AND OUTPUT
        prefix_model = PREFIX_MODEL.format(
            root_dir=args.derivatives_directory,
            sub=args.sub,
            ses=args.ses,
            hemi=hemi,
        )
        prefix_sub_template = PREFIX_SUB_TEMPLATE.format(
            root_dir=args.derivatives_directory,
            sub=args.sub,
            ses=args.ses,
            hemi=hemi,
        )

        # LOAD DATA
        visparc = nib.load(PATH_VISPARC.format(prefix_sub_template=prefix_sub_template))
        
        indices_v1 = get_indices_roi(LABELS_V1, visparc)
        indices_v2 = get_indices_roi(LABELS_V2, visparc)

        retinotopy_dict[(hemi, "benson", "eccentricity")] = surface.load_surf_data(
            PATH_ECCENTRICITY.format(prefix_sub_template=prefix_sub_template)
        )

        retinotopy_dict[(hemi, "benson", "polarangle")] = surface.load_surf_data(
            PATH_ANGLE.format(prefix_sub_template=prefix_sub_template)
        )

        for model in ["real"]:

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

            # fill only those vertices that pass the correlation threshold
            indices_v2_success = np.intersect1d(
                indices_v2,
                np.nonzero(np.abs(ccf_r) > float(args.threshold)),
                assume_unique=True,
            )
            # get indices into V1 for each area V2 vertex
            # indices refer to V1 vertex index inside V1, not whole brain!
            centers_v2 = ccf_v0[indices_v2_success]

            # indices to voxels in v1 array that are the centers for each v2 voxel
            # indexv1_centers_v2 = indices_v1[centers_v2]
            for param in ["eccentricity", "polarangle"]:

                retinotopy_out = in_out[param].format(
                    prefix_model=prefix_model,
                    model=model,
                    threshold=str(args.threshold).replace(".", ""),
                )

                # restrict to V1
                ret_v1 = retinotopy_dict[(hemi, "benson", param)][indices_v1]

                # project data from V1 vertices to V2
                ret_v2 = ret_v1[centers_v2]

                # make empty wholebrain
                ret_wholeb = np.zeros(
                    retinotopy_dict[(hemi, "benson", param)].shape)

                # fill retinotopy values into v2
                ret_wholeb[indices_v2_success] = ret_v2

                # keep data for model comparisons later
                retinotopy_dict[(hemi, model, param)] = ret_wholeb

            
                # save as gifti
                img_gifti = nib.gifti.GiftiImage(
                    darrays=[
                        nib.gifti.GiftiDataArray(
                            np.float32(np.squeeze(ret_wholeb)))
                    ]
                )

                nib.save(img_gifti, retinotopy_out)
                print(f"Saving {retinotopy_out}...")


        # CALCULATE DIFFERENCES BETWEEN MODELS PER HEMI
        # calculate mean square differences between models for eccentricity and angle
        # TODO decide how to deal with the thresholded vertices
        # might result in different numbers of 'active' vertices for real vs simulated data
        # can these two datasets be compared then?
        # 'inactive' vertices get filled with 0, this could bias the msd
        # use intersection of active vertices instead?
        #try:
        #    print_differences(retinotopy_dict, hemi, indices_v2)
        #except KeyError:
        #    print(
        #        f"Data missing for hemisphere {hemi}, not calculating differences for this..."
        #    )


if __name__ == "__main__":
    main()
