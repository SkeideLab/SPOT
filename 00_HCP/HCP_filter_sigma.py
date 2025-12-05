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

PREFIX_MODEL = (
    "{root_dir}/ccfmodel/{sub}/"
    "{sub}_hemi-{hemi}_mesh-native_dens-native"
)
PREFIX_SUB_TEMPLATE = (
    "{root_dir}/"
    "hemi-{hemi}_mesh-native_dens-native"
)

PATH_r = "{prefix_model}_desc-{model}_r.gii"
PATH_sigma = "{prefix_model}_desc-{model}_sigma.gii"

PATH_VISPARC = (
    "{prefix_sub_template}_desc-retinotbenson2014_label-visarea_dparc.label.gii"
)
LABELS_V1 = [1]
LABELS_V2 = [2, 3]
OUTPATH_sigma = "{prefix_model}_desc-{model}_sigma_{threshold}.gii"
OUTPATH_r = "{prefix_model}_desc-{model}_r_{threshold}.gii"

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


def main():
    args = parse_args()
    # -d /data/p_02495/dhcp_derivatives -sub CC00065XX08 -ses 18600
    for hemi in ["L", "R"]:
        # FORMAT PATHS FOR INPUT AND OUTPUT
        prefix_model = PREFIX_MODEL.format(
            root_dir=args.derivatives_directory,
            sub=args.sub,
            hemi=hemi,
        )
        prefix_sub_template = PREFIX_SUB_TEMPLATE.format(
            root_dir="/data/p_02915/SPOT/00_HCP/template",
            sub=args.sub,
            hemi=hemi,
        )
        for model in ["real", "simulated"]:

            ccf_r = surface.load_surf_data(
                        PATH_r.format(prefix_model=prefix_model, model=model)
                    )
            ccf_sigma = surface.load_surf_data(
                        PATH_sigma.format(prefix_model=prefix_model, model=model)
                    )

            # LOAD DATA
            visparc = nib.load(PATH_VISPARC.format(prefix_sub_template=prefix_sub_template))
            indices_v2 = get_indices_roi(LABELS_V2, visparc)

            # fill only those vertices that pass the correlation threshold
            indices_v2_success = np.intersect1d(
                indices_v2,
                np.nonzero(np.abs(ccf_r) > float(args.threshold)),
                assume_unique=True,
            )
            ret_wholeb = np.zeros(
                        ccf_r.shape)
                # get indices into V1 for each area V2 vertex
                # indices refer to V1 vertex index inside V1, not whole brain!
            ret_wholeb[indices_v2_success] = ccf_sigma[indices_v2_success]

            # save as gifti
            img_gifti = nib.gifti.GiftiImage(
                darrays=[
                    nib.gifti.GiftiDataArray(
                        np.float32(np.squeeze(ret_wholeb)))
                ]
            )
            retinotopy_out = OUTPATH_sigma.format(prefix_model=prefix_model, model=model, threshold = str(args.threshold).replace(".", ""))
            nib.save(img_gifti, retinotopy_out)
            print(f"Saving {retinotopy_out}...")

            ret_wholeb = np.zeros(
                        ccf_r.shape)

            ret_wholeb[indices_v2_success] = ccf_r[indices_v2_success]

            # save as gifti
            img_gifti = nib.gifti.GiftiImage(
                darrays=[
                    nib.gifti.GiftiDataArray(
                        np.float32(np.squeeze(ret_wholeb)))
                ]
            )
            retinotopy_out = OUTPATH_r.format(prefix_model=prefix_model, model=model, threshold = str(args.threshold).replace(".", ""))

            nib.save(img_gifti, retinotopy_out)
            print(f"Saving {retinotopy_out}...")


if __name__ == "__main__":
    main()
