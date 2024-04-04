"""
Project V1 retinotopy template values to V2 based on cortico-connective field results.
"""

from argparse import ArgumentParser

import nibabel as nib
import numpy as np

from nilearn import surface

# FIXED INPUT AND OUTPUT STRUCTURE
PREFIX_MODEL = "{root_dir}/ccfmodel/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}"
PREFIX_SUB_TEMPLATE = (
    "{root_dir}/dhcp_surface/"
    "sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
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


def get_indices_roi(labels_area, visparc_array):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (numpy.array): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        np.logical_or(visparc_array == labels_area[0], visparc_array == labels_area[1])
    )[0]
    return indices_area


def main():
    args = parse_args()
    # -d /data/p_02495/dhcp_derivatives -sub CC00065XX08 -ses 18600

    for model in ["real", "simulated"]:

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

            # list of corresponding input/output transformations
            in_out = [
                (
                    PATH_ECCENTRICITY.format(prefix_sub_template=prefix_sub_template),
                    OUTPUT_ECC.format(
                        prefix_model=prefix_model,
                        model=model,
                        threshold=str(args.threshold).replace(".", ""),
                    ),
                ),
                (
                    PATH_ANGLE.format(prefix_sub_template=prefix_sub_template),
                    OUTPUT_PANGLE.format(
                        prefix_model=prefix_model,
                        model=model,
                        threshold=str(args.threshold).replace(".", ""),
                    ),
                ),
            ]

            # LOAD DATA
            visparc = surface.load_surf_data(
                PATH_VISPARC.format(prefix_sub_template=prefix_sub_template)
            )

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

            indices_v1 = get_indices_roi(LABELS_V1, visparc)
            indices_v2 = get_indices_roi(LABELS_V2, visparc)

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
            for retinotopy_in, retinotopy_out in in_out:

                # load template data in sub space
                ret = surface.load_surf_data(retinotopy_in)
                # TODO check if these really are the indices into just  v1
                # restrict to V1
                ret_v1 = ret[indices_v1]

                # project data from V1 vertices to V2
                ret_v2 = ret_v1[centers_v2]

                # make empty wholebrain
                ret_wholeb = np.zeros(ret.shape)

                # fill retinotopy values into v2
                ret_wholeb[indices_v2_success] = ret_v2

                # save as gifti
                img_gifti = nib.gifti.GiftiImage(
                    darrays=[
                        nib.gifti.GiftiDataArray(np.float32(np.squeeze(ret_wholeb)))
                    ]
                )

                nib.save(img_gifti, retinotopy_out)
                print(f"Saving {retinotopy_out}...")


if __name__ == "__main__":
    main()
