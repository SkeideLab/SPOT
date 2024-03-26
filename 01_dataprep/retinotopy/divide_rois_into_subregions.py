"""Splits V2 into several subregions (foveal/middle/peripheral; dorsal/ventral).
Generates gifti files for each subregion with value 1 at each node that lies in the region and 0 otherwise.
Surface based, T2w space, native mesh.
Uses Benson2014 template to define subregions, based on eccentricity and polar angle.

Subregion definitions:
----------------------------
foveal: 0 < eccentricity <= (5.5 or 5)
middle: (5.5 or 5) < eccentricity <= (18.9 or 18.6)
peripheral: eccentricity > (18.9 or 18.6)
dorsal: polar angle < 90
ventral: polar angle > 90
-----------------------------

For more info see Bock 2015 JNeurosci"""

from argparse import ArgumentParser

import nibabel as nib
import numpy as np

from nilearn import surface

PREFIX_ALLDATA = "{path_output_data}/sub-{sub}/ses-{ses}/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native"

BENSON_ECC = "_desc-angleretinotbenson2014_seg.shape.gii"
BENSON_ANGLE = "_desc-angleretinotbenson2014_seg.shape.gii"
BENSON_VISAREAS = "_desc-retinotbenson2014_label-visarea_dparc.label.gii"

OUTSUFFIX_FOVEAL = "_desc-foveal_roi.shape.gii"
OUTSUFFIX_MIDDLE = "_desc-middle_roi.shape.gii"
OUTSUFFIX_PERIPHERAL = "_desc-peripheral_roi.shape.gii"
OUTSUFFIX_DORSAL = "_desc-dorsal_roi.shape.gii"
OUTSUFFIX_VENTRAL = "_desc-ventral_roi.shape.gii"


def parse_args():
    """Parses arguments from the command line."""

    parser = ArgumentParser()
    parser.add_argument(
        "-o",
        "--output_directory",
        required=True,
        help="Directory with subject-level output datasets (e.g. 'dhcp_surface')",
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
    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    # division into foveal, middle, peripheral, ventral, dorsal subregions
    # from Bock 2015
    boundaries = {"L": (5.5, 18.9), "R": (5, 18.6)}

    for hemi in ["L", "R"]:

        prefix_thisdata = PREFIX_ALLDATA.format(
            path_output_data=args.output_directory,
            sub=args.sub,
            ses=args.ses,
            hemi_upper=hemi,
        )
        ecc = surface.load_surf_data(prefix_thisdata + BENSON_ECC)
        pangle = surface.load_surf_data(prefix_thisdata + BENSON_ANGLE)
        visareas = surface.load_surf_data(prefix_thisdata + BENSON_VISAREAS)

        # restrict all masks to only V2 vertices
        mask_v2 = visareas == 2  # | (visareas == 3)  # 2 is code for V2
        mask_foveal = (ecc > 0) & (ecc <= boundaries[hemi][0]) & mask_v2
        mask_middle = (
            (ecc > boundaries[hemi][0]) & (ecc <= boundaries[hemi][1]) & mask_v2
        )
        mask_peripheral = (ecc > boundaries[hemi][1]) & mask_v2

        # TODO the masks have very different numbers of vertices (359, 977,11996).
        # Is this because there is no restriction to V2? change boundaries?
        # No, restriction to V2/V2+V3 does not change anything.
        # The reason must be different numbers of vertices on native vs fsaverage mesh in these regions.
        # Change the boundaries so that the numbers are similar? not sure

        for region, region_mask in [
            ("foveal", mask_foveal),
            ("middle", mask_middle),
            ("peripheral", mask_peripheral),
        ]:
            print(f"Size of {region} mask is: {region_mask.sum()} vertices.")

        # assert (
        #     np.abs(np.sum(mask_foveal) - np.sum(mask_middle)) < 10
        # ), "Middle and foveal mask have a difference of more than 10 vertices"
        # assert (
        #     np.abs(np.sum(mask_peripheral) - np.sum(mask_middle)) < 10
        # ), "Middle and peripheral mask have a difference of more than 10 vertices"

        mask_ventral = (pangle > 90) & mask_v2
        mask_dorsal = (pangle < 90) & (pangle > 0) & mask_v2

        for data, outpath in [
            (mask_foveal, OUTSUFFIX_FOVEAL),
            [mask_middle, OUTSUFFIX_MIDDLE],
            (mask_peripheral, OUTSUFFIX_PERIPHERAL),
            (mask_ventral, OUTSUFFIX_VENTRAL),
            (mask_dorsal, OUTSUFFIX_DORSAL),
        ]:
            img = nib.gifti.GiftiImage(
                darrays=[nib.gifti.GiftiDataArray(np.int32(data))]
            )
            nib.save(img, prefix_thisdata + outpath)
            print(f"Saving mask under {prefix_thisdata+outpath}...")


if __name__ == "__main__":
    main()
