import nibabel as nib
import numpy as np

from nilearn import surface

sub = "CC00058XX09"
ses = "11300"
path_output_data = "/data/p_02495/dhcp_derivatives/dhcp_surface"

benson_ecc = "{path_output_data}/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native_desc-angleretinotbenson2014_seg.shape.gii"
benson_angle = "{path_output_data}/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native_desc-angleretinotbenson2014_seg.shape.gii"

outpath_foveal = "{path_output_data}/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native_desc-foveal_roi.shape.gii"
outpath_middle = "{path_output_data}/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native_desc-middle_roi.shape.gii"
outpath_peripheral = "{path_output_data}/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native_desc-peripheral_roi.shape.gii"
outpath_dorsal = "{path_output_data}/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native_desc-dorsal_roi.shape.gii"
outpath_ventral = "{path_output_data}/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native_desc-ventral_roi.shape.gii"

# division into foveal, middle, peripheral, ventral, dorsal
boundaries = {"L": (5.5, 18.9), "R": (5, 18.6)}

for hemi_upper in ["L", "R"]:
    ecc = surface.load_surf_data(
        benson_ecc.format(
            path_output_data=path_output_data, sub=sub, ses=ses, hemi_upper=hemi_upper
        )
    )
    pangle = surface.load_surf_data(
        benson_angle.format(
            path_output_data=path_output_data, sub=sub, ses=ses, hemi_upper=hemi_upper
        )
    )
    mask_foveal = (ecc > 0) & (ecc <= boundaries[hemi_upper][0])
    mask_middle = (ecc > boundaries[hemi_upper][0]) & (ecc <= boundaries[hemi_upper][1])
    mask_peripheral = ecc > boundaries[hemi_upper][1]

    # TODO the masks have very different numbers of vertices (359, 977,11996).
    # Is this because there is no restriction to V2? change boundaries?
    # assert (
    #     np.abs(np.sum(mask_foveal) - np.sum(mask_middle)) < 10
    # ), "Middle and foveal mask have a difference of more than 10 vertices"
    # assert (
    #     np.abs(np.sum(mask_peripheral) - np.sum(mask_middle)) < 10
    # ), "Middle and peripheral mask have a difference of more than 10 vertices"

    mask_ventral = pangle > 90
    mask_dorsal = (pangle < 90) & (pangle > 0)

    for data, outpath in [
        (mask_foveal, outpath_foveal),
        [mask_middle, outpath_middle],
        (mask_peripheral, outpath_peripheral),
        (mask_ventral, outpath_ventral),
        (mask_dorsal, outpath_dorsal),
    ]:
        img = nib.gifti.GiftiImage(darrays=[nib.gifti.GiftiDataArray(np.int32(data))])
        nib.save(
            img,
            outpath.format(
                path_output_data=path_output_data,
                sub=sub,
                ses=ses,
                hemi_upper=hemi_upper,
            ),
        )
        print(
            f"Saving mask under {outpath.format(path_output_data=path_output_data,sub=sub,ses=ses,hemi_upper=hemi_upper)}..."
        )

    # TODO restrict to V2?? based on benson or wang??
