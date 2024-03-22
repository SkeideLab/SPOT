import nibabel as nib
import numpy as np

from nilearn import surface

sub = "CC00058XX09"
ses = "11300"
path_output_data = "/data/p_02495/dhcp_derivatives/dhcp_surface"
prefix_alldata = "{path_output_data}/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi_upper}_mesh-native_dens-native"

benson_ecc = "_desc-angleretinotbenson2014_seg.shape.gii"
benson_angle = "_desc-angleretinotbenson2014_seg.shape.gii"
benson_visareas = "_desc-retinotbenson2014_label-visarea_dparc.label.gii"

outpath_foveal = "_desc-foveal_roi.shape.gii"
outpath_middle = "_desc-middle_roi.shape.gii"
outpath_peripheral = "_desc-peripheral_roi.shape.gii"
outpath_dorsal = "_desc-dorsal_roi.shape.gii"
outpath_ventral = "_desc-ventral_roi.shape.gii"

# division into foveal, middle, peripheral, ventral, dorsal subregions
# from Bock 2015
boundaries = {"L": (5.5, 18.9), "R": (5, 18.6)}

for hemi in ["L", "R"]:

    prefix_thisdata = prefix_alldata.format(
        path_output_data=path_output_data, sub=sub, ses=ses, hemi_upper=hemi
    )
    ecc = surface.load_surf_data(prefix_thisdata + benson_ecc)
    pangle = surface.load_surf_data(prefix_thisdata + benson_angle)
    visareas = surface.load_surf_data(prefix_thisdata + benson_visareas)

    # restrict all masks to only V2 vertices
    mask_v2 = visareas == 2  # | (visareas == 3)  # 2 is code for V2
    mask_foveal = (ecc > 0) & (ecc <= boundaries[hemi][0]) & mask_v2
    mask_middle = (ecc > boundaries[hemi][0]) & (ecc <= boundaries[hemi][1]) & mask_v2
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
        (mask_foveal, outpath_foveal),
        [mask_middle, outpath_middle],
        (mask_peripheral, outpath_peripheral),
        (mask_ventral, outpath_ventral),
        (mask_dorsal, outpath_dorsal),
    ]:
        img = nib.gifti.GiftiImage(darrays=[nib.gifti.GiftiDataArray(np.int32(data))])
        nib.save(img, prefix_thisdata + outpath)
        print(f"Saving mask under {prefix_thisdata+outpath}...")
