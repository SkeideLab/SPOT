import nibabel as nib
from nilearn import surface
import numpy as np


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


model = "lc"  # 'lc' or 'ccf'
root_dir = "/data/p_02495/dhcp_derivatives"
file_v0 = "/data/p_02495/dhcp_derivatives/ccfmodel/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}_desc-{model}_v0i.gii"
file_eccentricity = "/data/p_02495/dhcp_derivatives/dhcp_surface/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-eccentretinotbenson2014_seg.shape.gii"
file_angle = "/data/p_02495/dhcp_derivatives/dhcp_surface/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-angleretinotbenson2014_seg.shape.gii"
file_visparc = "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/space-T2w/anat/sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
output_ecc = "/data/p_02495/dhcp_derivatives/ccfmodel/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}_label-eccentricity_desc-{model}_roi-v2_metric.gii"
output_angle = "/data/p_02495/dhcp_derivatives/ccfmodel/sub-{sub}/ses-{ses}/sub-{sub}_ses-{ses}_hemi-{hemi}_label-angle_desc-{model}_roi-v2_metric.gii"

sub = "CC00058XX09"
ses = "11300"

labels_v1 = [1, 2]  # see `ROIfiles_Labeling.txt` in this directory
labels_v2 = [3, 4]
for hemi in ["L"]:
    visparc = surface.load_surf_data(
        file_visparc.format(root_dir=root_dir, sub=sub, ses=ses, hemi=hemi)
    )
    indices_v1 = get_indices_roi(labels_v1, visparc)
    indices_v2 = get_indices_roi(labels_v2, visparc)
    ccf_v0 = surface.load_surf_data(
        file_v0.format(sub=sub, ses=ses, hemi=hemi, model=model)
    ).astype(int, copy=False)

    centers_v2 = ccf_v0[indices_v2]

    in_out = [
        (
            file_eccentricity.format(sub=sub, ses=ses, hemi=hemi),
            output_ecc.format(sub=sub, ses=ses, hemi=hemi, model=model),
        ),
        (
            file_angle.format(sub=sub, ses=ses, hemi=hemi),
            output_angle.format(sub=sub, ses=ses, hemi=hemi, model=model),
        ),
    ]

    # indices to voxels in v1 array that are the centers for each v2 voxel
    # indexv1_centers_v2 = indices_v1[centers_v2]
    for retinotopy_in, retinotopy_out in in_out:
        ret = surface.load_surf_data(retinotopy_in)
        ret_v1 = ret[indices_v1]
        ret_v2 = ret_v1[centers_v2]
        ret_wholeb = np.zeros(ret.shape)
        ret_wholeb[indices_v2] = ret_v2

        img_gifti = nib.gifti.GiftiImage(
            darrays=[nib.gifti.GiftiDataArray(np.float32(np.squeeze(ret_wholeb)))]
        )
        nib.save(img_gifti, retinotopy_out)
        print(f"Saving {retinotopy_out}...")
