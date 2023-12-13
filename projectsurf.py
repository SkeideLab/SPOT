from dhcp.func import surface
from dhcp.util import enums, mask
from pathlib import Path
import nipype.interfaces.fsl as fsl
import utils


base_dir = Path("/data/p_02495")
sub_ses = {
    (utils.get_sub(str(item)), utils.get_ses(str(item)))
    for item in base_dir.glob("dhcp*/**")
    if "sub-" in str(item) and "ses-" in str(item)
}
work_dir = base_dir / "dhcp_intermed"


for subid, sesid in sub_ses:
    print(f"Working on sub {subid} ses {sesid}")
    if (
        work_dir
        / f"sub-{subid}"
        / f"ses-{sesid}"
        / "surface"
        / "func_hemi-left_mesh-native.func.gii"
    ).exists():
        print(f"sub {subid} ses {sesid} already has functional surface data")
        continue

    anat_path = (
        base_dir / "dhcp_anat_pipeline" / f"sub-{subid}" / f"ses-{sesid}" / "anat"
    )
    functional_path = (
        base_dir / "dhcp_fmri_pipeline" / f"sub-{subid}" / f"ses-{sesid}" / "func"
    )

    # create mask for subcortical structures
    scmask_path = (
        work_dir
        / f"sub-{subid}"
        / f"ses-{sesid}"
        / "anat"
        / f"sub-{subid}_ses-{sesid}_desc-sc_space-T2w_mask.nii.gz"
    )

    if not scmask_path.exists():
        scmask_path.parent.mkdir(parents=True, exist_ok=True)
        dseg_type = enums.SegType(3)
        mask.create_mask(
            dseg=str(
                anat_path
                / f"sub-{subid}_ses-{sesid}_desc-drawem9_space-T2w_dseg.nii.gz"
            ),
            dseg_type=dseg_type,
            labels=["sc", "cb", "bs"],
            outname=str(scmask_path),
        )

    # invert bold to t2 affine
    transform_path = functional_path.parent / "xfm"
    if not transform_path.exists():
        print(f"no transform found for sub {subid} ses {sesid}")
        continue
    from_bold_to_t2w = str(list(transform_path.glob("*from-bold_to-T2w*"))[0])
    from_t2w_to_bold = from_bold_to_t2w.replace("from-bold_to-T2w", "from-T2w_to-bold")

    if not Path(from_t2w_to_bold).exists():
        invt = fsl.ConvertXFM()
        invt.inputs.in_file = from_bold_to_t2w
        invt.inputs.invert_xfm = True
        invt.inputs.out_file = from_t2w_to_bold

        invt.run()

    surface.sample_to_native(
        struct_white_surf_right=str(
            anat_path / f"sub-{subid}_ses-{sesid}_hemi-R_space-T2w_wm.surf.gii"
        ),
        struct_white_surf_left=str(
            anat_path / f"sub-{subid}_ses-{sesid}_hemi-L_space-T2w_wm.surf.gii"
        ),
        struct_pial_surf_right=str(
            anat_path / f"sub-{subid}_ses-{sesid}_hemi-R_space-T2w_pial.surf.gii"
        ),
        struct_pial_surf_left=str(
            anat_path / f"sub-{subid}_ses-{sesid}_hemi-L_space-T2w_pial.surf.gii"
        ),
        struct_mid_surf_right=str(
            anat_path
            / f"sub-{subid}_ses-{sesid}_hemi-R_space-T2w_midthickness.surf.gii"
        ),
        struct_mid_surf_left=str(
            anat_path
            / f"sub-{subid}_ses-{sesid}_hemi-L_space-T2w_midthickness.surf.gii"
        ),
        medialwall_shape_right=str(
            anat_path / f"sub-{subid}_ses-{sesid}_hemi-R_desc-medialwall_mask.shape.gii"
        ),
        medialwall_shape_left=str(
            anat_path / f"sub-{subid}_ses-{sesid}_hemi-L_desc-medialwall_mask.shape.gii"
        ),
        struct_to_func_affine=from_t2w_to_bold,
        func=str(
            functional_path
            / f"sub-{subid}_ses-{sesid}_task-rest_desc-preproc_bold.nii.gz"
        ),
        func_brainmask=str(
            functional_path
            / f"sub-{subid}_ses-{sesid}_task-rest_desc-preproc_space-bold_brainmask.nii.gz"
        ),
        struct=str(anat_path / f"sub-{subid}_ses-{sesid}_desc-restore_T2w.nii.gz"),
        workdir=str(work_dir / f"sub-{subid}" / f"ses-{sesid}" / "surface"),
        struct_subcortical_mask=str(scmask_path),
        do_sigloss=True,
    )
