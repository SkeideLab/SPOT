"""Sanity check for surface-projected data: does a seed correlation based on 
precuneus or posterior cingulate cortex show the default mode network?
"""
from pathlib import Path
import utils
import nibabel as nib
import numpy as np
from nilearn import surface
from scipy import stats
from nilearn import plotting

base_dir = Path("/data/p_02495")
sub_ses = {
    (utils.get_sub(str(item)), utils.get_ses(str(item)))
    for item in base_dir.glob("dhcp*/**")
    if "sub-" in str(item) and "ses-" in str(item)
}
work_dir = base_dir / "dhcp_intermed"

for sub, ses in sub_ses:
    for hemi in ["left", "right"]:
        func_surf_path = (
            base_dir
            / "dhcp_intermed"
            / f"sub-{sub}"
            / f"ses-{ses}"
            / "surface"
            / f"func_hemi-{hemi}_mesh-native.func.gii"
        )
        if not func_surf_path.exists():
            print("error")
            continue

        parcel_surf_path = list(
            (
                base_dir / "dhcp_anat_pipeline" / f"sub-{sub}" / f"ses-{ses}" / "anat"
            ).glob(f"*_hemi-{hemi[0].upper()}_desc-drawem_space-T2w_dparc.dlabel.gii")
        )[0]

        pial_path = list(
            (
                base_dir / "dhcp_anat_pipeline" / f"sub-{sub}" / f"ses-{ses}" / "anat"
            ).glob(f"*_hemi-{hemi[0].upper()}_space-T2w_pial.surf.gii")
        )[0]

        sulc_path = list(
            (
                base_dir / "dhcp_anat_pipeline" / f"sub-{sub}" / f"ses-{ses}" / "anat"
            ).glob(f"*_hemi-{hemi[0].upper()}_space-T2w_sulc.shape.gii")
        )[0]

        func_surf = surface.load_surf_data(func_surf_path)
        parcel_surf = surface.load_surf_data(parcel_surf_path)
        parcels = [34, 35]  # labels for L and R posterior cingulate gyrus
        mask = np.logical_or(*[parcel_surf == label for label in parcels])
        seed_ts = np.mean(func_surf[mask], axis=0)

        corr_img = np.zeros(func_surf.shape[0])

        for i in range(func_surf.shape[0]):
            corr_img[i] = stats.pearsonr(seed_ts, func_surf[i])[0]

        plotting.plot_surf_stat_map(
            pial_path,
            stat_map=corr_img,
            hemi=hemi,
            view="medial",
            colorbar=True,
            bg_map=sulc_path,
            bg_on_data=True,
            darkness=0.3,
            title="Correlation map",
            output_file=f"/data/u_kieslinger_software/code/flux/_outputs_sanitycheck/sub-{sub}_ses-{ses}_hemi-{hemi}_medial.png",
        )
        plotting.plot_surf_stat_map(
            pial_path,
            stat_map=corr_img,
            hemi=hemi,
            view="lateral",
            colorbar=True,
            bg_map=sulc_path,
            bg_on_data=True,
            darkness=0.3,
            title="Correlation map",
            output_file=f"/data/u_kieslinger_software/code/flux/_outputs_sanitycheck/sub-{sub}_ses-{ses}_hemi-{hemi}_lateral.png",
        )
