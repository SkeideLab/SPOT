#!/usr/bin/env python3
"""
Create dataframe of per-subject mean ROI differences across hemispheres.
Saves output to: roi_diff_results.csv
"""

import os
import glob
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn import surface

# Get the full path of the script
file_path = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------
# ROI utility
# ----------------------------------------------------
def get_indices_roi(labels_area, visparc):
    data = visparc.agg_data()
    return np.where((data == labels_area[0]) | (data == labels_area[1]))[0]

# ----------------------------------------------------
# Constants / paths
# ----------------------------------------------------
VISPARC_PATH1 = (
    "{root_dir}/hcp_surface/{sub}/anat/"
    "{sub}_hemi-{hemi}_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii"
)
VISPARC_PATH2 = (
    "{root_dir}/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii"
)
LABELS_V2 = (2, 3)

FIRST_SUFFIX  = "*_desc-real_r2_first_cross.gii"
SECOND_SUFFIX = "*_desc-real_r2_second_cross.gii"

GROUP_TO_CSV = {
    "2nd": '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv',
    "3rd": '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv',
    "preterm": '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_no_drop_v2.csv',
    "fullterm": '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_no_drop_v2.csv',
    "adolescent": '/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv',
    "adult": '/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv',
}

PREFIX_MODEL = (
    "{root_dir}/ccfmodel_var/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
)
PREFIX_MODEL_2 = (
    "{root_dir}/ccfmodel_var/{sub}/"
    "{sub}_hemi-{hemi}_mesh-native_dens-native"
)

# ----------------------------------------------------
# Output container
# ----------------------------------------------------
results = {}

# ----------------------------------------------------
# Iterate over hemispheres + groups
# ----------------------------------------------------
for hemi in ["L", "R"]:
    print(f"\n=== Hemisphere: {hemi} ===")
    for group, csv_path in GROUP_TO_CSV.items():
        print(f"\n-- Group: {group} --")
        subject_info = pd.read_csv(csv_path)

        for _, row in subject_info.iterrows():

            if group in ("adolescent", "adult"):
                sub_id = row["sub_id"]
                prefix_model = PREFIX_MODEL_2.format(
                    root_dir="/data/p_02915/dhcp_derivatives_SPOT/HCP-D",
                    sub=sub_id,
                    hemi=hemi,
                )
                visparc = nib.load(
                    VISPARC_PATH1.format(
                        root_dir="/data/p_02915/dhcp_derivatives_SPOT/HCP-D",
                        sub=sub_id,
                        hemi=hemi,
                    ))
                roi_idx = get_indices_roi(LABELS_V2, visparc)

            else:
                sub_id = row["sub_id"]
                sess_id = row["sess_id"]
                sub = sub_id.replace("sub-", "")
                ses = sess_id.replace("ses-", "")

                root_dir = (
                    "/data/p_02915/dhcp_derivatives_SPOT/fetal"
                    if group in ("2nd", "3rd")
                    else "/data/p_02915/dhcp_derivatives_SPOT/Neonates"
                )

                prefix_model = PREFIX_MODEL.format(
                    root_dir=root_dir, sub=sub, ses=ses, hemi=hemi
                )

                visparc = nib.load(
                    VISPARC_PATH2.format(
                        root_dir=root_dir, sub=sub, ses=ses, hemi=hemi
                    ))
                roi_idx = get_indices_roi(LABELS_V2, visparc)

            # Model file matching
            patt_first = prefix_model + "_" + FIRST_SUFFIX
            patt_second = prefix_model + "_" + SECOND_SUFFIX

            first_files = glob.glob(patt_first) or glob.glob(prefix_model + FIRST_SUFFIX)
            second_files = glob.glob(patt_second) or glob.glob(prefix_model + SECOND_SUFFIX)

            if len(first_files) == 0 or len(second_files) == 0:
                continue

            try:
                arr_first = surface.load_surf_data(first_files[0]).astype(float)
                arr_second = surface.load_surf_data(second_files[0]).astype(float)
            except:
                continue

            if arr_first.shape != arr_second.shape:
                continue

            roi_first_all = arr_first[roi_idx]
            roi_second_all = arr_second[roi_idx]

            valid_mask = roi_first_all > 0

            roi_first = roi_first_all[valid_mask]
            roi_second = roi_second_all[valid_mask]

            if len(roi_first) == 0:
                continue

            roi_diff = np.mean(roi_first - roi_second)

            results.setdefault((group, sub_id), []).append(roi_diff)

print("\nAll subjects processed.")

# ----------------------------------------------------
# Average across hemispheres
# ----------------------------------------------------
subject_avg_diff = {}

for (group, sub_id), diffs in results.items():
    subject_avg_diff.setdefault(group, []).append(np.mean(diffs))

# ----------------------------------------------------
# Save flattened dataframe
# ----------------------------------------------------
rows = []
for group, vals in subject_avg_diff.items():
    for v in vals:
        rows.append([group, v])

df = pd.DataFrame(rows, columns=["Group", "MeanROI_diff"])
df.to_csv(f"{file_path}/roi_diff_results.csv", index=False)

print("\nSaved: roi_diff_results.csv")
