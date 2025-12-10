import os
import glob
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn import surface
from scipy.stats import wilcoxon,rankdata

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
    "preterm": '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_no_drop_v2.csv',
    "fullterm": '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_no_drop_v2.csv',
    "2nd": '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv',
    "3rd": '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv',
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

        for index, row in subject_info.iterrows():

            # ------------------------------------------------
            # Build prefix path + ROI
            # ------------------------------------------------
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
                print(f"ROI size = {len(roi_idx)} vertices")

            else:
                sub_id = row["sub_id"]
                sess_id = row["sess_id"]
                sub = sub_id.replace("sub-", "")
                ses = sess_id.replace("ses-", "")

                if group in ("2nd", "3rd"):
                    root_dir = "/data/p_02915/dhcp_derivatives_SPOT/fetal"
                else:
                    root_dir = "/data/p_02915/dhcp_derivatives_SPOT/Neonates"

                prefix_model = PREFIX_MODEL.format(
                    root_dir=root_dir,
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )

                visparc = nib.load(
                    VISPARC_PATH2.format(
                        root_dir=root_dir,
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                    ))
                roi_idx = get_indices_roi(LABELS_V2, visparc)
                print(f"ROI size = {len(roi_idx)} vertices")

            # ------------------------------------------------
            # Find FIRST & SECOND map files
            # ------------------------------------------------
            patt_first  = prefix_model + "_" + FIRST_SUFFIX
            patt_second = prefix_model + "_" + SECOND_SUFFIX

            first_files  = glob.glob(patt_first)  or glob.glob(prefix_model + FIRST_SUFFIX)
            second_files = glob.glob(patt_second) or glob.glob(prefix_model + SECOND_SUFFIX)

            if len(first_files) == 0 or len(second_files) == 0:
                print(f"  Missing files for {sub_id} — skipped.")
                continue

            file_first  = first_files[0]
            file_second = second_files[0]

            try:
                arr_first  = surface.load_surf_data(file_first).astype(float)
                arr_second = surface.load_surf_data(file_second).astype(float)
            except Exception as e:
                print(f"  ERROR loading subject {sub_id}: {e}")
                continue

            if arr_first.shape != arr_second.shape:
                print(f"  Shape mismatch for {sub_id} — skipped.")
                continue

            # ------------------------------------------------
            # Extract ROI values
            # ------------------------------------------------
            roi_first  = arr_first[roi_idx]
            roi_second = arr_second[roi_idx]

            roi_diff = roi_second - roi_first
            roi_dir = (roi_second > roi_first).astype(int)

            # ------------------------------------------------
            # *** Wilcoxon signed-rank test (paired) ***
            # ------------------------------------------------
            # Hypothesis: first > second
            # Wilcoxon tests (first - second) > 0 → alternative="greater"
            try:
                stat, p_wilcoxon = wilcoxon(
                    roi_first, roi_second,
                    alternative="greater",
                    zero_method="wilcox"
                )
            except Exception as e:
                print(f"  Wilcoxon failed for {sub_id}: {e}")
                p_wilcoxon = np.nan
           
            # --------------------------------------------
            # Effect size: Rank-Biserial Correlation (RBC)
            # --------------------------------------------
            diff = roi_first - roi_second
            abs_diff = np.abs(diff)

            # remove zeros (ties) for proper RBC computation
            nonzero = abs_diff > 0
            ranks = rankdata(abs_diff[nonzero])
            signed_ranks = ranks * np.sign(diff[nonzero])

            W_pos = signed_ranks[signed_ranks > 0].sum()
            W_neg = -signed_ranks[signed_ranks < 0].sum()

            if (W_pos + W_neg) == 0:
                rbc = 0.0
            else:
                rbc = (W_pos - W_neg) / (W_pos + W_neg)

            # ------------------------------------------------
            # Store result
            # ------------------------------------------------
            results[(group, hemi, sub_id)] = {
                "roi_first": roi_first,
                "roi_second": roi_second,
                "roi_diff": roi_diff,
                "roi_direction": roi_dir,
                "p_value_wilcoxon": p_wilcoxon,
                'rbc': rbc
            }

            print(f"  {sub_id}: processed — Wilcoxon p = {p_wilcoxon:.4g}")

print("\nAll subjects processed.")
# ----------------------------------------------------
# Compute and print group-level averages
# ----------------------------------------------------
print("\n=== Group-level averages ===")
groups_hemi = set((g, h) for g, h, _ in results.keys())

for group, hemi in groups_hemi:
    group_vals = [vals for (g, h, _), vals in results.items() if g == group and h == hemi]

    p_vals = [v['p_value_wilcoxon'] for v in group_vals if not np.isnan(v['p_value_wilcoxon'])]
    rbcs   = [v['rbc'] for v in group_vals]

    mean_p = np.mean(p_vals) if p_vals else np.nan
    mean_rbc = np.mean(rbcs) if rbcs else np.nan

    print(f"Group {group} | Hemisphere {hemi} → mean Wilcoxon p = {mean_p:.4g}, mean RBC = {mean_rbc:.4g}")
# ----------------------------------------------------
# Save CSV
# ----------------------------------------------------
out_rows = []

for (group, hemi, subject_id), vals in results.items():

    roi_first  = vals["roi_first"]
    roi_second = vals["roi_second"]
    roi_diff   = vals["roi_diff"]
    roi_dir    = vals["roi_direction"]
    p_value    = vals["p_value_wilcoxon"]
    rbc = vals["rbc"]

    out_rows.append({
    "group": group,
    "hemisphere": hemi,
    "subject_id": subject_id,
    "roi_mean_first": np.mean(roi_first),
    "roi_mean_second": np.mean(roi_second),
    "roi_mean_diff": np.mean(roi_diff),
    "roi_prop_second_gt_first": np.mean(roi_dir),
    "roi_N": len(roi_first),
    "wilcoxon_stat": stat,
    "wilcoxon_p": p_value,
    "effect_size_rbc": rbc
})


df_out = pd.DataFrame(out_rows)

output_csv = "/data/p_02915/SPOT/Result/0._results_all_subjects_cross_validation.csv"
df_out.to_csv(output_csv, index=False)

print(f"\nSaved CSV → {output_csv}")
