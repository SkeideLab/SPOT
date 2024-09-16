"""
Statistical comparison of ccf results for six groups
"""
import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
from multiprocessing import Pool, cpu_count
import math


def mann_whitney_u_to_z(U, n1, n2):
    # Calculate mean (μ_U)
    mu_U = (n1 * n2) / 2

    # Calculate standard deviation (σ_U)
    sigma_U = math.sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)

    # Calculate z-score
    z = (U - mu_U) / sigma_U

    return z


def bootstrap_kruskal_resample(n_resamples, group_sizes, combined_data):
    bootstrap_H_statistics = []
    for _ in range(n_resamples):
        # Resample the combined data
        resampled_indices = np.random.choice(
            len(combined_data), size=len(combined_data), replace=True)
        resampled_data = combined_data[resampled_indices]

        # Split the resampled data back into groups of the original sizes
        resampled_groups = []
        start_idx = 0
        for size in group_sizes:
            resampled_groups.append(resampled_data[start_idx:start_idx + size])
            start_idx += size

        # Compute the Kruskal-Wallis statistic for the resampled groups
        bootstrap_H_statistic, _ = kruskal(*resampled_groups)
        bootstrap_H_statistics.append(bootstrap_H_statistic)
    return bootstrap_H_statistics


def bootstrap_resample_mannwhitney(n_resamples, combined_data, n_groupA, n_groupB):
    bootstrap_mannwhitney_statistics = []
    for _ in range(n_resamples):
        # Resample the combined data
        resample_indices = np.random.choice(
            len(combined_data), size=len(combined_data), replace=True)
        resampled_combined_data = combined_data[resample_indices]

        # Split the resampled data into two groups
        bootstrap_groupA = flatten(resampled_combined_data[:n_groupA])
        bootstrap_groupB = flatten(resampled_combined_data[n_groupA:])

        # Perform Mann-Whitney U test on the bootstrap sample
        bootstrap_mw_statistic, _ = mannwhitneyu(
            bootstrap_groupA, bootstrap_groupB, alternative="two-sided")
        bootstrap_mannwhitney_statistics.append(bootstrap_mw_statistic)
    return bootstrap_mannwhitney_statistics


def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        np.logical_or(
            visparc.agg_data() == labels_area[0], visparc.agg_data(
            ) == labels_area[1]
        )
    )[0]
    return indices_area


def flatten(arr):
    return arr.flatten()


VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")
LABELS_V2 = (3, 4)

print(bin)
df = pd.DataFrame()
for param in ["_dens-164k_desc-real_r", "_dens-164k_desc-real_sigma",
              "_dens-164k_label-eccentricity_desc-real_roi-v2th00_metric",
              "_dens-164k_label-polarangle_desc-real_roi-v2th00_metric"]:
    results_list = []
    print(param)
    for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]:
        if group == "neonates<37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage{param}.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "neonates>37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage{param}.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "fetal<29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage{param}.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "fetal>29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage{param}.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "12-16y":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
                "{sub}_hemi-{hemi}_mesh-fsaverage{param}.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "18-21y":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
                "{sub}_hemi-{hemi}_mesh-fsaverage{param}.gii"
            )
            subject_info = pd.read_csv(
                '/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
            sub_num = len(subject_info["sub_id"])

        for hemi in ["L", "R"]:
            visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
            indices_v2 = get_indices_roi(LABELS_V2, visparc)
            group_name = f"{group}_{hemi}"
            parameters = []
        # FORMAT PATHS FOR INPUT AND OUTPUT
            for index, row in subject_info.iterrows():
                if group == "12-16y" or group == "18-21y":
                    sub_id = subject_info.at[index, "sub_id"]
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub_id,
                        hemi=hemi,
                        param=param,
                    )
                    input_path = prefix_model

                    ccf = surface.load_surf_data(input_path)
                    ccf_v0 = ccf[indices_v2].astype(np.float64)
                    parameters.append(ccf_v0)
                else:
                    sub_id = subject_info.at[index, "sub_id"]
                    sess_id = subject_info.at[index, "sess_id"]
                    sub = sub_id.replace('sub-', '')
                    ses = sess_id.replace('ses-', '')
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                        param=param,
                    )
                    input_path = prefix_model

                    ccf = surface.load_surf_data(input_path)
                    ccf_v0 = ccf[indices_v2].astype(np.float64)
                    parameters.append(ccf_v0)

            if group == "neonates<37" and hemi == "L":
                neonates_y_L = np.array(parameters)
            elif group == "neonates>37" and hemi == "L":
                neonates_o_L = np.array(parameters)
            elif group == "fetal<29" and hemi == "L":
                fetal_y_L = np.array(parameters)
            elif group == "fetal>29" and hemi == "L":
                fetal_o_L = np.array(parameters)
            elif group == "12-16y" and hemi == "L":
                hcp_y_L = np.array(parameters)
            elif group == "18-21y" and hemi == "L":
                hcp_o_L = np.array(parameters)
            elif group == "neonates<37" and hemi == "R":
                neonates_y_R = np.array(parameters)
            elif group == "neonates>37" and hemi == "R":
                neonates_o_R = np.array(parameters)
            elif group == "fetal<29" and hemi == "R":
                fetal_y_R = np.array(parameters)
            elif group == "fetal>29" and hemi == "R":
                fetal_o_R = np.array(parameters)
            elif group == "12-16y" and hemi == "R":
                hcp_y_R = np.array(parameters)
            elif group == "18-21y" and hemi == "R":
                hcp_o_R = np.array(parameters)

            data = np.array(parameters).reshape(-1)
            data_nozero_flatten = data.flatten()
            data_r = pd.DataFrame(
                {'r': data_nozero_flatten, 'Group': group_name})
            df = pd.concat([df, data_r], ignore_index=True)
            med = np.max(data_nozero_flatten)

    n_bootstraps = 10000

    for hemi in ["Left", "Right"]:
        if hemi == "Left":
            groups_flatt = {
                # [fetal_y_L.flatten() != 0],
                'fetal_y_L': fetal_y_L.flatten(),
                # [fetal_o_L.flatten() != 0],
                'fetal_o_L': fetal_o_L.flatten(),
                # [neonates_y_L.flatten() != 0],
                'neonates_y_L': neonates_y_L.flatten(),
                # [neonates_o_L.flatten() != 0],
                'neonates_o_L': neonates_o_L.flatten(),
                'hcp_y_L': hcp_y_L.flatten(),  # [hcp_y_L.flatten() != 0],
                'hcp_o_L': hcp_o_L.flatten(),  # [hcp_o_L.flatten() != 0]
            }
            groups = {
                'fetal_y_L': fetal_y_L,  # [fetal_y_L != 0],
                'fetal_o_L': fetal_o_L,  # [fetal_o_L != 0],
                'neonates_y_L': neonates_y_L,  # [neonates_y_L != 0],
                'neonates_o_L': neonates_o_L,  # [neonates_o_L != 0],
                'hcp_y_L': hcp_y_L,  # [hcp_y_L != 0],
                'hcp_o_L': hcp_o_L,  # [hcp_o_L != 0]
            }
        elif hemi == "Right":

            groups_flatt = {
                # [fetal_y_R.flatten() != 0],
                'fetal_y_R': fetal_y_R.flatten(),
                # [fetal_o_R.flatten() != 0],
                'fetal_o_R': fetal_o_R.flatten(),
                # [neonates_y_R.flatten() != 0],
                'neonates_y_R': neonates_y_R.flatten(),
                # [neonates_o_R.flatten() != 0],
                'neonates_o_R': neonates_o_R.flatten(),
                'hcp_y_R': hcp_y_R.flatten(),  # [hcp_y_R.flatten() != 0],
                'hcp_o_R': hcp_o_R.flatten(),  # [hcp_o_R.flatten() != 0]
            }
            groups = {
                'fetal_y_R': fetal_y_R,  # [fetal_y_R != 0],
                'fetal_o_R': fetal_o_R,  # [fetal_o_R != 0],
                'neonates_y_R': neonates_y_R,  # [neonates_y_R != 0],
                'neonates_o_R': neonates_o_R,  # [neonates_o_R != 0],
                'hcp_y_R': hcp_y_R,  # [hcp_y_R != 0],
                'hcp_o_R': hcp_o_R,  # [hcp_o_R != 0]
            }

        all_groups_flatt = list(groups_flatt.values())
        i = 0
        original_H_value = []
        original_z_score = []
        original_p_value = []
        bootstrap_p_value = []
        bootstrap_ci_lows = []
        bootstrap_ci_highs = []
        index_label = []
        # Perform the original Kruskal-Wallis test
        kruskal_statistic, kruskal_p_value = kruskal(*all_groups_flatt)
        print(hemi)
        print(
            f"Original Kruskal-Wallis statistic: {kruskal_statistic}, p-value: {kruskal_p_value}")
        k = 6  # Number of groups
        # Total number of observations
        N = N = sum(len(group) for group in all_groups_flatt)

        eta_squared = np.sqrt((kruskal_statistic - k + 1) / (N - k))
        original_H_value.append(kruskal_statistic)
        original_p_value.append(kruskal_p_value)
        original_z_score.append(eta_squared)
        print(f"w = {eta_squared}")
        group_names = list(groups.keys())
        all_groups = [groups[name] for name in group_names]
        group_sizes = [len(group) for group in all_groups]

        # Combine all data for bootstrapping
        combined_data = np.concatenate(all_groups)

        n_cores = cpu_count()
        # print(n_cores)
        n_bootstraps_per_core = n_bootstraps // n_cores

        with Pool(n_cores) as pool:
            results = pool.starmap(bootstrap_kruskal_resample, [(
                n_bootstraps_per_core, group_sizes, combined_data)] * n_cores)

        # Combine results from all cores
        bootstrap_kruskal_statistics = np.concatenate(results)

        # Calculate the bootstrap p-value for Kruskal-Wallis test
        bootstrap_kruskal_statistics = np.array(bootstrap_kruskal_statistics)
        bootstrap_kruskal_p_value = np.sum(
            bootstrap_kruskal_statistics >= kruskal_statistic) / n_bootstraps

        print("\nBootstrapped Kruskal-Wallis H-statistic:", kruskal_statistic)
        print("Bootstrapped p-value:", bootstrap_kruskal_p_value)
        bootstrap_ci_low = np.percentile(bootstrap_kruskal_statistics, 2.5)
        bootstrap_ci_high = np.percentile(bootstrap_kruskal_statistics, 97.5)

        print(
            f"Bootstrap CI for Kruskal-Wallis H statistic - Low: {bootstrap_ci_low}, High: {bootstrap_ci_high}")

        bootstrap_p_value.append(bootstrap_kruskal_p_value)
        bootstrap_ci_lows.append(bootstrap_ci_low)
        bootstrap_ci_highs.append(bootstrap_ci_high)
        index_label.append("6-class")

        # Generate all pairs of groups
        from itertools import combinations
        group_pairs = list(combinations(groups.keys(), 2))

        # Perform Mann-Whitney U tests on each pair of groups
        mannwhitney_results = []

        for (groupA_name, groupB_name) in group_pairs:
            i = i+1
            groupA = groups[groupA_name]
            groupB = groups[groupB_name]
            groupA_flatt = groups_flatt[groupA_name]
            groupB_flatt = groups_flatt[groupB_name]

            # Combine data for Kruskal-Wallis test
            combined_data = np.concatenate([groupA, groupB])
            n_groupA = len(groupA)
            n_groupB = len(groupB)

            mw_statistic, mw_p_value = mannwhitneyu(
                groupA_flatt, groupB_flatt, alternative='two-sided')
            mannwhitney_results.append(
                (groupA_name, groupB_name, mw_statistic, mw_p_value))
            print(
                f"Original Mann-Whitney U test {groupA_name} vs. {groupB_name} - U statistic: {mw_statistic}, p-value: {mw_p_value}")
            z_score = mann_whitney_u_to_z(
                mw_statistic, len(groupA_flatt), len(groupB_flatt))
            index_label.append(f"{groupA_name} vs. {groupB_name}")
            original_H_value.append(mw_statistic)
            original_p_value.append(mw_p_value)
            original_z_score.append(z_score)

            with Pool(n_cores) as pool:
                # Split the work among the cores
                results = pool.starmap(bootstrap_resample_mannwhitney, [(
                    n_bootstraps_per_core, combined_data, n_groupA, n_groupB)] * n_cores)

            # Combine results from all cores
            bootstrap_mannwhitney_statistics = np.concatenate(results)
            # Calculate bootstrap p-value for Mann-Whitney U test
            bootstrap_mw_p_value = np.sum(
                bootstrap_mannwhitney_statistics >= mw_statistic) / n_bootstraps

            print(
                f"Bootstrapped Mann-Whitney U test {groupA_name} vs {groupB_name} - p-value: {bootstrap_mw_p_value}")
            # Compute confidence intervals for the bootstrap distribution
            bootstrap_ci_low = np.percentile(
                bootstrap_mannwhitney_statistics, 2.5)
            bootstrap_ci_high = np.percentile(
                bootstrap_mannwhitney_statistics, 97.5)

            print(
                f"Bootstrap CI for Mann-Whitney U statistic - Low: {bootstrap_ci_low}, High: {bootstrap_ci_high}")

            bootstrap_p_value.append(bootstrap_mw_p_value)
            bootstrap_ci_lows.append(bootstrap_ci_low)
            bootstrap_ci_highs.append(bootstrap_ci_high)

        results_dict = {
            "Original H": original_H_value,
            "original z": original_z_score,
            'Original p-value': original_p_value,
            'Bootstrapped p-value': bootstrap_p_value,
            'Bootstrap CI low': bootstrap_ci_lows,
            'Bootstrap CI high': bootstrap_ci_highs
        }

        # Create a DataFrame from the dictionary
        results_df = pd.DataFrame(results_dict, index=index_label)

        # Save the DataFrame to a CSV file
        results_df.to_csv(f'/data/p_02915/SPOT/01_results_{param}_{hemi}.csv')
