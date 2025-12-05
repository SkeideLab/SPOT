'''
Statistical comparison of MSD between emperical and simulated for each hemisphere
'''
import pandas as pd
from itertools import combinations
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
from multiprocessing import Pool, cpu_count
import math


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


def mann_whitney_u_to_z(U, n1, n2):
    # Calculate mean (μ_U)
    mu_U = (n1 * n2) / 2

    # Calculate standard deviation (σ_U)
    sigma_U = math.sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)

    # Calculate z-score
    z = (U - mu_U) / sigma_U

    return z


def flatten(arr):
    return arr.flatten()


original_H_value = []
original_z_score = []
original_p_value = []
bootstrap_p_value = []
bootstrap_ci_lows = []
bootstrap_ci_highs = []
index_label = []

df = pd.DataFrame()
for hemi in ["L", "R"]:
    for model in ["eccentricity", "polarangle"]:
        msd_real_value = []
        msd_simul_value = []
        for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]:
            if group == "neonates<37":
                MSD = pd.read_csv("/data/p_02915/SPOT/MSD_neonates_young.csv")
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37_v2.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "neonates>37":
                MSD = pd.read_csv("/data/p_02915/SPOT/MSD_neonates_old.csv")
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37_v2.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "fetal<29":
                MSD = pd.read_csv("/data/p_02915/SPOT/MSD_fetal_young.csv")
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "fetal>29":
                MSD = pd.read_csv("/data/p_02915/SPOT/MSD_fetal_old.csv")
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_old.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "12-16y":
                MSD = pd.read_csv("/data/p_02915/SPOT/MSD_HCP_young.csv")
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/hcp_subj_path_SPOT_young.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "18-21y":
                MSD = pd.read_csv("/data/p_02915/SPOT/MSD_HCP_old.csv")
                subject_info = pd.read_csv(
                    '/data/p_02915/SPOT/hcp_subj_path_SPOT_old.csv')
                sub_num = len(subject_info["sub_id"])

            msd_real = MSD[f"{hemi}_{model}_real_benson"]
            msd_simul = MSD[f"{hemi}_{model}_real_simulated"]
            msd_real_value.extend(msd_real)
            msd_simul_value.extend(msd_simul)

            if group == "neonates<37":
                neonates_y = msd_real-msd_simul
            elif group == "neonates>37":
                neonates_o = msd_real-msd_simul
            elif group == "fetal<29":
                fetal_y = msd_real-msd_simul
            elif group == "fetal>29":
                fetal_o = msd_real-msd_simul
            elif group == "12-16y":
                hcp_y = msd_real-msd_simul
            elif group == "18-21y":
                hcp_o = msd_real-msd_simul

        med = np.mean(msd_real_value)
        med1 = np.mean(msd_simul_value)
        print(f"{model} on {hemi} hemi real msd: {med}")
        print(f"{model} on {hemi} hemi simul msd: {med1}")

        n_bootstraps = 10000
        groups = {
            'fetal_y': fetal_y,  # [fetal_y_L != 0],
            'fetal_o': fetal_o,  # [fetal_o_L != 0],
            'neonates_y': neonates_y,  # [neonates_y_L != 0],
            'neonates_o': neonates_o,  # [neonates_o_L != 0],
            'hcp_y': hcp_y,  # [hcp_y_L != 0],
            'hcp_o': hcp_o,  # [hcp_o_L != 0]
        }

        all_groups_flatt = list(groups.values())
        i = 0
        original_H_value = []
        original_p_value = []
        original_z_score = []
        bootstrap_p_value = []
        bootstrap_ci_lows = []
        bootstrap_ci_highs = []
        index_label = []
        # Perform the original Kruskal-Wallis test
        kruskal_statistic, kruskal_p_value = kruskal(*all_groups_flatt)
        print(hemi)
        print(
            f"Original Kruskal-Wallis statistic: {kruskal_statistic}, p-value: {kruskal_p_value}")
        original_H_value.append(kruskal_statistic)
        original_z_score.append(kruskal_statistic)
        original_p_value.append(kruskal_p_value)
        group_names = list(groups.keys())
        all_groups = [groups[name] for name in group_names]
        group_sizes = [len(group) for group in all_groups]

        # Combine all data for bootstrapping
        combined_data = np.concatenate(all_groups)

        n_cores = cpu_count()
        # print(n_cores)
        n_bootstraps_per_core = n_bootstraps // n_cores
        # print(n_bootstraps_per_core)

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

        group_pairs = list(combinations(groups.keys(), 2))

        # Perform Mann-Whitney U tests on each pair of groups
        mannwhitney_results = []

        for (groupA_name, groupB_name) in group_pairs:
            i = i+1
            groupA = groups[groupA_name]
            groupB = groups[groupB_name]
            # Combine data for Kruskal-Wallis test
            combined_data = np.concatenate([groupA, groupB])
            n_groupA = len(groupA)
            n_groupB = len(groupB)

            mw_statistic, mw_p_value = mannwhitneyu(
                groupA, groupB, alternative='two-sided')
            mannwhitney_results.append(
                (groupA_name, groupB_name, mw_statistic, mw_p_value))
            print(
                f"Original Mann-Whitney U test {groupA_name} vs. {groupB_name} - U statistic: {mw_statistic}, p-value: {mw_p_value}")
            z_score = mann_whitney_u_to_z(mw_statistic, n_groupA, n_groupB)
            index_label.append(f"{groupA_name} vs. {groupB_name}")
            original_H_value.append(mw_statistic)
            original_z_score.append(z_score)
            original_p_value.append(mw_p_value)

            with Pool(n_cores) as pool:
                # Split the work among the cores
                results = pool.starmap(bootstrap_resample_mannwhitney, [(
                    n_bootstraps_per_core, combined_data, n_groupA, n_groupB)] * n_cores)

            # Combine results from all cores
            bootstrap_mannwhitney_statistics = np.concatenate(results)
            # Calculate bootstrap p-value for Mann-Whitney U test
            bootstrap_mw_p_value = np.sum(
                bootstrap_mannwhitney_statistics >= mw_statistic) / n_bootstraps
            if bootstrap_mw_p_value == 1:
                print("Sample of bootstrap statistics:",
                      bootstrap_mannwhitney_statistics[:10])
                # Check if all bootstrap statistics are >= the observed statistic
                comparison = bootstrap_mannwhitney_statistics >= mw_statistic
                print("Comparison array sample:", comparison[:20])
                print("Sum of comparison array:", np.sum(comparison))

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
        results_df.to_csv(f'/data/p_02915/SPOT/01_results_MSD_{model}_{hemi}.csv')
