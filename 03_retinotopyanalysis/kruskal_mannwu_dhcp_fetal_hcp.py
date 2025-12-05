"""
Statistical comparison of ccf results for six groups
"""
import numpy as np
import pandas as pd 
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
def flatten(arr):
    return arr.flatten()

for param in ["desc-real_r", "desc-real_sigma","label-eccentricity_desc-real_roi-v2th00_metric","label-polarangle_desc-real_roi-v2th00_metric"]: #
    if param =="desc-real_r":
        test_value = "r"
    elif param == "desc-real_sigma":
        test_value = "sigma"
    elif param == "label-eccentricity_desc-real_roi-v2th00_metric":
        test_value = "eccentricity"
    elif param =="label-polarangle_desc-real_roi-v2th00_metric":
        test_value = "polarangle"

    for area in ["V2_V3"]:
        if area == "V2":
            LABELS_V2 = [2]
        elif area == "V3":
            LABELS_V2 = [3]
        elif area == "V2_V3":
            LABELS_V2 = [2, 3]

        for hemi in ["L", "R"]:
            combat_harmonized = pd.read_csv(f"/data/p_02915/SPOT/Result/combat_hemi-{hemi}_area-{area}_{test_value}.csv", index_col=None).to_numpy()
            covars = pd.read_csv(f"/data/p_02915/SPOT/Result/covars_hemi-{hemi}.csv")
            print(type(combat_harmonized))
            groups_flatt = {               
                        '2nd': combat_harmonized[:, covars[covars["group"] == "2nd"].index].flatten(),
                        '3rd': combat_harmonized[:, covars[covars["group"] == "3rd"].index].flatten(),
                        'preterm': combat_harmonized[:, covars[covars["group"] == "preterm"].index].flatten(),
                        'fullterm': combat_harmonized[:, covars[covars["group"] == "fullterm"].index].flatten(),
                        'adolescent': combat_harmonized[:, covars[covars["group"] == "adolescent"].index].flatten(),
                        'adult': combat_harmonized[:, covars[covars["group"] == "adult"].index].flatten()
                    }
            groups = {
                        '2nd': combat_harmonized[:, covars[covars["group"] == "2nd"].index].T,
                        '3rd': combat_harmonized[:, covars[covars["group"] == "3rd"].index].T,
                        'preterm': combat_harmonized[:, covars[covars["group"] == "preterm"].index].T,
                        'fullterm': combat_harmonized[:, covars[covars["group"] == "fullterm"].index].T,
                        'adolescent': combat_harmonized[:, covars[covars["group"] == "adolescent"].index].T,
                        'adult': combat_harmonized[:, covars[covars["group"] == "adult"].index].T
                    }
            all_groups_flatt = list(groups_flatt.values())

            n_bootstraps = 10000
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
            #bootstrap_p_value.append(0)
            #bootstrap_ci_lows.append(0)
            #bootstrap_ci_highs.append(0)
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
            results_df.to_csv(f'/data/p_02915/SPOT/Result/01_results_hemi-{hemi}_area-{area}_{test_value}.csv')
