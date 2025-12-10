import pandas as pd
import numpy as np
import math
from itertools import combinations
from scipy.stats import kruskal, mannwhitneyu
from multiprocessing import Pool, cpu_count

def bootstrap_kruskal_resample(n_resamples, group_sizes, combined_data):
    """Performs bootstrap resampling for Kruskal-Wallis test."""
    bootstrap_H_statistics = []
    
    for _ in range(n_resamples):
        resampled_indices = np.random.choice(len(combined_data), size=len(combined_data), replace=True)
        resampled_data = combined_data[resampled_indices]

        # Split resampled data into groups
        resampled_groups = np.split(resampled_data, np.cumsum(group_sizes)[:-1])
        
        bootstrap_H_statistic, _ = kruskal(*resampled_groups)
        bootstrap_H_statistics.append(bootstrap_H_statistic)
    
    return bootstrap_H_statistics
 
def bootstrap_resample_mannwhitney(n_resamples, combined_data, n_groupA):
    """Performs bootstrap resampling for Mann-Whitney U test."""
    bootstrap_mannwhitney_statistics = []
    
    for _ in range(n_resamples):
        resampled_indices = np.random.choice(len(combined_data), size=len(combined_data), replace=True)
        resampled_data = combined_data[resampled_indices]
        
        bootstrap_groupA, bootstrap_groupB = resampled_data[:n_groupA], resampled_data[n_groupA:]
        bootstrap_mw_statistic, _ = mannwhitneyu(bootstrap_groupA, bootstrap_groupB, alternative="two-sided")
        
        bootstrap_mannwhitney_statistics.append(bootstrap_mw_statistic)
    
    return bootstrap_mannwhitney_statistics

def mann_whitney_u_to_z(U, n1, n2):
    """Converts Mann-Whitney U statistic to z-score."""
    mu_U = (n1 * n2) / 2
    sigma_U = math.sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
    
    return (U - mu_U) / sigma_U

def analyze_snr_data(param):
    results_list = []
    print(param)
    """Analyzes SNR data using Kruskal-Wallis and Mann-Whitney U tests."""
    combat_harmonized = pd.read_csv(f"/data/p_02915/SPOT/Result/raw_{param}snr.csv").to_numpy()
    covars = pd.read_csv("/data/p_02915/SPOT/Result/covars_hemi-L.csv")

    groups = {
        group: combat_harmonized[covars.loc[covars["group"] == group].index]  
        for group in ["2nd", "3rd", "preterm", "fullterm", "adolescent", "adult"]
    }


    group_names = list(groups.keys())
    all_groups = list(groups.values())
    group_sizes = [len(group) for group in all_groups]
    
    # Perform Kruskal-Wallis Test
    kruskal_statistic, kruskal_p_value = kruskal(*all_groups)
    print(f"Kruskal-Wallis H: {kruskal_statistic}, p-value: {kruskal_p_value}")

    # Compute effect size
    k = len(group_names)
    N = sum(group_sizes)
    eta_squared = np.sqrt((kruskal_statistic - k + 1) / (N - k))
    
    # Bootstrapping
    n_bootstraps = 10000
    combined_data = np.concatenate(all_groups)

    with Pool(cpu_count()) as pool:
        results = pool.starmap(bootstrap_kruskal_resample, [
            (n_bootstraps // cpu_count(), group_sizes, combined_data)
        ] * cpu_count())

    bootstrap_H_values = np.concatenate(results)
    bootstrap_p_value = np.sum(bootstrap_H_values >= kruskal_statistic) / n_bootstraps
    bootstrap_ci_low, bootstrap_ci_high = np.percentile(bootstrap_H_values, [2.5, 97.5])

    print(f"Bootstrapped H: {kruskal_statistic}, p-value: {bootstrap_p_value}")
    print(f"Bootstrap CI: {bootstrap_ci_low} - {bootstrap_ci_high}")
    results_list.append({
            "Comparison": "6-class",
            "Original U": kruskal_statistic[0],
            "Original p-value": kruskal_p_value[0],
            "Original Z": eta_squared,
            "Bootstrap p-value": bootstrap_p_value,
            "Bootstrap CI Low": bootstrap_ci_low,
            "Bootstrap CI High": bootstrap_ci_high
        })

    # Mann-Whitney U Tests
    for groupA, groupB in combinations(group_names, 2):
        groupA_data, groupB_data = groups[groupA], groups[groupB]
        combined_data = np.concatenate([groupA_data, groupB_data])
        n_groupA = len(groupA_data)

        mw_statistic, mw_p_value = mannwhitneyu(groupA_data, groupB_data, alternative='two-sided')
        z_score = mann_whitney_u_to_z(mw_statistic, len(groupA_data), len(groupB_data))

        # Bootstrap Mann-Whitney
        with Pool(cpu_count()) as pool:
            results = pool.starmap(bootstrap_resample_mannwhitney, [
                (n_bootstraps // cpu_count(), combined_data, n_groupA)
            ] * cpu_count())

        bootstrap_mw_values = np.concatenate(results)
        bootstrap_mw_p_value = np.sum(bootstrap_mw_values >= mw_statistic) / n_bootstraps
        bootstrap_ci_low, bootstrap_ci_high = np.percentile(bootstrap_mw_values, [2.5, 97.5])

        print(f"Mann-Whitney {groupA} vs {groupB}: U={mw_statistic}, p={mw_p_value}, Z={z_score}")
        print(f"Bootstrap MW p-value: {bootstrap_mw_p_value}, CI: {bootstrap_ci_low} - {bootstrap_ci_high}")

        results_list.append({
            "Comparison": f"{groupA} vs {groupB}",
            "Original U": mw_statistic[0],
            "Original p-value": mw_p_value[0],
            "Original Z": z_score[0],
            "Bootstrap p-value": bootstrap_mw_p_value,
            "Bootstrap CI Low": bootstrap_ci_low,
            "Bootstrap CI High": bootstrap_ci_high
        })

    # Save results to CSV
    results_df = pd.DataFrame(results_list)
    results_df.to_csv(f'/data/p_02915/SPOT/Result/01_results_{param}snr_raw_v2.csv', index=False)

if __name__ == "__main__":
    for param in ["t", "s"]:
        analyze_snr_data(param)
