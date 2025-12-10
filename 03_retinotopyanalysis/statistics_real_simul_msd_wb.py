'''
Statistical comparison of MSD between empirical and simulated across hemispheres
'''

import pandas as pd
import numpy as np
from scipy.stats import kruskal, wilcoxon
from multiprocessing import Pool, cpu_count
import math


def bootstrap_kruskal_resample(n_resamples, group_sizes, combined_data):
    bootstrap_H_statistics = []
    for _ in range(n_resamples):
        resampled_data = np.random.choice(
            combined_data, size=len(combined_data), replace=True)

        resampled_groups = []
        start_idx = 0
        for size in group_sizes:
            resampled_groups.append(resampled_data[start_idx:start_idx + size])
            start_idx += size
 
        bootstrap_H_statistic, _ = kruskal(*resampled_groups)
        bootstrap_H_statistics.append(bootstrap_H_statistic)
    return bootstrap_H_statistics


def bootstrap_resample_wilcoxon(n_resamples, combined_data, n_groupA, n_groupB):
    bootstrap_statistics = []
    for _ in range(n_resamples):
        resample_indices = np.random.choice(
            len(combined_data), size=len(combined_data), replace=True)
        resampled_data = combined_data[resample_indices]

        bootA = resampled_data[:n_groupA]
        bootB = resampled_data[n_groupA:]

        boot_diff = bootB - bootA
        boot_stat, _ = wilcoxon(boot_diff, alternative="two-sided")
        bootstrap_statistics.append(boot_stat)
    return bootstrap_statistics


LABELS_V1 = [1]
LABELS_V2 = [2, 3]

param = ["eccentricity", "polarangle"]

for param_index, param in enumerate(param):
    print(f"\n=== PARAM: {param} ===\n")

    real_mean = []
    simul_mean = []
    original_H_value = []
    original_eff = []
    original_p_value = []
    original_z_score = []
    bootstrap_p_value = []
    bootstrap_ci_lows = []
    bootstrap_ci_highs = []
    index_label = []

    for model in ["fetal<29", "fetal>29", "neonates<37", "neonates>37", "12-16y", "18-21y"]:

        if model =="neonates<37":
            msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_neonates_young.csv")
        elif model == "neonates>37":
            msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_neonates_old.csv")
        elif model == "fetal>29":
            msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_fetal_old.csv")
        elif model == "fetal<29":
            msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_fetal_young.csv")
        elif model == "12-16y":
            msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_HCP_old.csv")
        elif model == "18-21y":
            msd_value = pd.read_csv("/data/p_02915/SPOT/MSD_HCP_young.csv")

        print(model)

        groupA = msd_value[f"L_{param}_real_benson"] + msd_value[f"R_{param}_real_benson"]
        groupB = msd_value[f"L_{param}_real_simulated"] + msd_value[f"R_{param}_real_simulated"]

        groupA = np.array(groupA)
        groupB = np.array(groupB)

        real_mean.append(np.mean(groupA))
        simul_mean.append(np.mean(groupB))

        print("Mean real_benson:", np.mean(groupA))
        print("Mean real_simulated:", np.mean(groupB))

        # KEY CHANGE — use differences directly
        differences = groupB - groupA
        print("Difference = (simulated - benson)")
        print(differences)

        res = wilcoxon(differences, alternative='greater', method="approx")

        mw_statistic = res.statistic
        mw_p_value = res.pvalue
        z_value = res.zstatistic

        # NOW z sign reflects correct direction
        print(f"Wilcoxon signed-rank test: z = {z_value:.4f}, p = {mw_p_value:.5e}")

        effect_r = z_value / np.sqrt(len(differences))

        original_H_value.append(mw_statistic)
        original_z_score.append(z_value)
        original_p_value.append(mw_p_value)
        original_eff.append(effect_r)
        index_label.append(model)

        combined_data = np.concatenate([groupA, groupB])
        n_groupA = len(groupA)
        n_groupB = len(groupB)
        n_bootstraps = 10000
        n_cores = cpu_count()
        n_boot_core = n_bootstraps // n_cores

        with Pool(n_cores) as pool:
            results = pool.starmap(
                bootstrap_resample_wilcoxon,
                [(n_boot_core, combined_data, n_groupA, n_groupB)] * n_cores
            )

        bootstrap_stats = np.concatenate(results)
        boot_p = np.sum(bootstrap_stats >= mw_statistic) / n_bootstraps
        boot_low = np.percentile(bootstrap_stats, 2.5)
        boot_high = np.percentile(bootstrap_stats, 97.5)

        bootstrap_p_value.append(boot_p)
        bootstrap_ci_lows.append(boot_low)
        bootstrap_ci_highs.append(boot_high)

        print(f"Bootstrap p-value: {boot_p}")
        print(f"CI: [{boot_low}, {boot_high}]")

    results = {
        "Original statistics": original_H_value,
        "effect size r": original_eff,
        "original z": original_z_score,
        "Original p-value": original_p_value,
        "Bootstrapped p-value": bootstrap_p_value,
        "Bootstrap CI low": bootstrap_ci_lows,
        "Bootstrap CI high": bootstrap_ci_highs,
        "real_mean": real_mean,
        "simul_mean": simul_mean
    }

    results_df = pd.DataFrame(results, index=index_label)
    results_df.to_csv(f'/data/p_02915/SPOT/01_results_MSD_wb_{param}.csv')

    print("\nSAVED:", f"/data/p_02915/SPOT/01_results_MSD_wb_{param}.csv\n")
