'''
Statistical comparison of MSD between emperical and simulated across the hemispheres
'''
import pandas as pd
import numpy as np
from scipy.stats import kruskal, wilcoxon, norm
from multiprocessing import Pool, cpu_count
from nilearn import surface
import math
import nibabel as nib


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



def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
     # Ensure labels_area is a list
    if not isinstance(labels_area, list):
        labels_area = [labels_area]
    
    # Collect indices for all labels in labels_area
    indices_area = np.concatenate([
        np.nonzero(visparc.agg_data() == label)[0]
        for label in labels_area
    ])

    return indices_area


def calc_msd(param_area1, param_area2):
    # mean squared difference between two different model results in an area
    n1 = param_area1.size
    n2 = param_area2.size
    assert (
        n1 == n2
    ), "areas don't have the same number of vertices and cannot be compared!"

    return np.sum((param_area1 - param_area2) ** 2) / n1


def bootstrap_resample_mean_diff(n_resamples, combined_data, n_groupA, n_groupB):
    bootstrap_mean_diffs = []
    for _ in range(n_resamples):
        # Resample the combined data
        resample_indices = np.random.choice(
            len(combined_data), size=len(combined_data), replace=True)
        resampled_combined_data = combined_data[resample_indices]

        # Split the resampled data into two groups
        bootstrap_groupA = resampled_combined_data[:n_groupA]
        bootstrap_groupB = resampled_combined_data[n_groupA:]

        # Calculate the mean difference
        mean_diff = np.mean(bootstrap_groupA) - np.mean(bootstrap_groupB)
        bootstrap_mean_diffs.append(mean_diff)
    return bootstrap_mean_diffs

# Function to convert elements to float, setting non-convertible elements to np.nan


def to_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.nan


def flatten(arr):
    return arr.flatten()


LABELS_V1 = [1]
LABELS_V2 = [2, 3]

param = ["polarangle"]
hemispheres = ["L", "R"]

    
for param_index, param in enumerate(param):
    print(param)
    for hemi in hemispheres:
        for area in ["v", "d"]:
            i = 0
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
            df = pd.DataFrame()
            combat_harmonized = pd.read_csv(f"/data/p_02915/SPOT/raw_hemi-{hemi}_area-V2{area}_{param}_S.csv", index_col=None).to_numpy()
            raw_harmonized = pd.read_csv(f"/data/p_02915/SPOT/raw_hemi-{hemi}_area-V3{area}_{param}_S.csv", index_col=None).to_numpy()
            for model_idx, model in enumerate(["2nd", "3rd", "preterm", "fullterm", "adolescent", "adult"]):
                parameters=[]
                parameters2=[]       
                covars = pd.read_csv(f"/data/p_02915/SPOT/covars_hemi-L.csv")
                parameters = np.nanmean(combat_harmonized[:, covars[covars["group"] == model].index], axis=0)
                parameters2 = np.nanmean(raw_harmonized[:, covars[covars["group"] == model].index], axis=0)
                print(model)
                real_mean.append(np.nanmean(parameters))
                simul_mean.append(np.nanmean(parameters2))
                # Replace NaN values with 0
                parameters = np.nan_to_num(parameters, nan=0)
                parameters2 = np.nan_to_num(parameters2, nan=0)

                # Perform Mann-Whitney U tests on each pair of groups
                wilcoxon_results = []
                groupA = np.array(parameters)
                groupB = np.array(parameters2)
                combined_data = np.concatenate([groupA, groupB])
                n_groupA = len(groupA)
                n_groupB = len(groupB)

                res = wilcoxon(groupB, groupA, alternative='two-sided',
                            method="approx")
                print(res.zstatistic)
                mw_statistic = res.statistic
                mw_p_value = res.pvalue
                z_score = res.zstatistic
                wilcoxon_results.append(
                    ("V2", "V3", mw_statistic, mw_p_value))
                print(
                    f"V2 vs. V3 - U statistic: {mw_statistic}, p-value: {mw_p_value}")
                # z_score = mann_whitney_u_to_z(mw_statistic,n_groupA, n_groupB)
                eff_r = z_score / np.sqrt(n_groupA)
                index_label.append(f"{model}")
                original_H_value.append(mw_statistic)
                original_eff.append(eff_r)
                original_z_score.append(z_score)
                original_p_value.append(mw_p_value)
                n_bootstraps = 10000
                n_cores = cpu_count()
                # print(n_cores)
                n_bootstraps_per_core = n_bootstraps // n_cores
                # print(n_bootstraps_per_core)

                with Pool(n_cores) as pool:
                    # Split the work among the cores
                    results = pool.starmap(bootstrap_resample_mean_diff, [(
                        n_bootstraps_per_core, combined_data, n_groupA, n_groupB)] * n_cores)
                bootstrap_mean_diffs = np.concatenate(results)

                # Calculate observed mean difference
                observed_mean_diff = np.mean(combined_data[:n_groupA]) - np.mean(combined_data[n_groupA:])

                # Calculate bootstrap p-value for the mean difference test
                bootstrap_mean_diff_p_value = np.sum(np.abs(bootstrap_mean_diffs) >= np.abs(observed_mean_diff)) / n_bootstraps

                print(f"p_value: {mw_p_value}, bootstrapped p: {bootstrap_mean_diff_p_value}")
                # Calculate bootstrap p-value for Mann-Whitney U test
                
                bootstrap_p_value.append(bootstrap_mean_diff_p_value)
            results_dict = {
                "Original statistics": original_H_value,
                "effect size r": original_eff,
                "original z": original_z_score,
                'Original p-value': original_p_value,
                'Bootstrapped p-value': bootstrap_p_value,
                "V2 mean": real_mean,
                "V3_mean": simul_mean
            }

            # Create a DataFrame from the dictionary
            results_df = pd.DataFrame(results_dict, index=index_label)
            # Save the DataFrame to a CSV file
            results_df.to_csv(f'/data/p_02915/SPOT/01_results_gradient_hemi-{hemi}_area-{area}_{param}_S.csv')
