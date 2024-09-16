"""
Statistical tests of DVARS
"""
import pandas as pd
from itertools import combinations
import numpy as np
from scipy.stats import kruskal, mannwhitneyu, normaltest
from multiprocessing import Pool, cpu_count
import math
def bootstrap_kruskal_resample(n_resamples, group_sizes, combined_data):
    bootstrap_H_statistics = []
    for _ in range(n_resamples):
        # Resample the combined data
        resampled_indices = np.random.choice(len(combined_data), size=len(combined_data), replace=True)
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
        resample_indices = np.random.choice(len(combined_data), size=len(combined_data), replace=True)
        resampled_combined_data = combined_data[resample_indices]
        
        # Split the resampled data into two groups
        bootstrap_groupA = flatten(resampled_combined_data[:n_groupA])
        bootstrap_groupB = flatten(resampled_combined_data[n_groupA:])
        
        # Perform Mann-Whitney U test on the bootstrap sample
        bootstrap_mw_statistic,_ = mannwhitneyu(bootstrap_groupA, bootstrap_groupB, alternative="two-sided")
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

# Function to convert elements to float, setting non-convertible elements to np.nan
def to_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.nan
    
def flatten(arr):
    return arr.flatten()
 
original_H_value=[]
original_z_score=[]
original_p_value=[]
bootstrap_p_value=[]
bootstrap_ci_lows=[]
bootstrap_ci_highs=[]
index_label = []
i=0   
df = pd.DataFrame()
for hemi in ["L", "R"]:        
    for group in ["12-16y", "18-21y", "neonates<37", "neonates>37", "fetal<29", "fetal>29"]: 
        results_list = []    
        if group == "neonates<37":
            PREFIX_MODEL = ('/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_dvars.csv')
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
            sub_num = len(subject_info["sub_id"])   
        elif group == "neonates>37":
            PREFIX_MODEL = ('/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_dvars.csv')
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "fetal<29":
            PREFIX_MODEL = ('/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_dvars.csv')
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
            sub_num = len(subject_info["sub_id"])        
        elif group == "fetal>29":
            PREFIX_MODEL = ('/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{ses}_{hemi}_dvars.csv')
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "12-16y":
            PREFIX_MODEL = ('/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{hemi}_dvars.csv')
            subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "18-21y":
            PREFIX_MODEL = ('/data/p_02915/dhcp_derivatives_SPOT/statistics/{sub}_{hemi}_dvars.csv')
            subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
            sub_num = len(subject_info["sub_id"])
            
        group_name=f"{group}"
        parameters=[]
        parameters2=np.array([[]])
    # FORMAT PATHS FOR INPUT AND OUTPUT      
    
        for index, row in subject_info.iterrows():                        
            if group == "12-16y" or group == "18-21y":
                sub_id = subject_info.at[index, "sub_id"]
                prefix_model = PREFIX_MODEL.format(
                    sub=sub_id,
                    hemi=hemi,
                )
                input_path = prefix_model 
            else:
                sub_id = subject_info.at[index, "sub_id"]
                sess_id = subject_info.at[index, "sess_id"]
                sub = sub_id.replace('sub-','')
                ses = sess_id.replace('ses-','')
                prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                input_path = prefix_model

            both_data= pd.read_csv(prefix_model)
            parameters.append(np.mean(both_data))
            if parameters2.shape[1] == both_data.shape[1] or parameters2.shape[1] == 0:
                # If parameters2 is empty (shape (1, 0)), initialize it correctly
                if parameters2.shape[1] == 0:
                    parameters2 = both_data
                else:
                    parameters2 = np.vstack((parameters2, both_data))
        
        print(group_name)   
        print(np.mean(parameters2))
        print(np.std(parameters2))
        print(np.max(parameters2))

        if group == "neonates<37":
            neonates_y=np.array(parameters)
        elif group == "neonates>37":
            neonates_o=np.array(parameters)
        elif group == "fetal<29":
            fetal_y=np.array(parameters)
        elif group == "fetal>29":
            fetal_o=np.array(parameters)
        elif group == "12-16y":
            hcp_y=np.array(parameters)
        elif group == "18-21y":
            hcp_o=np.array(parameters)

    n_bootstraps=10000
    groups_flatt = {
    'fetal_y': fetal_y,#[fetal_y_L != 0],
    'fetal_o': fetal_o,#[fetal_o_L != 0],
    'neonates_y': neonates_y,#[neonates_y_L != 0],
    'neonates_o': neonates_o,#[neonates_o_L != 0],
    'hcp_y': hcp_y,#[hcp_y_L != 0],
    'hcp_o': hcp_o,#[hcp_o_L != 0]
    }

    all_groups_flatt = list(groups_flatt.values())      
    # Perform the original Kruskal-Wallis test
    kruskal_statistic, kruskal_p_value = kruskal(*all_groups_flatt)
    print(hemi)
    print(f"Original Kruskal-Wallis statistic: {kruskal_statistic}, p-value: {kruskal_p_value}")
    k = 6  # Number of groups
    N = N = sum(len(group) for group in all_groups_flatt)   # Total number of observations
    print(N)
    eta_squared = np.sqrt((kruskal_statistic - k + 1) / N -k)
    original_H_value.append(kruskal_statistic)
    original_z_score.append(eta_squared)
    original_p_value.append(kruskal_p_value)
    group_names = list(groups_flatt.keys())
    all_groups = [groups_flatt[name] for name in group_names]
    group_sizes = [len(group) for group in all_groups]

    # Combine all data for bootstrapping
    combined_data = np.concatenate(all_groups)

    n_cores = cpu_count()
    #print(n_cores)
    n_bootstraps_per_core=n_bootstraps // n_cores
    #print(n_bootstraps_per_core)

    with Pool(n_cores) as pool:
        results = pool.starmap(bootstrap_kruskal_resample, [(n_bootstraps_per_core, group_sizes, combined_data)] * n_cores)

    # Combine results from all cores
    bootstrap_kruskal_statistics = np.concatenate(results)

    # Calculate the bootstrap p-value for Kruskal-Wallis test
    bootstrap_kruskal_statistics = np.array(bootstrap_kruskal_statistics)
    bootstrap_kruskal_p_value = np.sum(bootstrap_kruskal_statistics >= kruskal_statistic) / n_bootstraps

    print("\nBootstrapped Kruskal-Wallis H-statistic:", kruskal_statistic)
    print("Bootstrapped p-value:", bootstrap_kruskal_p_value)
    bootstrap_ci_low = np.percentile(bootstrap_kruskal_statistics, 2.5)
    bootstrap_ci_high = np.percentile(bootstrap_kruskal_statistics, 97.5)

    print(f"Bootstrap CI for Kruskal-Wallis H statistic - Low: {bootstrap_ci_low}, High: {bootstrap_ci_high}")

    bootstrap_p_value.append(bootstrap_kruskal_p_value)
    bootstrap_ci_lows.append(bootstrap_ci_low)
    bootstrap_ci_highs.append(bootstrap_ci_high)
    index_label.append(f"{hemi}_6-class")
    i=i+1

    # Generate all pairs of groups

    group_pairs = list(combinations(groups_flatt.keys(), 2))

    # Perform Mann-Whitney U tests on each pair of groups
    mannwhitney_results = []

    for (groupA_name, groupB_name) in group_pairs:
        
        groupA = groups_flatt[groupA_name]
        groupB = groups_flatt[groupB_name]
        # Combine data for Kruskal-Wallis test
        combined_data = np.concatenate([groupA, groupB])
        n_groupA = len(groupA)
        n_groupB = len(groupB)

        mw_statistic, mw_p_value = mannwhitneyu(groupA, groupB, alternative='two-sided')
        mannwhitney_results.append((groupA_name, groupB_name, mw_statistic, mw_p_value))
        print(f"Original Mann-Whitney U test {groupA_name} vs. {groupB_name} - U statistic: {mw_statistic}, p-value: {mw_p_value}")
        z_score=mann_whitney_u_to_z(mw_statistic,n_groupA, n_groupB)
        index_label.append(f"{groupA_name} vs. {groupB_name}")
        original_H_value.append(mw_statistic)
        original_z_score.append(z_score)
        original_p_value.append(mw_p_value)

        with Pool(n_cores) as pool:
        # Split the work among the cores
            results = pool.starmap(bootstrap_resample_mannwhitney, [(n_bootstraps_per_core, combined_data, n_groupA, n_groupB)] * n_cores)

        # Combine results from all cores
        bootstrap_mannwhitney_statistics = np.concatenate(results)            
        # Calculate bootstrap p-value for Mann-Whitney U test
        bootstrap_mw_p_value = np.sum(bootstrap_mannwhitney_statistics >= mw_statistic) / n_bootstraps
        if bootstrap_mw_p_value == 1: 
            print("Sample of bootstrap statistics:", bootstrap_mannwhitney_statistics[:10])
            # Check if all bootstrap statistics are >= the observed statistic
            comparison = bootstrap_mannwhitney_statistics >= mw_statistic
            print("Comparison array sample:", comparison[:20])
            print("Sum of comparison array:", np.sum(comparison))               

        print(f"Bootstrapped Mann-Whitney U test {groupA_name} vs {groupB_name} - p-value: {bootstrap_mw_p_value}")
        # Compute confidence intervals for the bootstrap distribution
        bootstrap_ci_low = np.percentile(bootstrap_mannwhitney_statistics, 2.5)
        bootstrap_ci_high = np.percentile(bootstrap_mannwhitney_statistics, 97.5)

        print(f"Bootstrap CI for Mann-Whitney U statistic - Low: {bootstrap_ci_low}, High: {bootstrap_ci_high}")

        bootstrap_p_value.append(bootstrap_mw_p_value)
        bootstrap_ci_lows.append(bootstrap_ci_low)
        bootstrap_ci_highs.append(bootstrap_ci_high)
        i=i+1   


    
    # Check lengths of each array in 'groups_flatt'
    for key, value in groups_flatt.items():
        print(f"Length of {key}: {len(value)}")

    # Determine the maximum length of the arrays
    max_length = max(len(value) for value in groups_flatt.values())

    # Pad the arrays with None (or np.nan) to ensure all arrays have the same length
    for key, value in groups_flatt.items():
        if len(value) < max_length:
            padding = np.full(max_length - len(value), np.nan)  # Using np.nan for padding
            groups_flatt[key] = np.concatenate((value, padding))

    # Convert to DataFrame
    results_groups = pd.DataFrame(groups_flatt)
    results_groups.to_csv(f'/data/p_02915/SPOT/02_averaged_hemi-{hemi}_ROIs_dvars.csv')

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
results_df.to_csv(f'/data/p_02915/SPOT/01_results_ROI_dvars.csv')
