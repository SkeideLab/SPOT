'''
Statistical comparison of correlation coefficient between emperical and simulated 
'''
import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from scipy.stats import kruskal, mannwhitneyu
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
        bootstrap_H_statistic,_ = kruskal(*resampled_groups)
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
            visparc.agg_data() == labels_area[0], visparc.agg_data() == labels_area[1]
        )
    )[0]
    return indices_area

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


VISPARC_PATH = ("/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")
LABELS_V2 = (3, 4)
 
i=0
original_H_value=[]
original_z_score=[]
original_p_value=[]
bootstrap_p_value=[]
bootstrap_ci_lows=[]
bootstrap_ci_highs=[]
index_label = []
df = pd.DataFrame()
for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]:
    results_list = []    
    for hemi in ["L", "R"]: 
        for model in ["real", "simulated"]:
            if group == "neonates<37":
                PREFIX_MODEL = (
                    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
                    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-{model}_r.gii"
                )
                subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
                sub_num = len(subject_info["sub_id"])   
            elif group == "neonates>37":
                PREFIX_MODEL = (
                    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
                    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-{model}_r.gii"
                )
                subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "fetal<29":
                PREFIX_MODEL = (
                    "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
                    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-{model}_r.gii"
                )
                subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
                sub_num = len(subject_info["sub_id"])        
            elif group == "fetal>29":
                PREFIX_MODEL = (
                    "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
                    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-{model}_r.gii"
                )
                subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "12-16y":
                PREFIX_MODEL = (
                    "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
                    "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-{model}_r.gii"
                )
                subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
                sub_num = len(subject_info["sub_id"])
            elif group == "18-21y":
                PREFIX_MODEL = (
                    "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
                    "{sub}_hemi-{hemi}_mesh-fsaverage_dens-164k_desc-{model}_r.gii"
                )
                subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
                sub_num = len(subject_info["sub_id"])

        
            visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
            indices_v2 = get_indices_roi(LABELS_V2, visparc)
            group_name=f"{group}_{hemi}"
            parameters=[]
        # FORMAT PATHS FOR INPUT AND OUTPUT       
            for index, row in subject_info.iterrows():
                if group == "12-16y" or group == "18-21y":
                    sub_id = subject_info.at[index, "sub_id"]
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub_id,
                        hemi=hemi,
                        model=model,
                    )
                    input_path = prefix_model

                    ccf = surface.load_surf_data(input_path)
                    ccf_v0 = ccf[indices_v2].astype(np.float64)
                    parameters.append(ccf_v0)           
                else:
                    sub_id = subject_info.at[index, "sub_id"]
                    sess_id = subject_info.at[index, "sess_id"]
                    sub = sub_id.replace('sub-','')
                    ses = sess_id.replace('ses-','')
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                        model=model,
                    )
                    input_path = prefix_model

                    ccf = surface.load_surf_data(input_path)
                    ccf_v0 = ccf[indices_v2].astype(np.float64)
                    parameters.append(ccf_v0)           

            if model == "real":
                corr_real = np.array(parameters)
            elif model == "simulated":
                corr_simul = np.array(parameters)

        n_bootstraps=10000

        mannwhitney_results = []
        real_flatt = corr_real.flatten()
        simul_flatt = corr_simul.flatten()
        
        # Combine data for Kruskal-Wallis test
        combined_data = np.concatenate([corr_real, corr_simul])
        n_groupA = len(corr_real)
        n_groupB = len(corr_simul)
        
        groupA_name = f"{hemi}_{group}_real"
        groupB_name = f"{hemi}_{group}_simul"
        mw_statistic, mw_p_value = mannwhitneyu(real_flatt, simul_flatt, alternative='two-sided')
        mannwhitney_results.append((groupA_name, groupB_name, mw_statistic, mw_p_value))
        print(f"Original Mann-Whitney U test {groupA_name} vs. {groupB_name} - U statistic: {mw_statistic}, p-value: {mw_p_value}")
        index_label.append(f"{groupA_name} vs. {groupB_name}")
        z_score=mann_whitney_u_to_z(mw_statistic, len(real_flatt), len(simul_flatt))
        print(f"z-score: {z_score}")
        original_H_value.append(mw_statistic)
        original_z_score.append(z_score)
        original_p_value.append(mw_p_value)

        n_cores = cpu_count()
        #print(n_cores)
        n_bootstraps_per_core=n_bootstraps // n_cores

        with Pool(n_cores) as pool:
        # Split the work among the cores
            results = pool.starmap(bootstrap_resample_mannwhitney, [(n_bootstraps_per_core, combined_data, n_groupA, n_groupB)] * n_cores)

        # Combine results from all cores
        bootstrap_mannwhitney_statistics = np.concatenate(results)
        # Calculate bootstrap p-value for Mann-Whitney U test
        bootstrap_mw_p_value = np.sum(bootstrap_mannwhitney_statistics >= mw_statistic) / n_bootstraps
        
        print(f"Bootstrapped Mann-Whitney U test {groupA_name} vs {groupB_name} - p-value: {bootstrap_mw_p_value}")
        # Compute confidence intervals for the bootstrap distribution
        bootstrap_ci_low = np.percentile(bootstrap_mannwhitney_statistics, 2.5)
        bootstrap_ci_high = np.percentile(bootstrap_mannwhitney_statistics, 97.5)

        print(f"Bootstrap CI for Mann-Whitney U statistic - Low: {bootstrap_ci_low}, High: {bootstrap_ci_high}")

        bootstrap_p_value.append(bootstrap_mw_p_value)
        bootstrap_ci_lows.append(bootstrap_ci_low)
        bootstrap_ci_highs.append(bootstrap_ci_high)
        i=i+1

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
results_df.to_csv(f'/data/p_02915/SPOT/01_results_correlation_coefficient.csv')