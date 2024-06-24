'''
Statistical comparison of MSD between emperical and simulated across the hemispheres
'''
import pandas as pd
import numpy as np
from scipy.stats import kruskal, wilcoxon, norm
from multiprocessing import Pool, cpu_count
from nilearn import surface
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

def get_indices_roi(labels_area, visparc_array):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (numpy.array): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
        np.logical_or(visparc_array == labels_area[0], visparc_array == labels_area[1])
    )[0]
    return indices_area


def calc_msd(param_area1, param_area2):
    # mean squared difference between two different model results in an area
    n1 = param_area1.size
    n2 = param_area2.size
    assert (
        n1 == n2
    ), "areas don't have the same number of vertices and cannot be compared!"

    return np.sum((param_area1 - param_area2) ** 2) / n1

def bootstrap_resample_wilcoxon(n_resamples, combined_data, n_groupA, n_groupB):
    bootstrap_mannwhitney_statistics = []
    for _ in range(n_resamples):
        # Resample the combined data
        resample_indices = np.random.choice(len(combined_data), size=len(combined_data), replace=True)
        resampled_combined_data = combined_data[resample_indices]
        
        # Split the resampled data into two groups
        bootstrap_groupA = flatten(resampled_combined_data[:n_groupA])
        bootstrap_groupB = flatten(resampled_combined_data[n_groupA:])
        
        # Perform Mann-Whitney U test on the bootstrap sample
        bootstrap_mw_statistic,_ = wilcoxon(bootstrap_groupA, bootstrap_groupB, alternative="two-sided")
        bootstrap_mannwhitney_statistics.append(bootstrap_mw_statistic)
    return bootstrap_mannwhitney_statistics

# Function to convert elements to float, setting non-convertible elements to np.nan
def to_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.nan
    
def flatten(arr):
    return arr.flatten()

LABELS_V1 = (1, 2)
LABELS_V2 = (3, 4)

df = pd.DataFrame()
for param in ["eccentricity", "polarangle"]:
    print(param)
    i=0    
    real_mean=[]
    simul_mean=[]
    original_H_value=[]
    original_p_value=[]
    original_z_score=[]
    bootstrap_p_value=[]
    bootstrap_ci_lows=[]
    bootstrap_ci_highs=[]
    index_label = []
    for group in ["12-16y", "18-21y", "neonates<37", "neonates>37", "fetal<29", "fetal>29"]: 
        results_list = []    
        if group == "neonates<37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-{model}_v0i.gii"
            )
            PREFIX_SUB_TEMPLATE = (
                "/data/p_02915/dhcp_derivatives_SPOT/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
            )

            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
            sub_num = len(subject_info["sub_id"])   
        elif group == "neonates>37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-{model}_v0i.gii"
            )
            PREFIX_SUB_TEMPLATE = (
                "/data/p_02915/dhcp_derivatives_SPOT/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
            )
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "fetal<29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-{model}_v0i.gii"
            )
            PREFIX_SUB_TEMPLATE = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
            )
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
            sub_num = len(subject_info["sub_id"])        
        elif group == "fetal>29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native_desc-{model}_v0i.gii"
            )
            PREFIX_SUB_TEMPLATE = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/dhcp_surface/sub-{sub}/ses-{ses}/anat/"
                "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-native_dens-native"
            )
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "12-16y":       
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
                "{sub}_hemi-{hemi}_mesh-native_dens-native_desc-{model}_v0i.gii"
            )
            PREFIX_SUB_TEMPLATE = (
                "/data/p_02915/SPOT/00_HCP/template/"
                "hemi-{hemi}_mesh-native_dens-native"
            )
            subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
            sub_num = len(subject_info["sub_id"])
        elif group == "18-21y":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/{sub}/"
                "{sub}_hemi-{hemi}_mesh-native_dens-native_desc-{model}_v0i.gii"
            )
            PREFIX_SUB_TEMPLATE = (
                "/data/p_02915/SPOT/00_HCP/template/"
                "hemi-{hemi}_mesh-native_dens-native"
            )            
            subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
            sub_num = len(subject_info["sub_id"])
                
        PATH_ECCENTRICITY = "{prefix_sub_template}_desc-eccentretinotbenson2014_seg.shape.gii"
        PATH_ANGLE = "{prefix_sub_template}_desc-angleretinotbenson2014_seg.shape.gii"
        PATH_VISPARC = (
            "{prefix_sub_template}_desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
        )
        OUTPUT_ECC = (
            "{prefix_model}_label-eccentricity_desc-{model}_roi-v2th{threshold}_metric.gii"
        )
        OUTPUT_PANGLE = (
            "{prefix_model}_label-polarangle_desc-{model}_roi-v2th{threshold}_metric.gii"
        )

        group_name=f"{group}"
        parameters=[]
        parameters2=[]
        # FORMAT PATHS FOR INPUT AND OUTPUT       
        for index, row in subject_info.iterrows():                
            if group == "12-16y" or group == "18-21y":                
                for model in ["real", "simulated"]:
                    for hemi in ["L", "R"]:
                        hemi_results = {}  # Results for the current hemisphere
                        sub_id = subject_info.at[index, "sub_id"]
                        sub = sub_id.replace('sub-','')
                        # FORMAT PATHS FOR INPUT AND OUTPUT
                        prefix_model = PREFIX_MODEL.format(sub=sub, hemi=hemi, model = model, param = param)
                        prefix_sub_template = PREFIX_SUB_TEMPLATE.format(sub=sub, hemi=hemi)
                        # LOAD DATA
                        visparc = surface.load_surf_data(
                            PATH_VISPARC.format(prefix_sub_template=prefix_sub_template)
                        )
                        indices_v1 = get_indices_roi(LABELS_V1, visparc)
                        indices_v2 = get_indices_roi(LABELS_V2, visparc)
                        if hemi=="L":
                            if param == "eccentricity":
                                benson = surface.load_surf_data(PATH_ECCENTRICITY.format(prefix_sub_template=prefix_sub_template))
                                
                            elif param == "polarangle":
                                benson = surface.load_surf_data(PATH_ANGLE.format(prefix_sub_template=prefix_sub_template))
                                
                            ccf_v0 = surface.load_surf_data(prefix_model).astype(int, copy=False)
                            centers_v2 = ccf_v0[indices_v2]
                             # indices to voxels in v1 array that are the centers for each v2 voxel
                            # restrict to V1
                            ret_v1 = benson[indices_v1]                      
                            # project data from V1 vertices to V2
                            ret_v2 = ret_v1[centers_v2]
                            # make empty wholebrain
                            ret_wholeb = np.zeros(benson.shape)
                            # fill retinotopy values into v2
                            ret_wholeb[indices_v2] = ret_v2
                            # keep data for model comparisons later
                            benson_L = benson[indices_v2]
                            if model == "real":
                                real_L = ret_wholeb[indices_v2]
                            elif model == "simulated":
                                simulated_L = ret_wholeb[indices_v2]

                        elif hemi=="R":
                            if param == "eccentricity":
                                benson = surface.load_surf_data(PATH_ECCENTRICITY.format(prefix_sub_template=prefix_sub_template))
                               
                            elif param == "polarangle":
                                benson = surface.load_surf_data(PATH_ANGLE.format(prefix_sub_template=prefix_sub_template))
                                
                            ccf_v0 = surface.load_surf_data(prefix_model).astype(int, copy=False)
                            centers_v2 = ccf_v0[indices_v2]
                             # indices to voxels in v1 array that are the centers for each v2 voxel
                            # restrict to V1
                            ret_v1 = benson[indices_v1]                      
                            # project data from V1 vertices to V2
                            ret_v2 = ret_v1[centers_v2]
                            # make empty wholebrain
                            ret_wholeb = np.zeros(benson.shape)
                            # fill retinotopy values into v2
                            ret_wholeb[indices_v2] = ret_v2
                            # keep data for model comparisons later
                            benson_R = benson[indices_v2]
                            if model == "real":
                                real_R = ret_wholeb[indices_v2]
                            elif model == "simulated":
                                simulated_R = ret_wholeb[indices_v2]
                        
                benson_wb = np.hstack((benson_L, benson_R))
                real_wb = np.hstack((real_L, real_R))
                simulated_wb = np.hstack((simulated_L, simulated_R))
                msd_r_b = calc_msd(real_wb, benson_wb)
                msd_r_s = calc_msd(real_wb, simulated_wb)
                
                parameters.append(msd_r_b)        
                parameters2.append(msd_r_s)   
            else:   
                for model in ["real", "simulated"]:
                    for hemi in ["L", "R"]:
                        hemi_results = {}  # Results for the current hemisphere
                        sub_id = subject_info.at[index, "sub_id"]
                        sess_id = subject_info.at[index, "sess_id"]
                        sub = sub_id.replace('sub-','')
                        ses = sess_id.replace('ses-','')
                        # FORMAT PATHS FOR INPUT AND OUTPUT
                        prefix_model = PREFIX_MODEL.format(sub=sub, ses=ses, hemi=hemi, model = model, param = param)
                        prefix_sub_template = PREFIX_SUB_TEMPLATE.format(sub=sub,ses=ses, hemi=hemi)
                        # LOAD DATA
                        visparc = surface.load_surf_data(
                            PATH_VISPARC.format(prefix_sub_template=prefix_sub_template)
                        )
                        indices_v1 = get_indices_roi(LABELS_V1, visparc)
                        indices_v2 = get_indices_roi(LABELS_V2, visparc)
                        if hemi=="L":
                            if param == "eccentricity":
                                benson = surface.load_surf_data(PATH_ECCENTRICITY.format(prefix_sub_template=prefix_sub_template))
                                
                            elif param == "polarangle":
                                benson = surface.load_surf_data(PATH_ANGLE.format(prefix_sub_template=prefix_sub_template))
                                
                            ccf_v0 = surface.load_surf_data(prefix_model).astype(int, copy=False)
                            centers_v2 = ccf_v0[indices_v2]
                             # indices to voxels in v1 array that are the centers for each v2 voxel
                            # restrict to V1
                            ret_v1 = benson[indices_v1]                      
                            # project data from V1 vertices to V2
                            ret_v2 = ret_v1[centers_v2]
                            # make empty wholebrain
                            ret_wholeb = np.zeros(benson.shape)
                            # fill retinotopy values into v2
                            ret_wholeb[indices_v2] = ret_v2
                            # keep data for model comparisons later
                            benson_L = benson[indices_v2]
                            if model == "real":
                                real_L = ret_wholeb[indices_v2]
                            elif model == "simulated":
                                simulated_L = ret_wholeb[indices_v2]

                        elif hemi=="R":
                            if param == "eccentricity":
                                benson = surface.load_surf_data(PATH_ECCENTRICITY.format(prefix_sub_template=prefix_sub_template))
                               
                            elif param == "polarangle":
                                benson = surface.load_surf_data(PATH_ANGLE.format(prefix_sub_template=prefix_sub_template))
                                
                            ccf_v0 = surface.load_surf_data(prefix_model).astype(int, copy=False)
                            centers_v2 = ccf_v0[indices_v2]
                             # indices to voxels in v1 array that are the centers for each v2 voxel
                            # restrict to V1
                            ret_v1 = benson[indices_v1]                      
                            # project data from V1 vertices to V2
                            ret_v2 = ret_v1[centers_v2]
                            # make empty wholebrain
                            ret_wholeb = np.zeros(benson.shape)
                            # fill retinotopy values into v2
                            ret_wholeb[indices_v2] = ret_v2
                            # keep data for model comparisons later
                            benson_R = benson[indices_v2]
                            if model == "real":
                                real_R = ret_wholeb[indices_v2]
                            elif model == "simulated":
                                simulated_R = ret_wholeb[indices_v2]
                        
                benson_wb = np.hstack((benson_L, benson_R))
                real_wb = np.hstack((real_L, real_R))
                simulated_wb = np.hstack((simulated_L, simulated_R))
                msd_r_b = calc_msd(real_wb, benson_wb)
                msd_r_s = calc_msd(real_wb, simulated_wb)
                
                parameters.append(msd_r_b)        
                parameters2.append(msd_r_s)  
        print(group_name)   
        print(np.mean(parameters))  
        print(np.mean(parameters2))           
        real_mean.append(np.mean(parameters))
        simul_mean.append(np.mean(parameters2))
        
    
        # Perform Mann-Whitney U tests on each pair of groups
        wilcoxon_results = []
        groupA = parameters
        groupB = parameters2
        # Combine data for Kruskal-Wallis test
        combined_data = np.concatenate([groupA, groupB])
        n_groupA = len(groupA)
        n_groupB = len(groupB)
        
        res = wilcoxon(groupA, groupB, alternative='two-sided', method="approx")
        print(res.zstatistic)
        mw_statistic = res.statistic
        mw_p_value = res.pvalue        
        z_score=res.zstatistic
        wilcoxon_results.append(("real_benson", "real_simulated", mw_statistic, mw_p_value))
        print(f"Original wilcoxon test real_benson vs real_simulated - U statistic: {mw_statistic}, p-value: {mw_p_value}")
        #z_score = mann_whitney_u_to_z(mw_statistic,n_groupA, n_groupB)
        eff_r = z_score / np.sqrt(n_groupA)
        index_label.append(f"{group_name}")
        original_H_value.append(eff_r)
        original_z_score.append(z_score)
        original_p_value.append(mw_p_value)
        n_bootstraps=10000
        n_cores = cpu_count()
        #print(n_cores)
        n_bootstraps_per_core=n_bootstraps // n_cores
        #print(n_bootstraps_per_core)


        with Pool(n_cores) as pool:
        # Split the work among the cores
            results = pool.starmap(bootstrap_resample_wilcoxon, [(n_bootstraps_per_core, combined_data, n_groupA, n_groupB)] * n_cores)
        
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
        
        print(f"Bootstrapped Mann-Whitney U test benson vs simulated - p-value: {bootstrap_mw_p_value}")
        # Compute confidence intervals for the bootstrap distribution
        bootstrap_ci_low = np.percentile(bootstrap_mannwhitney_statistics, 2.5)
        bootstrap_ci_high = np.percentile(bootstrap_mannwhitney_statistics, 97.5)

        print(f"Bootstrap CI for Mann-Whitney U statistic - Low: {bootstrap_ci_low}, High: {bootstrap_ci_high}")

        bootstrap_p_value.append(bootstrap_mw_p_value)
        bootstrap_ci_lows.append(bootstrap_ci_low)
        bootstrap_ci_highs.append(bootstrap_ci_high)

    results_dict = {
    "effect size r": original_H_value,
    "original z": original_z_score,
    'Original p-value': original_p_value,
    'Bootstrapped p-value': bootstrap_p_value,
    'Bootstrap CI low': bootstrap_ci_lows,
    'Bootstrap CI high': bootstrap_ci_highs,
    "real_mean": real_mean,
    "simul_mean": simul_mean
    }

    # Create a DataFrame from the dictionary
    results_df = pd.DataFrame(results_dict, index=index_label)
    # Save the DataFrame to a CSV file
    results_df.to_csv(f'/data/p_02915/SPOT/01_results_MSD_wb_{param}.csv')
