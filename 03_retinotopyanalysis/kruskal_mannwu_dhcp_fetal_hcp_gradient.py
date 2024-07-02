import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from scipy.stats import kruskal
from multiprocessing import Pool, cpu_count
from nilearn import surface


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

VISPARC_PATH = ("/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")

def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
    indices_area = np.nonzero(
            visparc.agg_data() == labels_area
    )[0]
    return indices_area

def extract_roi_mesh(brain_mesh, indices):
    # Extract the vertices and faces corresponding to the ROI
    mask = np.zeros(brain_mesh.n_points, dtype=bool)
    mask[indices] = True

    # Extract the ROI mesh
    roi_mesh = brain_mesh.extract_points(mask, include_cells=True)
    return roi_mesh

# Create a plotter
#plotter = pv.Plotter(shape=(6, 4),border=False, off_screen=True)

df = pd.DataFrame()
for param in ["eccentricity", "polarangle"]:
    results_list = []    
    print(param)            
    gradients={}
    for group in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]:
        original_H_value=[]
        original_z_score=[]
        original_p_value=[]
        bootstrap_p_value=[]
        bootstrap_ci_lows=[]
        bootstrap_ci_highs=[]
        index_label = []
        if group == "neonates<37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/"
                "Averaged_{hemi}_label-{param}_desc-real_roi-v2th00_metric_less_37.gii"
            )
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
            sub_num = len(subject_info["sub_id"])   
            row_p = 2
        elif group == "neonates>37":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/"
                "Averaged_{hemi}_label-{param}_desc-real_roi-v2th00_metric_less_37.gii"
            )
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 3
        elif group == "fetal<29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/"
                "Averaged_younger_fetal_{hemi}_label-{param}_desc-real_roi-v2th00_metric.gii"
            )
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
            sub_num = len(subject_info["sub_id"])        
            row_p = 0
        elif group == "fetal>29":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{ses}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/"
                "Averaged_older_fetal_{hemi}_label-{param}_desc-real_roi-v2th00_metric.gii"
            )
            subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 1
        elif group == "12-16y":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/"
                "Averaged_{hemi}_label-{param}_desc-real_roi-v2th00_metric_young.gii"
            )
            subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_young.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 4
        elif group == "18-21y":
            PREFIX_MODEL = (
                "/data/p_02915/dhcp_derivatives_SPOT/statistics/"
                "{sub}_{hemi}_{param}_gradient.gii"
            )
            AVER_MAP = (
                "/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ccfmodel/"
                "Averaged_{hemi}_label-{param}_desc-real_roi-v2th00_metric_old.gii"
            )
            subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
            sub_num = len(subject_info["sub_id"])
            row_p = 5
        
        for axis in ["R","A","S"]:
            if axis == "R":
                ax = 0
            elif axis == "A":
                ax = 1
            elif axis == "S":
                ax = 2

            for hemi in ["L", "R"]:
                for region in ["v2v", "v2d"]:    
                    if region == "v2v":
                        LABELS_V2 = 3
                    elif region == "v2d":
                        LABELS_V2 = 4            
                    group_name=f"{group}_{hemi}"
                    parameters=[]
                    visparc = nib.load(VISPARC_PATH.format(hemi=hemi))
                    indices_v2 = get_indices_roi(LABELS_V2, visparc)
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
                            ccf_v0 = ccf[ax,indices_v2].astype(np.float64)
                            ccf_v0 = ccf_v0[~np.isnan(ccf_v0)]
                            parameters.append(np.mean(ccf_v0))     
                        else:
                            sub_id = subject_info.at[index, "sub_id"]
                            sess_id = subject_info.at[index, "sess_id"]
                            sub = sub_id.replace('sub-','')
                            ses = sess_id.replace('ses-','')
                            prefix_model = PREFIX_MODEL.format(
                                sub=sub,
                                ses=ses,
                                hemi=hemi,
                                param=param,
                            )                            
                            input_path = prefix_model
                            ccf = surface.load_surf_data(input_path)                            
                            ccf_v0 = ccf[ax,indices_v2].astype(np.float64)
                            ccf_v0 = ccf_v0[~np.isnan(ccf_v0)]
                            parameters.append(np.mean(ccf_v0))    
                    gradients[f"{group}_{hemi}_{region}_{axis}"] = np.array(parameters)

    for hemi in ["L", "R"]:
        for region in ["v2v", "v2d"]:
            for axis in ["R","A","S"]:
                n_bootstraps=10000
                groups_flatt = {
                'fetal_y': gradients[f"fetal<29_{hemi}_{region}_{axis}"],
                'fetal_o': gradients[f"fetal>29_{hemi}_{region}_{axis}"],
                'neonates_y': gradients[f"neonates<37_{hemi}_{region}_{axis}"],
                'neonates_o': gradients[f"neonates>37_{hemi}_{region}_{axis}"],
                'hcp_y': gradients[f"12-16y_{hemi}_{region}_{axis}"],
                'hcp_o': gradients[f"18-21y_{hemi}_{region}_{axis}"],
                }

                all_groups_flatt = list(groups_flatt.values())      
                # Perform the original Kruskal-Wallis test
                kruskal_statistic, kruskal_p_value = kruskal(*all_groups_flatt)
                print(hemi)
                print(f"Original Kruskal-Wallis statistic: {kruskal_statistic}, p-value: {kruskal_p_value}")
                k = 6  # Number of groups
                N = N = sum(len(group) for group in all_groups_flatt)   # Total number of observations
                print(N)
                eta_squared = np.sqrt(kruskal_statistic / N)
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
                index_label.append(f"{hemi}_{region}_{axis}")
    
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
    results_df.to_csv(f'/data/p_02915/SPOT/01_results_{param}_gradient.csv')


