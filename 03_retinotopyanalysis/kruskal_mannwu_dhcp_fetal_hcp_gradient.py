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


VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii")


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
# plotter = pv.Plotter(shape=(6, 4),border=False, off_screen=True)


df = pd.DataFrame()
for param in ["eccentricity", "polarangle"]:    
    i = 0
    original_H_value = []
    original_z_score = []
    original_p_value = []
    bootstrap_p_value = []
    bootstrap_ci_lows = []
    bootstrap_ci_highs = []
    index_label = []

    for hemi in ["L", "R"]:
        for area in ["V2v", "V2d", "V3v", "V3d"]:
            for axis in ["R", "A", "S"]:
                n_bootstraps = 10000
                combat_harmonized = pd.read_csv(f"/data/p_02915/SPOT/combat_hemi-{hemi}_area-{area}_{param}_{axis}.csv", index_col=None).to_numpy()
                covars = pd.read_csv(f"/data/p_02915/SPOT/covars_hemi-{hemi}.csv")
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
                
                # Perform the original Kruskal-Wallis test
                kruskal_statistic, kruskal_p_value = kruskal(*all_groups_flatt)
                print(hemi)
                print(
                    f"Original Kruskal-Wallis statistic: {kruskal_statistic}, p-value: {kruskal_p_value}")
                k = 6  # Number of groups
                # Total number of observations
                N = N = sum(len(group) for group in all_groups_flatt)
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
                # print(n_cores)
                n_bootstraps_per_core = n_bootstraps // n_cores
                # print(n_bootstraps_per_core)

                with Pool(n_cores) as pool:
                    results = pool.starmap(bootstrap_kruskal_resample, [(
                        n_bootstraps_per_core, group_sizes, combined_data)] * n_cores)

                # Combine results from all cores
                bootstrap_kruskal_statistics = np.concatenate(results)

                # Calculate the bootstrap p-value for Kruskal-Wallis test
                bootstrap_kruskal_statistics = np.array(
                    bootstrap_kruskal_statistics)
                bootstrap_kruskal_p_value = np.sum(
                    bootstrap_kruskal_statistics >= kruskal_statistic) / n_bootstraps

                print("\nBootstrapped Kruskal-Wallis H-statistic:",
                      kruskal_statistic)
                print("Bootstrapped p-value:", bootstrap_kruskal_p_value)
                bootstrap_ci_low = np.percentile(
                    bootstrap_kruskal_statistics, 2.5)
                bootstrap_ci_high = np.percentile(
                    bootstrap_kruskal_statistics, 97.5)

                print(
                    f"Bootstrap CI for Kruskal-Wallis H statistic - Low: {bootstrap_ci_low}, High: {bootstrap_ci_high}")

                bootstrap_p_value.append(bootstrap_kruskal_p_value)
                bootstrap_ci_lows.append(bootstrap_ci_low)
                bootstrap_ci_highs.append(bootstrap_ci_high)
                index_label.append(f"{hemi}_{area}_{axis}")

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
