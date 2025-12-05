import pandas as pd
import numpy as np

# Ensure combat_harmonized is defined or loaded here
# combat_harmonized = ...

# Loop over each area (V1, V2, V3)
for area in ["v1", "v2", "v3"]:
    # Load left and right hemisphere data for each area
    left = pd.read_csv(f"/data/p_02915/SPOT/covars_hemi-L_tsnr_{area}.csv",index_col=None)
    left = np.mean(left, axis=1)
    
    right = pd.read_csv(f"/data/p_02915/SPOT/covars_hemi-R_tsnr_{area}.csv", index_col=None)
    right = np.mean(right, axis=1)

    # Combine left and right hemisphere data into one DataFrame
    combat_harmonized = pd.concat([left, right], axis=1)
    combat_harmonized=combat_harmonized.T.to_numpy()

    # Load covariate data
    covars = pd.read_csv(f"/data/p_02915/SPOT/covars_hemi-L.csv")

    # Group data into categories
    groups_flatt = {
        'tri-2 prenatal': combat_harmonized[:, covars[covars["group"] == "2nd"].index].flatten(),
        'tri-3 prenatal': combat_harmonized[:, covars[covars["group"] == "3rd"].index].flatten(),
        'preterm neonatal': combat_harmonized[:, covars[covars["group"] == "preterm"].index].flatten(),
        'fullterm neonatal': combat_harmonized[:, covars[covars["group"] == "fullterm"].index].flatten(),
        'adolescent': combat_harmonized[:, covars[covars["group"] == "adolescent"].index].flatten(),
        'adult': combat_harmonized[:, covars[covars["group"] == "adult"].index].flatten()
    }

    # Initialize lists to hold means and standard deviations
    means = {}
    std_devs = {}

    # Calculate mean and standard deviation for each group
    for group_name, values in groups_flatt.items():
        means[group_name] = np.mean(values)
        std_devs[group_name] = np.std(values)

    # Convert results to DataFrame for better readability
    stats_df = pd.DataFrame({
        'Region': [area] * len(means),  # Same region for all groups
        'Group': list(means.keys()),  # Group names
        'Mean': list(means.values()),  # Mean values
        'Standard Deviation': list(std_devs.values())  # Standard deviations
    })

    # Display the statistics DataFrame
    print(stats_df)
