import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

msd_value = pd.read_csv("/data/p_02915/SPOT/MSD.csv")

for hemi in ["L", "R"]:
    for param in ["eccentricity", "polarangle"]:
        Benson = msd_value.loc[:,f"{hemi}_{param}_real_benson"]
        Simulated = msd_value.loc[:,f"{hemi}_{param}_real_simulated"]
        differences = Benson - Simulated
        num_bins = len(differences)

        fig, ax = plt.subplots()

        # the histogram of the data
        n, bins, patches = ax.hist(differences, num_bins)

        ax.set_xlabel('MSD differences', fontsize = 20)
        ax.set_ylabel('Number of subject', fontsize = 20)
        ax.set_title(f"Histogram of MSD_Benson - MSD_Simulated data on hemisphere {hemi} for {param}", fontsize = 20)
        ax.tick_params(axis='x', labelsize=15)  # Adjust the fontsize as needed
        ax.tick_params(axis='y', labelsize=15)  # Adjust the fontsize as needed
   
        plt.show()