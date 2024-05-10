import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("/data/p_02915/SPOT/median_r.csv")

for hemi in ["left", "right"]:
        
    datasource = "real"
    median_correlation = np.median(df[f"{hemi}_{datasource}"])
    mean_median = np.mean(df[f"{hemi}_{datasource}"])
    print(f"Median correlation of {datasource} data on hemisphere {hemi} is {median_correlation}.")
    num_bins = len(df[f"{hemi}_{datasource}"])

    fig, ax = plt.subplots()

    # the histogram of the data
    n, bins, patches = ax.hist(df[f"{hemi}_{datasource}"], num_bins, density=True)

    ax.set_xlabel('Median correlation')
    ax.set_ylabel('Number of subject')
    ax.set_title(f"Histogram of {datasource} on {hemi} hemi: median correlation = {median_correlation}.")
    
    plt.show()
