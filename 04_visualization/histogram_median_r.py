import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
# Get the screen dimensions using tkinter
import tkinter as tk
root = tk.Tk()
screen_width_px = root.winfo_screenwidth()
screen_height_px = root.winfo_screenheight()
root.destroy()

# Convert pixels to inches (assuming 100 dpi for conversion)
dpi = 100
screen_width_in = screen_width_px / dpi
screen_height_in = screen_height_px / dpi

# Set desired figure size to half of screen width and full screen height
fig_width_in = screen_width_in / 2
fig_height_in = screen_height_in

colors = {'neonates<37': 'limegreen', 'neonates>37': 'mediumaquamarine', 'fetal<29': 'gold', 'fetal>29': 'yellowgreen', '12-16y': 'cadetblue', '18-21y': 'steelblue'}
text_label = {'neonates<37': 'Preterm neonatal', 'neonates>37': 'Full-term neonatal', 'fetal<29': 'Fetal second trimester', 'fetal>29': 'Fetal third trimester', '12-16y': 'Adolescent', '18-21y': 'Adult'}
for hemi in ["left", "right"]:  
    for param in ["real", 'simulated']:
        
        df = pd.DataFrame()
        for model in ["neonates<37", "neonates>37", "fetal<29", "fetal>29", "12-16y", "18-21y"]:
        
            if model =="neonates<37":
                median_r = pd.read_csv("/data/p_02915/SPOT/median_r_less_37.csv")
            elif model == "neonates>37":
                median_r = pd.read_csv("/data/p_02915/SPOT/median_r_over_37.csv")
            elif model == "fetal>29w":
                median_r = pd.read_csv("/data/p_02915/SPOT/median_r_fetal_old.csv")
            elif model == "fetal<29w":
                median_r = pd.read_csv("/data/p_02915/SPOT/median_r_fetal_young.csv")    
            elif model == "12-16y":
                median_r = pd.read_csv("/data/p_02915/SPOT/median_r_HCP_old.csv")
            elif model == "18-21y":
                median_r = pd.read_csv("/data/p_02915/SPOT/median_r_HCP_young.csv")            
                    
            
            data_median_r = pd.DataFrame({'median_r': median_r[f"{hemi}_{param}"], 'Group': model})
            df = pd.concat([df,data_median_r], ignore_index=True)
            
         # Calculate the common bin edges for all groups
        all_differences = df['median_r']
        bin_edges = np.histogram_bin_edges(all_differences, bins=15)

        # Create subplots
        fig, axes = plt.subplots(6, 1, sharex=True, figsize=(fig_width_in, fig_height_in))

        # Plot normalized histograms for each group
        for ax, group in zip(axes, df['Group'].unique()):
            print(group)
            sns.histplot(data=df[df['Group'] == group], x='median_r', ax=ax, bins=bin_edges, stat='density', color=colors[group])
            ax.text(-0.3, 0.5, text_label[group], transform=ax.transAxes, 
                fontsize=10, verticalalignment='center', rotation='horizontal')
            ax.set_ylabel('')     

        # Set common x-axis label
        axes[-1].set_xlabel('Median correlation (r)')

        fig.suptitle(f'Median correlation of {param} on {hemi} hemisphere', fontsize=16)
        
        # Adjust layout and spacing between subplots
        plt.tight_layout(h_pad=0.5)

        plt.show()