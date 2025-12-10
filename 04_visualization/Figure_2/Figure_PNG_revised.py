import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from neuromaps.datasets import fetch_fsaverage
import matplotlib.pyplot as plt
from nilearn import plotting
import seaborn as sns
from matplotlib import cm
import scipy.stats as stats
import os

# Get the full path of the script
file_path = os.path.dirname(os.path.abspath(__file__))

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

def r_squared_significance(R2, n, k):
    """
    Test the significance of R-squared using an F-test.
    
    Parameters:
    R2 : float - R-squared value from the regression
    n  : int   - Number of observations (sample size)
    k  : int   - Number of predictors (independent variables)
    
    Returns:
    F-statistic and p-value
    """
    if R2 == 1:  # To avoid division by zero error
        return float('inf'), 0.0  
    
    # Compute the F-statistic
    F_stat = (R2 / k) / ((1 - R2) / (n - k - 1))
    
    # Compute p-value
    p_value = 1 - stats.f.cdf(F_stat, dfn=k, dfd=n - k - 1)
    
    return F_stat, p_value

VISPARC_PATH = (
    "{file_path}/hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii")
LABELS_V2 = [2, 3]


left_temp = nib.load(VISPARC_PATH.format(file_path=file_path, hemi="L"))
indices_left = get_indices_roi([2, 3], left_temp)
indices_left_v2 = get_indices_roi([2], left_temp)
indices_left_v3 = get_indices_roi([3], left_temp)
right_temp = nib.load(VISPARC_PATH.format(file_path=file_path, hemi="R"))
indices_right = get_indices_roi([2, 3], right_temp)
indices_right_v2 = get_indices_roi([2], right_temp)
indices_right_v3 = get_indices_roi([3], right_temp)

df = pd.DataFrame()
surfaces = fetch_fsaverage(density="164k")
lh, rh = surfaces["inflated"]
# Load sulcal depth data
# Load the left hemisphere sulcal depth data
lh_sulc = surface.load_surf_data(
    f"{file_path}/lh.curv")
# Load the right hemisphere sulcal depth data
rh_sulc = surface.load_surf_data(
    f"{file_path}/rh.curv")
lh_sulc_map_binary = np.where(lh_sulc < 0.0, 0.25, 0.6)
rh_sulc_map_binary = np.where(rh_sulc < 0.0, 0.25, 0.6)# Load the surface files for left and right hemispheres
indices_left_v1 = get_indices_roi([1, 2, 3], left_temp)
indices_right_v1 = get_indices_roi([1, 2, 3], right_temp)

# Assuming `bin` is defined and `df` contains your data
fig = plt.figure(figsize=(11, 8))
groups = ["2nd", "3rd", "preterm",
          "fullterm", 'adolescent', 'adult', "benson"]
text_label = {'preterm': 'preterm neonatal', 'fullterm': 'full-term neonatal',
              '2nd': 'tri-2nd prenatal', '3rd': 'tri-3rd prenatal', 'adolescent': 'adolescent', 'adult': 'adult', "benson": "template"}
left_label = ["a", "b", "c", "d", "e", "f", "g"]
colors = {'preterm': 'limegreen', 'fullterm': 'mediumaquamarine', '2nd': 'gold',
          '3rd': 'yellowgreen', 'adolescent': 'cadetblue', 'adult': 'steelblue', "benson":"#656364"}
vmax = 180
# Load the GIFTI file (atlas file)
atlas_lh = nib.load(f'{file_path}/hemi-L_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii')
atlas_rh = nib.load(f'{file_path}/hemi-R_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii')
template_lh = nib.load(lh)  # Load the left hemisphere surface
template_rh = nib.load(rh)
fit = pd.read_csv(f'{file_path}/fitting_linear_quad.csv')
# Extract the label data
atlas_data_l = atlas_lh.darrays[0].data
atlas_data_r = atlas_rh.darrays[0].data  # Extract the data array
# Specify the labels of interest
left_coords, left_faces = nib.freesurfer.io.read_geometry(f"{file_path}/lh.sphere")
right_coords, right_faces = nib.freesurfer.io.read_geometry(f"{file_path}/rh.sphere")
def divide_indices_by_z(coords, indices):
    z_values = coords[indices][:, 2]
    middle_z = np.median(z_values)
    below_indices = indices[z_values < middle_z]
    above_indices = indices[z_values >= middle_z]
    return below_indices, above_indices
# Compute R^2 values
def r2_score(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)  # Residual sum of squares
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)  # Total sum of squares
    return 1 - (ss_res / ss_tot)

# Assuming `left_temp` and `right_temp` contain the coordinates for left and right hemispheres
# For left hemisphere V2
below_middle_indices_left_v2, above_middle_indices_left_v2 = divide_indices_by_z(left_coords, indices_left_v2)

# For left hemisphere V3
below_middle_indices_left_v3, above_middle_indices_left_v3 = divide_indices_by_z(left_coords, indices_left_v3)

# For right hemisphere V2
below_middle_indices_right_v2, above_middle_indices_right_v2 = divide_indices_by_z(right_coords, indices_right_v2)

# For right hemisphere V3
below_middle_indices_right_v3, above_middle_indices_right_v3 = divide_indices_by_z(right_coords, indices_right_v3)

distance_left= pd.read_csv(f"{file_path}/L_polarangle_distance.csv")   
distance_right= pd.read_csv(f"{file_path}/R_polarangle_distance.csv") 

for i, group in enumerate(groups):    
    if group == "preterm":
        average_left = f"{file_path}/Averaged_younger_n_L_label-polarangle_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_younger_n_R_label-polarangle_desc-real_roi-v2th00_metric.gii"
    elif group == "fullterm":
        average_left = f"{file_path}/Averaged_older_n_L_label-polarangle_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_older_n_R_label-polarangle_desc-real_roi-v2th00_metric.gii"
    elif group == "2nd":
        average_left = f"{file_path}/Averaged_younger_fetal_L_label-polarangle_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_younger_fetal_R_label-polarangle_desc-real_roi-v2th00_metric.gii"
    elif group == "3rd":
        average_left = f"{file_path}/Averaged_older_fetal_L_label-polarangle_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_older_fetal_R_label-polarangle_desc-real_roi-v2th00_metric.gii"
    elif group == "adolescent":
        average_left = f"{file_path}/Averaged_young_L_label-polarangle_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_young_R_label-polarangle_desc-real_roi-v2th00_metric.gii"
    elif group == "adult":
        average_left = f"{file_path}/Averaged_old_L_label-polarangle_desc-real_roi-v2th00_metric.gii"
        average_right = f"{file_path}/Averaged_old_R_label-polarangle_desc-real_roi-v2th00_metric.gii"
    elif group == "benson":
        template_left = f"{file_path}/hemi-L_space-fsaverage_dens-164k_desc-angleretinotbenson2014_seg.shape.gii"
        template_right = f"{file_path}/hemi-R_space-fsaverage_dens-164k_desc-angleretinotbenson2014_seg.shape.gii"
        ccf_l = surface.load_surf_data(template_left)             
        average_left = np.zeros(ccf_l.shape)
        average_left2 = np.zeros(ccf_l.shape)
        average_left[indices_left_v1]= ccf_l[indices_left_v1].astype(np.float64)
        average_left2[indices_left]= ccf_l[indices_left].astype(np.float64)
        ccf_r = surface.load_surf_data(template_right)
        average_right = np.zeros(ccf_r.shape)
        average_right2 = np.zeros(ccf_r.shape)
        average_right[indices_right_v1]= ccf_r[indices_right_v1].astype(np.float64)
        average_right2[indices_right]= ccf_r[indices_right].astype(np.float64)
                    
    # Histogram positions
    if i==0:
        hist_L_pos_0 = [0.1, 1 +0.07 - i*0.13, 0.12, 0.05]
        hist_L_pos_1 = [0.1, 1- i*0.13, 0.12, 0.05]
        # Adjusted position for the second histogram
        hist_R_pos_0 = [0.59, 1 +0.07 - i*0.13, 0.12, 0.05]
        hist_R_pos_1 = [0.59, 1- i*0.13, 0.12, 0.05]
    else:
        hist_L_pos_0 = [0.1, 1 +0.059 - i*0.13, 0.12, 0.05]
        hist_L_pos_1 = [0.1, 1-0.005- i*0.13, 0.12, 0.05]
        # Adjusted position for the second histogram
        hist_R_pos_0 = [0.59, 1 +0.059 - i*0.13, 0.12, 0.05]
        hist_R_pos_1 = [0.59, 1-0.005- i*0.13, 0.12, 0.05]

    inner_left = [hist_L_pos_0, hist_L_pos_1]
    inner_right = [hist_R_pos_0, hist_R_pos_1]
    
    # Load polarangle values
    average_left = surface.load_surf_data(average_left)
    average_right = surface.load_surf_data(average_right)
    # Define processing function
    def process_data(ax, distance_data, indices, x_labels, y_ticks, ylim, hemi, area, group, average):
        distance = np.concatenate([distance_data.loc[distance_data["Vertex"].isin(indices[0]), "Distance_to_Border"],
                                   distance_data.loc[distance_data["Vertex"].isin(indices[1]), "Distance_to_Border"]])
        indices = np.concatenate(indices)
        if hemi == "L":
            log_values = average[indices]
        elif hemi == "R":
            log_values = average[indices]
        df = pd.DataFrame({"distance": distance, "polarangle": log_values}).sort_values(by="distance").reset_index(drop=True)
        # Bin and compute stats
        num_bins = 30
        bins = np.linspace(df["distance"].min(), df["distance"].max(), num_bins + 1)
        df["distance_bin"] = pd.cut(df["distance"], bins=bins, labels=0.5 * (bins[:-1] + bins[1:]))
        stats = df.groupby("distance_bin")["polarangle"].agg(["mean", "std"]).reset_index()
        stats["distance_bin"] = stats["distance_bin"].astype(float)        
        # Plot
        #ax.scatter(df["distance"], df["polarangle"], color="k", s=0.5, alpha=0.1, label="Raw Data")  # Scatter plot
        ax.plot(stats["distance_bin"], stats["mean"], label="Averaged Polar angle", color="k")
        ax.fill_between(stats["distance_bin"],
                         stats["mean"] - stats["std"],
                         stats["mean"] + stats["std"],
                         color="gray", alpha=0.3, label="Standard Deviation")  
        ax.set_xlim(-15, 15)
        ax.set_ylim(*ylim)
        if group == "benson":            
            ax.set_yticks([]) 
            ax.set_yticklabels([])
            if area == "dorsal":                
                ax.set_xticks([])
                ax.set_xticklabels([])
            else:
                ax.set_ylabel("polar angle", labelpad=5, loc='bottom')
                ax.set_xticks([-15, 0, 15])
                ax.set_xticklabels([-15, 0, 15])
                ax.set_xlabel("distance", labelpad=2)
        elif group =="2nd":
            ax.set_xticks([-12, 0, 12])
            ax.set_xticklabels(x_labels, fontsize=7)
            ax.xaxis.set_tick_params(pad=1)
            ax.set_yticks(y_ticks)
            ax.tick_params(axis="y")  
        else:
            ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_yticks([])
            ax.set_yticklabels([])

        x = np.linspace(-15, 15, 31)
        fit = pd.read_csv(f'{file_path}/fitting_linear_quad.csv')

        subset = fit.loc[(fit["Group"] == group) & (fit["hemi"] == hemi) & (fit["area"] == area)].reset_index(drop=True)            
        if group == "benson":                
            slope = subset["slope"].values[0]
            y = subset["y_shift"].values[0]
        else:            
            # Perform bootstrapping for the slope and y-shift
            bootstrapped_slope = []
            bootstrapped_y = []
            
            for _ in range(10000):  # Bootstrapping iterations
                slope_sample = np.random.choice(subset["slope"], size=len(subset["slope"]), replace=True)
                bootstrapped_slope.append(np.mean(slope_sample))
                y_sample = np.random.choice(subset["y_shift"], size=len(subset["y_shift"]), replace=True)
                bootstrapped_y.append(np.mean(y_sample))
            
            # Calculate the mean of the bootstrapped samples
            slope = np.mean(bootstrapped_slope)
            y = np.mean(bootstrapped_y)
        
        # Restrict data to distance between -15 and 15
        valid_indices = (stats["distance_bin"] >= -15) & (stats["distance_bin"] <= 15)
        filtered_distances = stats["distance_bin"][valid_indices]
        filtered_means = stats["mean"][valid_indices]

        # Compute linear fit only for the restricted range
        linear_values = np.polyval([slope, y], abs(filtered_distances))

        # Compute R² only for the restricted range
        r2_linear = r2_score(filtered_means, linear_values)        

        #print(linear_values, stats["mean"])
        # Calculate the linear fit values using the mean slope and y-shift
        #recon_values = np.polyval(linear_fit, abs(x))
        F_stat, p_val = r_squared_significance(r2_linear, len(linear_values), 1)
        print(f"F-statistic: {F_stat:.4f}")
        print(f"P-value: {p_val:.4f}")
        
        # Plot the linear fit on the corresponding subplot
        #ax.plot(stats["distance_bin"], linear_values, color="green", linestyle="--")
        ax.plot(filtered_distances, linear_values, label="Reconstructed", color="k", linestyle="--")
        print(f"Group: {group}, Hemi: {hemi}, Area: {area}, R²: {r2_linear}, slpe: {slope}")
        if area == "dorsal":
            if p_val <0.01:
                ax.text(-5, 160,  # Adjusted position (closer to bottom-right)
                f"R²={r2_linear:.2f}*", 
                fontsize=7)
            else:
                ax.text(-5, 160,  # Adjusted position (closer to bottom-right)
                f"R²={r2_linear:.2f}", 
                fontsize=7)
        else:
            if p_val <0.01:
                ax.text(-5, 5,  # Adjusted position (closer to bottom-right)
                f"R²={r2_linear:.2f}*", 
                fontsize=7)
            else:
                ax.text(-5, 5,  # Adjusted position (closer to bottom-right)
                f"R²={r2_linear:.2f}", 
                fontsize=7)
        
    process_data(fig.add_axes(inner_left[0]), distance_left, [above_middle_indices_left_v2, above_middle_indices_left_v3], ['V3d', 0, 'V2d'], [90, 180], (90, 180), "L", "dorsal", group, average_left)
    process_data(fig.add_axes(inner_left[1]), distance_left, [below_middle_indices_left_v2, below_middle_indices_left_v3], ['V3v', 0, 'V2v'], [0, 90], (0, 90), "L", "ventral", group, average_left)
    process_data(fig.add_axes(inner_right[0]), distance_right, [above_middle_indices_right_v2, above_middle_indices_right_v3], ["V3d", 0, "V2d"], [90, 180], (90, 180), "R", "dorsal", group, average_right)
    process_data(fig.add_axes(inner_right[1]), distance_right, [below_middle_indices_right_v2, below_middle_indices_right_v3], ['V3v', 0, "V2v"], [0, 90], (0, 90), "R", "ventral", group, average_right)

    fig.text(0.395, 1 - (i * 0.13) + 0.11,
             text_label[group], ha='center', va='center', fontsize=14)
    #plt.show()
    # Create 3D axes for brain surface images and adjust positions to overlap by one-third
    ax1 = fig.add_axes([0.205, 1 - i*0.13-0.01, 0.12, 0.11], projection='3d')
    ax2 = fig.add_axes([0.29, 1 - i*0.13, 0.12, 0.1], projection='3d')
    ax3 = fig.add_axes([0.375, 1 - i*0.13, 0.12, 0.1], projection='3d')
    ax4 = fig.add_axes([0.455, 1 - i*0.13-0.01, 0.12, 0.11], projection='3d')
    
    
    log_vmin = np.log10(1)
    log_vmax = np.log10(vmax)

    plotting.plot_surf_stat_map(
        lh, average_left, bg_map=lh_sulc_map_binary,
        hemi='left', view=(0, 290), axes=ax1, darkness=1.0,
        colorbar=False, cmap="jet", vmax=vmax, vmin=0,
        threshold=0.00001)    
    
    plotting.plot_surf_contours(
        surf_mesh=lh, roi_map=atlas_data_l, hemi='left',
        view=(0, 290), axes=ax1, levels=[2, 1, 3], colors=['k', 'k', 'k'],
    )
    display=plotting.plot_surf_stat_map(
        lh, average_left, bg_map=lh_sulc_map_binary,
        hemi='left', view=(-10, -70.0), axes=ax2, darkness=1.0,
        colorbar=False, cmap="jet", vmax=vmax, vmin=0,
        threshold=0.00001)    
    plotting.plot_surf_contours(
        surf_mesh=lh, roi_map=atlas_data_l, hemi='left',
        view=(-10, -70), axes=ax2, levels=[2, 1, 3], colors=['k', 'k', 'k'],
    )

    plotting.plot_surf_stat_map(
        rh, average_right, bg_map=rh_sulc_map_binary,
        hemi='right', view=(-10, -110.0), axes=ax3, darkness=1.0,
        colorbar=False, cmap="jet", vmax=vmax, vmin=0,
        threshold=0.00001)
    plotting.plot_surf_contours(
        surf_mesh=rh, roi_map=atlas_data_r, hemi='right',
        view=(-10, -110), axes=ax3, levels=[2, 1, 3], colors=['k', 'k', 'k'],
    )
    display1=plotting.plot_surf_stat_map(
        rh, average_right, bg_map=rh_sulc_map_binary,
        hemi='right', view=(0, 250), axes=ax4, darkness=1.0,
        colorbar=False, cmap="jet", vmax=vmax, vmin=0,
        threshold=0.00001)
    plotting.plot_surf_contours(
        surf_mesh=rh, roi_map=atlas_data_r, hemi='right',
        view=(0, 250), axes=ax4, levels=[2, 1, 3], colors=['k', 'k', 'k'],
    )
    
    # Remove axis ticks and labels for brain surface plots
    ax1.set_axis_off()
    ax2.set_axis_off()
    ax3.set_axis_off()
    ax4.set_axis_off()

    ax1.set_xlim([-50, 50])
    ax1.set_ylim([-50, 50])
    ax1.set_zlim([-50, 50])
    ax2.set_xlim([-50, 10])
    ax2.set_ylim([-35, 35])
    ax2.set_zlim([-30, 35])
    ax3.set_xlim([-10, 50])
    ax3.set_ylim([-35, 35])
    ax3.set_zlim([-30, 35])
    ax4.set_xlim([-50, 50])
    ax4.set_ylim([-50, 50])
    ax4.set_zlim([-50, 50])
        

# Adjust layout to leave space for the colorbar at the bottom
plt.subplots_adjust(bottom=0.15, hspace=0.5)
#for i in range(7):
#        fig.text(0.085, 1 - (i * 0.13)+0.108,
#                left_label[i], ha='center', va='center', fontsize=12, weight='bold')
        
fig.text(0.24, 1.09, "L",  ha='center', va='center', fontsize=12)
fig.text(0.54, 1.09, "R",  ha='center', va='center', fontsize=12)

fig.text(0.392, 1.09, "———V3———",  ha='center', va='center', fontsize=5,)
fig.text(0.392, 1.075, "———V2———",  ha='center', va='center', fontsize=5,)
fig.text(0.392, 0.31, "———V3———",  ha='center', va='center', fontsize=5,)
fig.text(0.392, 0.295, "———V2———",  ha='center', va='center', fontsize=5,)
fig.text(0.392, 0.27, "———V1———",  ha='center', va='center', fontsize=5, )
# Add colorbar
# Adjust the position and height of the colorbar axis
cbar_ax = fig.add_axes([0.24, 0.18, 0.3, 0.02])
cbar = plt.colorbar(cm.ScalarMappable(cmap='jet', norm=plt.Normalize(
    vmin=0, vmax=vmax)), cax=cbar_ax, orientation='horizontal', fraction=.1)
#cbar.set_label('polar angle')
# Set the ticks of the colorbar to include the min and max values
cbar.set_ticks([0, vmax/2, vmax])

# Set the tick labels to show min and max values (you can format them as needed)
cbar.set_ticklabels([ 0, "polar angle",f'{vmax:.0f}'])

#plt.show()
# Save the figure
fig.savefig(f"{file_path}/Figure2.png",
            dpi=300, bbox_inches='tight')
