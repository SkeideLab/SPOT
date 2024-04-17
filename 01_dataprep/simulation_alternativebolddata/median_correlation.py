import numpy as np
from nilearn import plotting, surface

# Load your fMRI data (example data, replace with your own data)
# Here, assume data is a 4D numpy array with dimensions (time, x, y, z)
FUNC_PATH = ("/data/p_02915/dhcp_derivatives_SPOT/dhcp_surface/sub-CC00056XX07/ses-10700/func/sub-CC00056XX07_ses-10700_hemi-L_mesh-native_bold.func.gii")

surface_data = surface.load_surf_data(FUNC_PATH)
num_vertices, num_time_points = surface_data.shape
surface_data = surface_data.astype(float)

surface_data_2d = surface_data

# Check for NaN and infinite values
nan_mask = np.isnan(surface_data)
inf_mask = np.isfinite(surface_data)

# Count NaN and infinite values
num_nan = np.sum(nan_mask)
num_inf = np.sum(~inf_mask)

print("Number of NaN values:", num_nan)
print("Number of infinite values:", num_inf)

print("Type of surface_data:", type(surface_data))
print("Shape of surface_data:", surface_data.shape)
# Calculate correlation coefficient voxel by voxel
correlation_map = np.corrcoef(surface_data)
# masked_correlation_map = np.ma.masked_array(correlation_map, mask=np.isnan(correlation_map))
nan_mask = np.isnan(correlation_map)

# Invert the mask to get a mask for non-NaN values
non_nan_mask = ~nan_mask

# Use the mask to select only the non-NaN values from the array
masked_correlation_map = correlation_map[non_nan_mask]
print(masked_correlation_map)
# Find the median value of all correlation coefficients
median_correlation = np.median(masked_correlation_map)
print("Median correlation coefficient:", median_correlation)

