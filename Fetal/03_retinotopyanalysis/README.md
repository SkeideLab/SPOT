## 03_retinotopyanalysis
Retinotopy analysis consists of projecting V1 retinotopy template values to V2 based on cortico-connective field result, preparation results for statictical analysis, and codes for statistical analysis.

# Retinotopy analysis of dhcp dataset

- `analyse_retinotopy.py`: project V1 retinotopy template values to V2 based on cortico-connective field result
- `averaged_retinotopy.py`: averaged eccentricity and polar angle maps (fsaverage space)
- `project_ccf_to_fsaverage.sh`: transform eccentricity and polar angle values from native space to fsaverage space
- `project_correlation_to_fsaverage.sh`: transform r, r^2, sigma, v0 from native space to fsaverage space
- `msd_without_threshold.py`: calculate mean square difference

# Calculate parameters for group comparison

- `calculate_dvars.py`: calculate dvars 
- `calculate_gradient.py`: calculate gradient
