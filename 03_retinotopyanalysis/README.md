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
- `calculate_tSNR.py`: calculate tSNR
- `calculate_slice_SNR.py`: calculate sSNR

# Statistical analysis across the groups (kruskal (6-class) and mann-whitney-u tests (group pairs))

- `kruskal_mannwu_dhcp_fetal_hcp.py`: group comparison of r, sigma, eccentricity, and polar angle
- `kruskal_mannwu_dhcp_fetal_hcp_FD.py`: group comparison of dvars
- `kruskal_mannwu_dhcp_fetal_hcp_gradient.py`: group comparison of gradient
- `kruskal_mannwu_dhcp_fetal_hcp_tSNR.py`: group comparison of tSNR and sSNR

# Statistical analysis between real and simulated (Wilcoxon signed-rank test)

- `statistics_real_simul_corr.py`: comparison of r between real and simulated data (each hemisphere)
- `statistics_real_simul_msd.py`: comparison of msd between real and simulated data (each hemisphere)
- `statistics_real_simul_msd_wb`: comparison of averaged msd in both hemispheres
