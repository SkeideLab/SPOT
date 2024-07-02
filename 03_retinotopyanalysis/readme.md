# Retinotopy analysis of dhcp dataset

- analyse_retinotopy: project V1 retinotopy template values to V2 based on cortico-connective field result
- averaged_retinotopy: averaged eccentricity and polar angle maps (fsaverage space)
- project_ccf_to_fsaverage: transform eccentricity and polar angle values from native space to fsaverage space
- project_correlation_to_fsaverage: transform r, r^2, sigma, v0 from native space to fsaverage space
- msd_without_threshold: calculate mean square difference

# Calculate parameters for group comparison

- calculate_dvars
- calculate_gradient

# Statistical analysis across the groups (kruskal (6-class) and mann-whitney-u tests (group pairs))

- kruskal_mannwu_dhcp_fetal_hcp: group comparison of r, sigma, eccentricity, and polar angle
- kruskal_mannwu_dhcp_fetal_hcp_FD: group comparison of dvars
- kruskal_mannwu_dhcp_fetal_hcp_gradient: group comparison of gradient

# Statistical analysis between real and simulated

- statistics_real_simul_corr
- statistics_real_simul_msd
- statistics_real_simul_msd_wb: comparison of averaged msd in both hemispheres.
