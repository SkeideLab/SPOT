Codes in this folder is to analyze Human connectome project dataset (adolescent and adult)
1. HCP_run_all
   - HCP_file_prep: convert cifti to gifti
   - HCP_warp_templates to native: transformation from fsaverage to fsLR
   - HCP_simulation_model: generate simulated data on V1 and V2
   - HCP_ccf_model: run cortical connective field modeling
   - HCP_retinotopy_analysis: project V1 retinotopy template values to V2 based on cortico-connective field results
   - HCP_native_to fsaverage: transform eccentricity and polar angle maps from fsLR space to fsaverage space
3. HCP_parameter_to_fsaverage: transform r, r^2, sigma, v0 from fsLR space to fsaverage space
4. averaged_retinotopy: average eccentricity and polar angle maps
5. HCP_median_correlation: calculate median correlation
6. HCP_MSD: calculate mean square difference between real and simulated data

