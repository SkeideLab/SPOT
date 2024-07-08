# 00_HCP

00_HCP includes all codes for the retinotopy analysis (dataprep, ccfanalysis, and retinotopyanalysis) for Human connectome project dataset (adolescent and adult)
Note, preprocessed HCP dataset is already in the fs_LR space so you only need a transformation from fs_LR to fsaverage.

# Folders

`template´ : This folder has retinotopy templates (Benson's eccentricity and polarangle and Wang's template) on the fs_LR space and the transformation from fsaverage to fs_LR space.
- 'ccf_model': This folder has subcodes are included in 'HCP_ccf_model', it is same with 02_ccfanalysis>ccf_model

# Sequence of code
1. 'HCP_run_all' : Contains following codes 
   - 'HCP_file_prep': convert cifti to gifti
   - 'HCP_warp_templates to native': transformation from fsaverage to fsLR
   - 'HCP_simulation_model': generate simulated data on V1 and V2
   - 'HCP_ccf_model': run cortical connective field modeling
   - 'HCP_retinotopy_analysis': project V1 retinotopy template values to V2 based on cortico-connective field results
   - 'HCP_native_to fsaverage': transform eccentricity and polar angle maps from fsLR space to fsaverage space

2. 'HCP_parameter_to_fsaverage': transform r, r^2, sigma, v0 from fsLR space to fsaverage space
3. 'averaged_retinotopy': average eccentricity and polar angle maps
4. 'HCP_median_correlation': calculate median correlation
5. 'HCP_MSD': calculate mean square difference between real and simulated data
