#!/bin/bash

AFNI SurfFWHM \
    -input /data/p_02495/dhcp_derivatives/dhcp_surface/sub-CC00058XX09/ses-11300/func/sub-CC00058XX09_ses-11300_hemi-L_space-T2w_desc-simulatedsmoothed_bold_tmp.func.gii \
    -i_gii /data/p_02495/dhcp_derivatives/dhcp_anat_pipeline/sub-CC00058XX09/ses-11300/anat/sub-CC00058XX09_ses-11300_hemi-L_space-T2w_pial.surf.gii \
    -i_gii /data/p_02495/dhcp_derivatives/dhcp_anat_pipeline/sub-CC00058XX09/ses-11300/anat/sub-CC00058XX09_ses-11300_hemi-L_space-T2w_sphere.surf.gii \
    -detrend 0 \
    -hood -1

# AFNI SurfFWHM \
#     -input /data/p_02495/dhcp_derivatives/dhcp_surface/sub-CC00058XX09/ses-11300/func/sub-CC00058XX09_ses-11300_hemi-L_space-T2w_desc-simulated_bold_tmp.func.gii \
#     -i_gii /data/p_02495/dhcp_derivatives/dhcp_anat_pipeline/sub-CC00058XX09/ses-11300/anat/sub-CC00058XX09_ses-11300_hemi-L_space-T2w_pial.surf.gii
