#!/bin/bash
# transforms the Wang template to individual subject surface space
set -x -u -e
#---------------SUBJECT INFO---------------------------------
subject_label="CC00062XX05"
session_label="13801"
subject_age_weeks="41"
#---------------PATHS TO DATA--------------------------------------------
path_bids_data=/data/p_02495/dhcp_derivatives/dhcp_anat_pipeline
path_output_data=/data/p_02495/dhcp_derivatives/dhcp_surface
#---------------PATHS TO TEMPLATES---------------------------------------
file_volume_template_40wks=/data/p_02495/templates/template_augmentedvolumetricatlas_dhcp/atlas/T2/template-40.nii.gz
name_volume_template_40wks="dhcp40wk"
path_surface_template=/data/p_02495/templates/template_corticalsurfaceneonatessym_williams2023_dhcp/dhcpSym_template
name_surface_template="dhcpSym"
#--------------PATHS TO SOFTWARE--------------------------------------------------
path_newmsm="/data/u_kieslinger_software/fsldevdir/bin/newmsm"
path_wbcommand="/bin/wb_command"
path_mirtk="/afs/cbs.mpg.de/software/mirtk/0.20231123/debian-bullseye-amd64/bin/mirtk"
#---------------------------------------------------------------------------------

path_script=$(dirname $0)

# Create registration from individual space to surface template space
$path_script/alignment/align_to_template_2nd_release.sh \
    $path_bids_data \
    $subject_label \
    $session_label \
    $subject_age_weeks \
    $file_volume_template_40wks \
    $name_volume_template_40wks \
    $path_surface_template \
    $name_surface_template \
    $path_script/alignment/rotational_transforms/week40_toFS_LR_rot.%hemi%.txt \
    $path_output_data \
    $path_script/alignment/configs/config_subject_to_40_week_template \
    $path_newmsm \
    $path_wbcommand \
    $path_mirtk

# $(make_fsavg_to_native_warp.sh)
# $(hcp_surface.sh)
