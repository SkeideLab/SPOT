#!/bin/bash
# transforms the Wang template to individual subject surface space
set -x -u -e

########################INPUTS######################################
#---------------SUBJECT INFO---------------------------------
subject_label="CC00066XX09"
session_label="19200"
subject_age_weeks="42"
#---------------PATHS TO DATA--------------------------------------------
path_derivatives="/data/p_02495/dhcp_derivatives"
path_anat_data="$path_derivatives/dhcp_anat_pipeline"
path_func_data="$path_derivatives/dhcp_fmri_pipeline"
path_output_data="$path_derivatives/dhcp_surface"
#---------------PATHS TO TEMPLATES---------------------------------------
file_volume_template_40wks="/data/p_02495/templates/template_augmentedvolumetricatlas_dhcp/atlas/T2/template-40.nii.gz"
name_volume_template_40wks="dhcp40wk"
path_surface_template="/data/p_02495/templates/template_corticalsurfaceneonatessym_williams2023_dhcp/dhcpSym_template"
name_surface_template="dhcpSym"
path_HCPtemplates_standardmeshatlases="/data/u_kieslinger_software/code/HCPpipelines/global/templates/standard_mesh_atlases"
path_fsaverage="/data/p_02495/templates/template_fsaverage/fsaverage"
#--------------PATHS TO SOFTWARE--------------------------------------------------
path_newmsm="/data/u_kieslinger_software/fsldevdir/bin/newmsm"
path_wbcommand="/bin/wb_command"
path_mirtk="/afs/cbs.mpg.de/software/mirtk/0.20231123/debian-bullseye-amd64/bin/mirtk"
#---------------------------------------------------------------------------------
#######################################################################################
path_script=$(dirname $0)
logfile=$path_output_data/sub-${subject_label}/ses-${session_label}/logs/surfaceprep_$(date +%Y-%m-%dT%H%M).log
mkdir -p $(dirname $logfile)
echo "Redirecting outputs to logfile $logfile..."
{
    ###########################REGISTRATION##############################################
    # Create registration from individual space to surface template space
    MIRTK $path_script/alignment/align_to_template_2nd_release.sh \
        $path_anat_data \
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

    ###############################RESAMPLING RETINOTOPY TEMPLATE#############################
    $path_script/retinotopy/warp_templates_to_native.sh \
        $subject_label \
        $session_label \
        $path_anat_data \
        $path_output_data \
        $path_HCPtemplates_standardmeshatlases \
        $path_surface_template \
        $path_fsaverage

    ##############################DEFINING ROIS###########################################
    python $path_script/retinotopy/divide_rois_into_subregions.py \
        -o $path_output_data \
        -sub $subject_label \
        -ses $session_label

    ###############################PROJECTING FUNCTIONAL DATA TO SURFACE##########################
    FSL $path_script/surfaceprojection/project_bold_to_surface.sh \
        $subject_label \
        $session_label \
        $path_wbcommand \
        $path_anat_data \
        $path_func_data \
        $path_output_data

    ########################SIMULATING ALTERNATIVE DATA###################################
    python $path_script/simulation_alternativemodel/alternativemodel.py \
        -d $path_derivatives \
        -sub $subject_label \
        -ses $session_label

} >>$logfile 2>&1
