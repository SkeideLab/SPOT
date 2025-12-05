#!/bin/bash
# transforms the Wang template to individual subject surface space
set -x -u -e
echo "$0" "$@" # print function call

subject_label=$1
session_label=$2
subject_age_weeks=$3
path_derivatives=$4
path_anat_data=$5
path_func_data=$6
path_output_data=$7
file_volume_template_40wks=$8
name_volume_template_40wks=$9
path_surface_template=${10}
name_surface_template=${11}
path_HCPtemplates_standardmeshatlases=${12}
path_fsaverage=${13}
path_newmsm=${14}
path_wbcommand=${15}
path_mirtk=${16}

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

    # ###############################PROJECTING FUNCTIONAL DATA TO SURFACE##########################
     FSL $path_script/surfaceprojection/project_bold_to_surface.sh \
         $subject_label \
         $session_label \
         $path_wbcommand \
         $path_anat_data \
         $path_func_data \
         $path_output_data

    ########################SIMULATING ALTERNATIVE DATA###################################
    /data/u_yoos_software/miniforge3/envs/SPOT/bin/python $path_script/simulation_alternativebolddata/simulated_model_smooth_surf.py \
        -da $path_anat_data \
        -d $path_derivatives \
        -sub $subject_label \
        -ses $session_label

} >>$logfile 2>&1
