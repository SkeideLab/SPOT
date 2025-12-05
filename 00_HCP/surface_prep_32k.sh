#!/bin/bash
# transforms the Wang template to individual subject surface space
set -x -u -e
echo "$0" "$@" # print function call

subject_label=$1
path_derivatives=$2
path_anat_data=$3
path_func_data=$4
path_output_data=$5
path_HCPtemplates_standardmeshatlases=${6}
path_fsaverage=${7}
path_newmsm=${8}
path_wbcommand=${9}
path_mirtk=${10}

#######################################################################################
path_script=$(dirname $0)
logfile=$path_output_data/${subject_label}/logs/surfaceprep_$(date +%Y-%m-%dT%H%M).log
mkdir -p $(dirname $logfile)
echo "Redirecting outputs to logfile $logfile..."
{
    MIRTK $path_script/align_native_fsaverage.sh \
         $path_anat_data \
         $subject_label \
         $path_output_data \
         $path_script/config_subject_to_fsaverage \
         $path_newmsm \
         $path_wbcommand \
         $path_mirtk

    ###############################RESAMPLING RETINOTOPY TEMPLATE#############################
    $path_script/HCP_warp_templates_to_native_32k.sh \
        $subject_label \
        $path_anat_data \
        $path_output_data \
        $path_HCPtemplates_standardmeshatlases \
        $path_fsaverage


} >>$logfile 2>&1
