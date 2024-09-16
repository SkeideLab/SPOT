#!/bin/bash
# script for HCP-D data release
# script to seprate cifti file to gifti
set -x -u -e

path_script=$(dirname $0)
subjid=$1
path_anat_data=$2
path_func_data=$3
outdir=$4
WB_BIN=$5

########## DEFINE PATHS TO VARIABLES ##########
# TODO define mesh and space suffixes for all surfaces
#inputs
native_func_mesh=${path_func_data}/rfMRI_REST_Atlas_MSMAll_hp0_clean_vn.dscalar.nii
native_func_data=${path_func_data}/rfMRI_REST_Atlas_MSMAll_hp0_clean.dtseries.nii

#outputs
sub_output_dir=${outdir}/${subjid}
mkdir -p $sub_output_dir/func $sub_output_dir/anat

outname_left=$sub_output_dir/func/${subjid}_hemi-L-Atlas-MSMAll_hp0_clean_bold.func.gii
outname_right=$sub_output_dir/func/${subjid}_hemi-R-Atlas-MSMAll_hp0_clean_bold.func.gii
outname_left_mesh=$sub_output_dir/func/${subjid}_hemi-L-Atlas-MSMAll_hp0_clean_bold.shape.gii
outname_right_mesh=$sub_output_dir/func/${subjid}_hemi-R-Atlas-MSMAll_hp0_clean_bold.shape.gii


# Separate one cifti file into two gifti files
${WB_BIN} -cifti-separate \
    $native_func_data \
    COLUMN -metric \
    CORTEX_LEFT $outname_left \

${WB_BIN} -cifti-separate \
    $native_func_data \
    COLUMN -metric \
    CORTEX_RIGHT $outname_right

${WB_BIN} -cifti-separate \
    $native_func_mesh \
    COLUMN -metric \
    CORTEX_LEFT $outname_left_mesh \

${WB_BIN} -cifti-separate \
    $native_func_mesh \
    COLUMN -metric \
    CORTEX_RIGHT $outname_right_mesh