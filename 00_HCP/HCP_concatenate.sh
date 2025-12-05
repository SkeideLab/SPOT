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
# Loop through the resting-state fMRI parameters
for param in rfMRI_REST1_AP rfMRI_REST1_PA rfMRI_REST2_AP rfMRI_REST2_PA; do
    native_func_mesh=${path_func_data}/${param}/${param}_Atlas_MSMAll_hp0_vn.dscalar.nii
    native_func_data=${path_func_data}/${param}/${param}_Atlas_MSMAll_hp0_clean.dtseries.nii

    # Check if both native_func_mesh and native_func_data exist
    if [[ -f "$native_func_mesh" && -f "$native_func_data" ]]; then

        # Outputs
        sub_output_dir=${outdir}/${subjid}
        mkdir -p $sub_output_dir/func $sub_output_dir/anat

        outname_left=$sub_output_dir/func/${subjid}_${param}_hemi-L-Atlas-MSMAll_hp0_clean_bold.func.gii
        outname_right=$sub_output_dir/func/${subjid}_${param}_hemi-R-Atlas-MSMAll_hp0_clean_bold.func.gii
        outname_left_mesh=$sub_output_dir/func/${subjid}_${param}_hemi-L-Atlas-MSMAll_hp0_clean_bold.shape.gii
        outname_right_mesh=$sub_output_dir/func/${subjid}_${param}_hemi-R-Atlas-MSMAll_hp0_clean_bold.shape.gii

        # Separate one cifti file into two gifti files
        ${WB_BIN} -cifti-separate $native_func_data COLUMN -metric CORTEX_LEFT $outname_left
        ${WB_BIN} -cifti-separate $native_func_data COLUMN -metric CORTEX_RIGHT $outname_right
        ${WB_BIN} -cifti-separate $native_func_mesh COLUMN -metric CORTEX_LEFT $outname_left_mesh
        ${WB_BIN} -cifti-separate $native_func_mesh COLUMN -metric CORTEX_RIGHT $outname_right_mesh

    else
        echo "File $native_func_mesh or $native_func_data does not exist, skipping $param."
    fi
done

# Now, concatenate the resulting GIFTI files for each hemisphere

# Define the output for concatenated files
concatenated_left_output=$sub_output_dir/func/${subjid}_hemi-L_concatenated.func.gii
concatenated_right_output=$sub_output_dir/func/${subjid}_hemi-R_concatenated.func.gii

# Initialize the command with the output file name for left and right hemisphere
left_merge_cmd="${WB_BIN} -metric-merge $concatenated_left_output"
right_merge_cmd="${WB_BIN} -metric-merge $concatenated_right_output"

# Loop again to add only the existing files to the merge commands
for param in rfMRI_REST1_AP rfMRI_REST1_PA rfMRI_REST2_AP rfMRI_REST2_PA; do
    left_file=${outdir}/${subjid}/func/${subjid}_${param}_hemi-L-Atlas-MSMAll_hp0_clean_bold.func.gii
    right_file=${outdir}/${subjid}/func/${subjid}_${param}_hemi-R-Atlas-MSMAll_hp0_clean_bold.func.gii

    # Check if left and right hemisphere files exist before adding to merge command
    if [[ -f "$left_file" ]]; then
        left_merge_cmd="$left_merge_cmd -metric $left_file"
    fi
    if [[ -f "$right_file" ]]; then
        right_merge_cmd="$right_merge_cmd -metric $right_file"
    fi
done

# Run the merge commands if there are any files to merge
if [[ "$left_merge_cmd" != "${WB_BIN} -metric-merge $concatenated_left_output" ]]; then
    eval $left_merge_cmd
else
    echo "No left hemisphere files to concatenate."
fi

if [[ "$right_merge_cmd" != "${WB_BIN} -metric-merge $concatenated_right_output" ]]; then
    eval $right_merge_cmd
else
    echo "No right hemisphere files to concatenate."
fi

# Calculate the mean across the concatenated GIFTI metrics for Left Hemisphere
mean_left_output=$sub_output_dir/func/${subjid}_hemi-L_mean.func.gii
if [[ -f "$concatenated_left_output" ]]; then
    ${WB_BIN} -metric-reduce $concatenated_left_output MEAN $mean_left_output
else
    echo "No concatenated left hemisphere file to compute mean."
fi

# Calculate the mean across the concatenated GIFTI metrics for Right Hemisphere
mean_right_output=$sub_output_dir/func/${subjid}_hemi-R_mean.func.gii
if [[ -f "$concatenated_right_output" ]]; then
    ${WB_BIN} -metric-reduce $concatenated_right_output MEAN $mean_right_output
else
    echo "No concatenated right hemisphere file to compute mean."
fi
