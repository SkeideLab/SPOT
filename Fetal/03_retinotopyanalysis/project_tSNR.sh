#!/bin/bash
# project ccf results from native space to fsaverage space
set -u -x -e

sub=$1
ses=$2
path_anat_data=$3
path_output_dir=$4
path_HCPtemplates_standardmeshatlases=$5
path_fsaverage=$6
path_wbcommand=$7

threshold="00"

for hemi in L R; do
    if [ $hemi = "L" ]; then
        hemi_down="left"
    elif [ $hemi = "R" ]; then
        hemi_down="right"
    fi
 
    sphere_native="${path_anat_data}/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi_down}_sphere.surf.gii"
    registration_from_native_to_fsaverage="${path_output_dir}/dhcp_surface/sub-$sub/ses-$ses/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii"
    registration_from_fsaverage_to_native="${path_output_dir}/dhcp_surface/sub-$sub/ses-$ses/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi}_from-fsaverage_to-nativeT2w_dens-164k_mode-sphere_reg.surf.gii"

    res_file_native="/data/p_02915/dhcp_derivatives_SPOT/statistics/${sub}_${ses}_${hemi}_sliceSNR.gii"
    res_file_fsaverage="/data/p_02915/dhcp_derivatives_SPOT/statistics/${sub}_${ses}_${hemi}_sliceSNR_fsaverage.gii"

    $path_wbcommand -metric-resample \
        $res_file_native `# label-in: ` \
        $sphere_native`# current sphere: ` \
        $registration_from_fsaverage_to_native`# new sphere: sphere in register with current sphere and desired mesh` \
        ADAP_BARY_AREA \
        $res_file_fsaverage \
        -area-surfs \
        $path_anat_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi_down}_midthickness.surf.gii \
        $path_fsaverage/surf/$(echo $hemi | tr '[:upper:]' '[:lower:]')h_midthickness_surf.gii

done
