#!/bin/bash
set -u -x -e

sub=$1
ses=$2
hemi=$3
path_anat_data=$4
path_output_data=$5
path_HCPtemplates_standardmeshatlases=$6
path_fsaverage=$7
path_wbcommand=$8

#FIXME
# WARNING It does not make sense to transform vertices with
# indices to V1 to fsaverage without transforming the indices themselves.
# This needs to be reworked to transform the retinotopy values themselves after assignment.

file_ccf_v0="/data/p_02495/dhcp_derivatives/ccfmodel/sub-$sub/ses-$ses/sub-${sub}_ses-${ses}_hemi-${hemi}_desc-test_v0i.gii"
sphere_native="/data/p_02495/dhcp_derivatives/dhcp_anat_pipeline/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_space-T2w_sphere.surf.gii"
registration=$path_output_data/sub-$sub/ses-$ses/space-dhcpSym/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii

mkdir -p $path_output_data/sub-$sub/ses-$ses/space-T2w/surface_transforms $path_output_data/sub-$sub/ses-$ses/space-fsaverage/anat/
#TODO reverse registration
$path_wbcommand -surface-sphere-project-unproject \
    $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi}.164k_fsavg_${hemi}.surf.gii \
    $registration `# sphere-project-to: sphere that aligns with sphere-in` \
    $path_anat_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_space-T2w_sphere.surf.gii `# sphere-unproject-from: sphere-project-to deformed to desired output space: native (not rotated)` \
    $path_output_data/sub-$sub/ses-$ses/space-T2w/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi}_from-fsaverage_to-nativeT2w_dens-164k_mode-sphere_reg.surf.gii `# sphere-out: output sphere`

$path_wbcommand -label-resample \
    $file_ccf_v0 `# label-in: ` \
    $sphere_native`# current sphere: replace with rotated sphere??? ` \
    $path_output_data/sub-$sub/ses-$ses/space-T2w/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi}_from-fsaverage_to-nativeT2w_dens-164k_mode-sphere_reg.surf.gii`# new sphere: sphere in register with current sphere and desired mesh` \
    ADAP_BARY_AREA \
    $path_output_data/sub-$sub/ses-$ses/space-fsaverage/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_mesh-fsaverage_dens-164k_desc-v2_label-v0i_dparc.label.gii \
    -area-surfs \
    $path_anat_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_space-T2w_midthickness.surf.gii \
    $path_fsaverage/surf/$(echo $hemi | tr '[:upper:]' '[:lower:]')h_midthickness_surf.gii
