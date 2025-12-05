#!/bin/bash
set -u -x -e

# GOAL: registration from fsaverage to hcp fs_lr to dhcpSym40 to native
# Use this registration to resample the wang template to individual surfaces

path_script=$(dirname $0)
sub=$1
ses=$2
path_bids_data=$3
path_output_data=$4
path_HCPtemplates_standardmeshatlases=$5
path_surfacetemplate=$6
path_fsaverage=$7

# prepare output directory for subject in T2w space
path_out_native=$path_output_data/sub-$sub/ses-$ses/anat
mkdir -p $path_out_native

for hemi in left right; do
    if [ $hemi = "left" ]; then
        hemi_upper="L"
    elif [ $hemi = "right" ]; then
        hemi_upper="R"
    fi

    # registrations necessary to make the transformation from fsaverage to native
    # will be created if they don't exist already
    registration_from_dhcpSym_to_fsaverage="$path_script/standard_registrations/hemi-${hemi_upper}_from-dhcpSym_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii"
    registration_from_native_to_fsaverage="$path_output_data/sub-$sub/ses-$ses/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii"

    # output of previous step, see alignment folder
    # check if this exists, terminate otherwise
    # should be changed for fetal
    registration_from_native_to_dhcpSym="$path_output_data/sub-$sub/ses-$ses/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-dhcpSym40_dens-32k_mode-sphere.reg40.surf.gii"
    if [ ! -f $registration_from_native_to_dhcpSym ]; then
        echo "Registration from native to dhcpSym is missing. Please execute alignment folder main script first."
        echo "Terminating..."
        exit 1
    fi

    # templates in fsaverage space
    # come with the repo, for more info see templates_retinotopy README
    path_templates_fsaverage="$path_script/templates_retinotopy/hemi-${hemi_upper}_space-fsaverage_dens-164k"
    sigma_template="${path_templates_fsaverage}_desc-sigmaretinotbenson2014_seg.shape.gii"

    # output
    path_out_native_hemi="$path_out_native/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_mesh-native_dens-native"
    sigma_native="${path_out_native_hemi}_desc-sigmaretinotbenson2014_seg.shape.gii"
    if [ ! -f $sigma_native ]; then
        wb_command -metric-resample \
            $sigma_template \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
            $registration_from_native_to_fsaverage \
            ADAP_BARY_AREA \
            $sigma_native -area-surfs \
            $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
            $path_bids_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_midthickness.surf.gii
    fi

done
