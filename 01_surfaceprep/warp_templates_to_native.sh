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

for hemi in left right; do
    if [ $hemi = "left" ]; then
        hemi_upper="L"
    elif [ $hemi = "right" ]; then
        hemi_upper="R"
    fi

    visareas_wang_template=$path_script/templates_retinotopy/hemi-${hemi_upper}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii
    visareas_benson_template=$path_script/templates_retinotopy/hemi-${hemi_upper}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii
    ecc_template=$path_script/templates_retinotopy/hemi-${hemi_upper}_space-fsaverage_dens-164k_desc-eccentretinotbenson2014_seg.shape.gii
    angl_template=$path_script/templates_retinotopy/hemi-${hemi_upper}_space-fsaverage_dens-164k_desc-angleretinotbenson2014_seg.shape.gii

    # 1. Combine registration from dhcpSym to fs_lr to fsaverage
    registration_from_dhcpSym_to_fsaverage=$path_script/registrations/hemi-${hemi_upper}_from-dhcpSym_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii
    if [ ! -f $registration_from_dhcpSym_to_fsaverage ]; then
        wb_command -surface-sphere-project-unproject \
            $path_script/registrations/dHCP_HCP-YA.MSMStrain.sphere.reg.surf.gii `# sphere-in: sphere with desired output mesh: dhcpsym registered to fslr` \
            $path_HCPtemplates_standardmeshatlases/${hemi_upper}.sphere.32k_fs_LR.surf.gii `# sphere-project-to: sphere that aligns with sphere-in: fs_lr` \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemi_upper}.sphere.32k_fs_LR.surf.gii `# sphere-unproject-from: sphere-project-to deformed to desired output space: fs_lr registered to fsaverage` \
            $registration_from_dhcpSym_to_fsaverage`# sphere-out: output sphere`
    fi

    # 2. Combine registration from native to dhcpSym to fsaverage
    registration_from_native_to_fsaverage=$path_output_data/sub-$sub/ses-$ses/space-dhcpSym/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii
    if [ ! -f $registration_from_native_to_fsaverage ]; then
        wb_command -surface-sphere-project-unproject \
            $path_output_data/sub-$sub/ses-$ses/space-dhcpSym/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-dhcpSym40_dens-32k_mode-sphere.reg40.surf.gii `# sphere-in: sphere with desired output mesh: native registered to dhcpSym` \
            $path_surfacetemplate/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii `# sphere-project-to: sphere that aligns with sphere-in: dhcpSym` \
            $path_script/registrations/hemi-${hemi_upper}_from-dhcpSym_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii`# sphere-unproject-from: sphere-project-to deformed to desired output space: dhcpSym registered to fsaverage` \
            $registration_from_native_to_fsaverage `# sphere-out: output sphere`
    fi

    # needed here as output: native sphere registered to fsaverage
    # 4. Resample retinotopy templates onto native
    mkdir -p $path_output_data/sub-$sub/ses-$ses/space-T2w/anat/
    wb_command -label-resample \
        $visareas_wang_template `# label-in: visareas_wang_template` \
        $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
        $path_output_data/sub-$sub/ses-$ses/space-dhcpSym/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii `# new sphere: sphere in register with current sphere and desired mesh` \
        ADAP_BARY_AREA \
        $path_output_data/sub-$sub/ses-$ses/space-T2w/anat/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_mesh-native_dens-native_desc-visualtopographywang2015_label-maxprob_dparc.label.gii \
        -area-surfs \
        $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
        $path_bids_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_space-T2w_midthickness.surf.gii

    wb_command -label-resample \
        $visareas_benson_template `# label-in: visareas_wang_template` \
        $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
        $path_output_data/sub-$sub/ses-$ses/space-dhcpSym/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii `# new sphere: sphere in register with current sphere and desired mesh` \
        ADAP_BARY_AREA \
        $path_output_data/sub-$sub/ses-$ses/space-T2w/anat/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_mesh-native_dens-native_desc-retinotbenson2014_label-visarea_dparc.label.gii \
        -area-surfs \
        $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
        $path_bids_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_space-T2w_midthickness.surf.gii

    wb_command -metric-resample \
        $ecc_template \
        $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
        $path_output_data/sub-$sub/ses-$ses/space-dhcpSym/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii `# new sphere: sphere in register with current sphere and desired mesh` \
        ADAP_BARY_AREA \
        $path_output_data/sub-$sub/ses-$ses/space-T2w/anat/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_mesh-native_dens-native_desc-eccentretinotbenson2014_seg.shape.gii \
        -area-surfs \
        $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
        $path_bids_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_space-T2w_midthickness.surf.gii

    wb_command -metric-resample \
        $angl_template \
        $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
        $path_output_data/sub-$sub/ses-$ses/space-dhcpSym/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii `# new sphere: sphere in register with current sphere and desired mesh` \
        ADAP_BARY_AREA \
        $path_output_data/sub-$sub/ses-$ses/space-T2w/anat/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_mesh-native_dens-native_desc-angleretinotbenson2014_seg.shape.gii \
        -area-surfs \
        $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
        $path_bids_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_space-T2w_midthickness.surf.gii

done
