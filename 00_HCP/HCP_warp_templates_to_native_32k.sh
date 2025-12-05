#!/bin/bash
set -u -x -e

# GOAL: registration from fsaverage to hcp fs_lr to dhcpSym40 to native
# Use this registration to resample the wang template to individual surfaces

path_script=$(dirname $0)
sub=$1
path_bids_data=$2
path_output_data=$3
path_HCPtemplates_standardmeshatlases=$4
path_fsaverage=$5

# prepare output directory for subject in T2w space
path_out_native=$path_output_data/$sub/anat
mkdir -p $path_out_native $path_output_data/$sub/surface_transforms

for hemi in left right; do
    if [ $hemi = "left" ]; then
        hemi_upper="L"
    elif [ $hemi = "right" ]; then
        hemi_upper="R"
    fi

    # registrations necessary to make the transformation from fsaverage to native
    # will be created if they don't exist already    
    registration_from_native_to_fsaverage="$path_output_data/$sub/surface_transforms/${sub}_hemi-${hemi_upper}_from-native_to-fsaverage164k_dens-32k_mode-sphere_reg.surf.gii"

    # templates in fsaverage space
    # come with the repo, for more info see templates_retinotopy README
    path_templates_fsaverage="$path_script/template/hemi-${hemi_upper}_space-fsaverage_dens-164k"
    visareas_wang_template="${path_templates_fsaverage}_desc-visualtopographywang2015_label-maxprob_seg.label.gii"
    visareas_benson_template="${path_templates_fsaverage}_desc-retinotbenson2014_label-visarea_seg.label.gii"
    ecc_template="${path_templates_fsaverage}_desc-eccentretinotbenson2018_seg.shape.gii"
    angl_template="${path_templates_fsaverage}_desc-angleretinotbenson2018_seg.shape.gii"

    # output
    path_out_native_hemi="$path_out_native/${sub}_hemi-${hemi_upper}_mesh-native_dens-32k-native"
    visareas_wang_native="${path_out_native_hemi}_desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
    visareas_benson_native="${path_out_native_hemi}_desc-retinotbenson2014_label-visarea_dparc.label.gii"
    ecc_native="${path_out_native_hemi}_desc-eccentretinotbenson2018_seg.shape.gii"
    pangle_native="${path_out_native_hemi}_desc-angleretinotbenson2018_seg.shape.gii"

    # 1. Combine registration from native to fs_lr to fsaverage
    
    wb_command -surface-sphere-project-unproject \
        "$path_output_data/$sub/surface_transforms/${sub}_hemi-${hemi_upper}_from-native_to-fsLR_dens-32k_mode-sphere.reg.surf.gii" `# sphere-in: sphere with desired output mesh: dhcpsym registered to fslr` \
        $path_HCPtemplates_standardmeshatlases/${hemi_upper}.sphere.32k_fs_LR.surf.gii `# sphere-project-to: sphere that aligns with sphere-in: fs_lr` \
        $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemi_upper}.sphere.32k_fs_LR.surf.gii `# sphere-unproject-from: sphere-project-to deformed to desired output space: fs_lr registered to fsaverage` \
        $registration_from_native_to_fsaverage`# sphere-out: output sphere`
    
    # needed here as output: native sphere registered to fsaverage
    # 2. Resample retinotopy templates onto native
    
    wb_command -label-resample \
        $visareas_wang_template `# label-in: visareas_wang_template` \
        $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
        $registration_from_native_to_fsaverage \
        ADAP_BARY_AREA \
        $visareas_wang_native \
        -area-surfs \
        $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
        "/data/pt_02880/HCP_D/fmriresults01/${sub}/MNINonLinear/Native/${sub}.${hemi_upper}.midthickness.native.surf.gii"
        
    wb_command -label-resample \
        $visareas_benson_template `# label-in: visareas_wang_template` \
        $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
        $registration_from_native_to_fsaverage \
        ADAP_BARY_AREA \
        $visareas_benson_native \
        -area-surfs \
        $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
        "/data/pt_02880/HCP_D/fmriresults01/${sub}/MNINonLinear/Native/${sub}.${hemi_upper}.midthickness.native.surf.gii"
    wb_command -metric-resample \
        $ecc_template \
        $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
        $registration_from_native_to_fsaverage \
        ADAP_BARY_AREA \
        $ecc_native \
        -area-surfs \
        $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
        "/data/pt_02880/HCP_D/fmriresults01/${sub}/MNINonLinear/Native/${sub}.${hemi_upper}.midthickness.native.surf.gii"
    wb_command -metric-resample \
        $angl_template \
        $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
        $registration_from_native_to_fsaverage \
        ADAP_BARY_AREA \
        $pangle_native -area-surfs \
        $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
        "/data/pt_02880/HCP_D/fmriresults01/${sub}/MNINonLinear/Native/${sub}.${hemi_upper}.midthickness.native.surf.gii"
    
done
