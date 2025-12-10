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
    registration_from_native_to_dhcpSym="$path_output_data/sub-$sub/ses-$ses/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-fetal_to-dhcpSym40_dens-32k_mode-sphere.reg40.surf.gii"
    if [ ! -f $registration_from_native_to_dhcpSym ]; then
        echo "Registration from native to dhcpSym is missing. Please execute alignment folder main script first."
        echo "Terminating..."
        exit 1
    fi

    # templates in fsaverage space
    # come with the repo, for more info see templates_retinotopy README
    path_templates_fsaverage="$path_script/templates_retinotopy/hemi-${hemi_upper}_space-fsaverage_dens-164k"
    visareas_wang_template="${path_templates_fsaverage}_desc-visualtopographywang2015_label-maxprob_seg.label.gii"
    visareas_benson_template="${path_templates_fsaverage}_desc-retinotbenson2014_label-visarea_seg.label.gii"
    ecc_template="${path_templates_fsaverage}_desc-eccentretinotbenson2018_seg.shape.gii"
    angl_template="${path_templates_fsaverage}_desc-angleretinotbenson2018_seg.shape.gii"

    # output
    path_out_native_hemi="$path_out_native/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_mesh-native_dens-native"
    visareas_wang_native="${path_out_native_hemi}_desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
    visareas_benson_native="${path_out_native_hemi}_desc-retinotbenson2014_label-visarea_dparc.label.gii"
    ecc_native="${path_out_native_hemi}_desc-eccentretinotbenson2018_seg.shape.gii"
    pangle_native="${path_out_native_hemi}_desc-angleretinotbenson2018_seg.shape.gii"

    # 1. Combine registration from dhcpSym to fs_lr to fsaverage
    if [ ! -f $registration_from_dhcpSym_to_fsaverage ]; then
        wb_command -surface-sphere-project-unproject \
            $path_script/standard_registrations/dHCP_HCP-YA.MSMStrain.${hemi_upper}.sphere.reg.surf.gii `# sphere-in: sphere with desired output mesh: dhcpsym registered to fslr` \
            $path_HCPtemplates_standardmeshatlases/${hemi_upper}.sphere.32k_fs_LR.surf.gii `# sphere-project-to: sphere that aligns with sphere-in: fs_lr` \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemi_upper}.sphere.32k_fs_LR.surf.gii `# sphere-unproject-from: sphere-project-to deformed to desired output space: fs_lr registered to fsaverage` \
            $registration_from_dhcpSym_to_fsaverage`# sphere-out: output sphere`
    fi

    # 2. Combine registration from native to dhcpSym to fsaverage
    if [ ! -f $registration_from_native_to_fsaverage ]; then
        wb_command -surface-sphere-project-unproject \
            $registration_from_native_to_dhcpSym `# sphere-in: sphere with desired output mesh: native registered to dhcpSym` \
            /data/p_02915/templates/template_corticalsurfaceneonatessym_williams2023_dhcp/dhcpSym_template/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii `# sphere-project-to: sphere that aligns with sphere-in: dhcpSym` \
            $registration_from_dhcpSym_to_fsaverage`# sphere-unproject-from: sphere-project-to deformed to desired output space: dhcpSym registered to fsaverage` \
            $registration_from_native_to_fsaverage `# sphere-out: output sphere`
    fi

    # needed here as output: native sphere registered to fsaverage
    # 4. Resample retinotopy templates onto native
    if [ ! -f $visareas_wang_native ]; then
        wb_command -label-resample \
            $visareas_wang_template `# label-in: visareas_wang_template` \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
            $registration_from_native_to_fsaverage \
            ADAP_BARY_AREA \
            $visareas_wang_native \
            -area-surfs \
            $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
            $path_bids_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_midthickness.surf.gii
    fi
    if [ ! -f $visareas_benson_native ]; then
        wb_command -label-resample \
            $visareas_benson_template `# label-in: visareas_wang_template` \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
            $registration_from_native_to_fsaverage \
            ADAP_BARY_AREA \
            $visareas_benson_native \
            -area-surfs \
            $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
            $path_bids_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_midthickness.surf.gii
    fi
    if [ ! -f $ecc_native ]; then
        wb_command -metric-resample \
            $ecc_template \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
            $registration_from_native_to_fsaverage \
            ADAP_BARY_AREA \
            $ecc_native \
            -area-surfs \
            $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
            $path_bids_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_midthickness.surf.gii
    fi
    if [ ! -f $pangle_native ]; then
        wb_command -metric-resample \
            $angl_template \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
            $registration_from_native_to_fsaverage \
            ADAP_BARY_AREA \
            $pangle_native -area-surfs \
            $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
            $path_bids_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_midthickness.surf.gii
    fi

done
