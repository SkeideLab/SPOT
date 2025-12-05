#!/bin/bash
set -u -x -e

# GOAL: registration from fsaverage to hcp fs_lr 
# Use this registration to resample the wang template to individual surfaces

#sub_id=$1
#path_bids_data=$2
#path_output_data=$3
#path_HCPtemplates_standardmeshatlases=$4
#path_fsaverage=$5

path_output_data="/data/p_02915/SPOT/00_HCP/template"
path_HCPtemplates_standardmeshatlases="/data/u_yoos_software/HCPpipelines/global/templates/standard_mesh_atlases/"
path_fsaverage="/data/p_02915/templates/template_fsaverage/fsaverage/"
# prepare output directory for subject in T2w space
#path_out_native=$path_output_data/${sub_id}/anat
#mkdir -p $path_out_native

for hemi in left right; do
    if [ $hemi = "left" ]; then
        hemi_upper="L"
    elif [ $hemi = "right" ]; then
        hemi_upper="R"
    fi

    # registrations necessary to make the transformation from fsaverage to native
    # will be created if they don't exist already    
    registration_from_fsLR_to_fsaverage="${path_HCPtemplates_standardmeshatlases}/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemi_upper}.sphere.32k_fs_LR.surf.gii"
    #registration_from_fsaverage_to_fsLR="/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/${sub_id}/surface_transforms/${sub_id}_hemi-${hemi_upper}_from-fsaverage_to-fsLR_dens-32k-sphere_reg.surf.gii"
    registration_from_fsaverage_to_fsLR="${path_output_data}/${hemi_upper}_from-fsaverage_to-fsLR_dens-32k-sphere_reg.surf.gii"
    #midthickness_from_fsLR="${path_bids_data}/${sub_id}.${hemi_upper}.midthickness_MSMAll.32k_fs_LR.surf.gii"
    #fs_LR_sphere="${path_bids_data}/${sub_id}.${hemi_upper}.sphere.32k_fs_LR.surf.gii"
    midthickness_from_fsLR="/data/p_02915/templates/fs_LR_32-master/fs_LR.32k.${hemi_upper}.midthickness.surf.gii"
    fs_LR_sphere="/data/p_02915/templates/fs_LR_32-master/fs_LR.32k.${hemi_upper}.sphere.surf.gii"

    # templates in fsaverage space
    # come with the repo, for more info see templates_retinotopy README
    path_templates_fsaverage="/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-${hemi_upper}_space-fsaverage_dens-164k"
    visareas_wang_template="${path_templates_fsaverage}_desc-visualtopographywang2015_label-maxprob_seg.label.gii"
    visareas_benson_template="${path_templates_fsaverage}_desc-retinotbenson2014_label-visarea_seg.label.gii"
    ecc_template="${path_templates_fsaverage}_desc-eccentretinotbenson2014_seg.shape.gii"
    angl_template="${path_templates_fsaverage}_desc-angleretinotbenson2014_seg.shape.gii"    

    # output
    path_out_native_hemi="${path_output_data}/hemi-${hemi_upper}_mesh-native_dens-native"
    visareas_wang_native="${path_out_native_hemi}_desc-visualtopographywang2015_label-maxprob_dparc.label.gii"
    visareas_benson_native="${path_out_native_hemi}_desc-retinotbenson2014_label-visarea_dparc.label.gii"
    ecc_native="${path_out_native_hemi}_desc-eccentretinotbenson2014_seg.shape.gii"
    pangle_native="${path_out_native_hemi}_desc-angleretinotbenson2014_seg.shape.gii"

    wb_command -surface-sphere-project-unproject \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# sphere-in: fsaverage` \
            $registration_from_fsLR_to_fsaverage `# sphere-project-to: sphere that aligns with sphere-in: fs_lr to fsaverage` \
            $fs_LR_sphere`# sphere-unproject-from: fsLR` \
            $registration_from_fsaverage_to_fsLR`# sphere-out: output sphere`

    
    # needed here as output: native sphere registered to fsaverage
    # 4. Resample retinotopy templates onto native
    if [ ! -f $visareas_wang_native ]; then
        wb_command -label-resample \
            $visareas_wang_template `# label-in: visareas_wang_template` \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
            $registration_from_fsLR_to_fsaverage \
            ADAP_BARY_AREA \
            $visareas_wang_native \
            -area-surfs \
            $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
            $midthickness_from_fsLR
    fi
    if [ ! -f $visareas_benson_native ]; then
        wb_command -label-resample \
            $visareas_benson_template `# label-in: visareas_wang_template` \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
            $registration_from_fsLR_to_fsaverage \
            ADAP_BARY_AREA \
            $visareas_benson_native \
            -area-surfs \
            $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
            $midthickness_from_fsLR
    fi
    if [ ! -f $ecc_native ]; then
        wb_command -metric-resample \
            $ecc_template \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
            $registration_from_fsLR_to_fsaverage \
            ADAP_BARY_AREA \
            $ecc_native \
            -area-surfs \
            $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
            $midthickness_from_fsLR
    fi
    if [ ! -f $pangle_native ]; then
        wb_command -metric-resample \
            $angl_template \
            $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi_upper}.164k_fsavg_${hemi_upper}.surf.gii `# current sphere: fsavg 164k ` \
            $registration_from_fsLR_to_fsaverage \
            ADAP_BARY_AREA \
            $pangle_native -area-surfs \
            $path_fsaverage/surf/${hemi:0:1}h_midthickness_surf.gii \
            $midthickness_from_fsLR
    fi

done
