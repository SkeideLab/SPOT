#!/bin/bash
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

    for model in real; do
        for param in eccentricity polarangle; do
            sphere_native="/data/p_02915/templates/fs_LR_32-master/fs_LR.32k.${hemi}.sphere.surf.gii"
            registration_from_native_to_fsaverage="/data/p_02915/SPOT/00_HCP/template/${hemi}_from-fsaverage_to-fsLR_dens-32k-sphere_reg.surf.gii"
            native_midthickness="/data/p_02915/templates/fs_LR_32-master/fs_LR.32k.${hemi}.midthickness.surf.gii"

            res_file_native="${path_output_dir}/ccfmodel_fsLR/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_hemi-${hemi}_mesh-native_dens-native_label-${param}_desc-${model}_roi-v2th${threshold}_metric.gii"
            res_file_fsaverage="${path_output_dir}/ccfmodel_fsLR/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_hemi-${hemi}_mesh-fsaverage_dens-164k_label-${param}_desc-${model}_roi-v2th${threshold}_metric.gii"
            
            $path_wbcommand -label-resample \
                $res_file_native `# label-in: ` \
                $sphere_native`# current sphere:  ` \
                $registration_from_native_to_fsaverage`# new sphere: sphere in register with current sphere and desired mesh` \
                ADAP_BARY_AREA \
                $res_file_fsaverage \
                -area-surfs \
                $native_midthickness \
                $path_fsaverage/surf/$(echo $hemi | tr '[:upper:]' '[:lower:]')h_midthickness_surf.gii
        done
    done
done
