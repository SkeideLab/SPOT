#!/bin/bash
# Resample ccf results from fsLR space to fsaverage space
set -u -x -e

sub=$1
path_anat_data=$2
path_output_dir=$3
path_HCPtemplates_standardmeshatlases=$4
path_fsaverage=$5
path_wbcommand=$6

threshold="00"
 
for hemi in L R; do

    for model in real simulated; do
        for param in r rss sigma v0i; do

            sphere_native="${path_anat_data}/${sub}.${hemi}.sphere.native.surf.gii"

            registration_from_native_to_fsaverage="/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/${sub}/surface_transforms/${sub}_hemi-${hemi}_from-native_to-fsaverage164k_dens-164k_mode-sphere_reg.surf.gii"
            registration_from_fsaverage_to_native="/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/${sub}/surface_transforms/${sub}_hemi-${hemi}_from-fsaverage_to-nativeT2w_dens-164k_mode-sphere_reg.surf.gii"
            native_midthickness="${path_anat_data}/${sub}.${hemi}.midthickness.native.surf.gii"

            res_file_native="${path_output_dir}/ccfmodel_var/${sub}/${sub}_hemi-${hemi}_mesh-native_dens-native_desc-${model}_${param}.gii"
            res_file_fsaverage="${path_output_dir}/ccfmodel_var/${sub}/${sub}_hemi-${hemi}_mesh-fsaverage_dens-164k_desc-${model}_${param}.gii"
            
            $path_wbcommand -surface-sphere-project-unproject \
                $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi}.164k_fsavg_${hemi}.surf.gii \
                $registration_from_native_to_fsaverage `# sphere-project-to: sphere that aligns with sphere-in` \
                $sphere_native `# sphere-unproject-from: sphere-project-to deformed to desired output space: native (not rotated)` \
                $registration_from_fsaverage_to_native `# sphere-out: output sphere`

            $path_wbcommand -metric-resample \
                $res_file_native `# label-in: ` \
                $sphere_native`# current sphere:  ` \
                $registration_from_fsaverage_to_native`# new sphere: sphere in register with current sphere and desired mesh` \
                ADAP_BARY_AREA \
                $res_file_fsaverage \
                -area-surfs \
                $native_midthickness \
                $path_fsaverage/surf/$(echo $hemi | tr '[:upper:]' '[:lower:]')h_midthickness.surf.gii
        done
    done
done



