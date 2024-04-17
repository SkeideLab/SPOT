#!/bin/bash
set -u -x -e

sub=$1
ses=$2
path_anat_data=$3
path_output_dir=$4
path_HCPtemplates_standardmeshatlases=$5
path_fsaverage=$6
path_wbcommand=$7

threshold="01"

for hemi in L R; do
    if [ $hemi = "L" ]; then
        hemi_down="left"
    elif [ $hemi = "R" ]; then
        hemi_down="right"
    fi

    for model in real simulated; do
        for param in eccentricity polarangle; do

            #TODO should we use rotated sphere here?
            # should not matter since the result is just a label file without geometric info? not sure

            sphere_native="${path_anat_data}/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi_down}_sphere.surf.gii"
            registration_from_native_to_fsaverage="${path_output_dir}/dhcp_surface/sub-$sub/ses-$ses/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii"
            registration_from_fsaverage_to_native="${path_output_dir}/dhcp_surface/sub-$sub/ses-$ses/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi}_from-fsaverage_to-nativeT2w_dens-164k_mode-sphere_reg.surf.gii"

            res_file_native="${path_output_dir}/ccfmodel/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_hemi-${hemi}_mesh-native_dens-native_label-${param}_desc-${model}_roi-v2th${threshold}_metric.gii"
            res_file_fsaverage="${path_output_dir}/ccfmodel/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_hemi-${hemi}_mesh-fsaverage_dens-164k_label-${param}_desc-${model}_roi-v2th${threshold}_metric.gii"

            $path_wbcommand -surface-sphere-project-unproject \
                $path_HCPtemplates_standardmeshatlases/resample_fsaverage/fsaverage_std_sphere.${hemi}.164k_fsavg_${hemi}.surf.gii \
                $registration_from_native_to_fsaverage `# sphere-project-to: sphere that aligns with sphere-in` \
                $sphere_native `# sphere-unproject-from: sphere-project-to deformed to desired output space: native (not rotated)` \
                $registration_from_fsaverage_to_native `# sphere-out: output sphere`

            $path_wbcommand -label-resample \
                $res_file_native `# label-in: ` \
                $sphere_native`# current sphere: replace with rotated sphere??? ` \
                $registration_from_fsaverage_to_native`# new sphere: sphere in register with current sphere and desired mesh` \
                ADAP_BARY_AREA \
                $res_file_fsaverage \
                -area-surfs \
                $path_anat_data/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi_down}_midthickness.surf.gii \
                $path_fsaverage/surf/$(echo $hemi | tr '[:upper:]' '[:lower:]')h_midthickness_surf.gii
        done
    done
done
