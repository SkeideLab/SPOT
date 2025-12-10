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
path_wbcommand=$8
 
# prepare output directory for subject in T2w space
path_out_native=$path_output_data/sub-$sub/ses-$ses/func
mkdir -p $path_out_native
logfile=$path_output_data/sub-$sub/ses-$ses/logs/surfaceprep_$(date +%Y-%m-%dT%H%M).log
mkdir -p $(dirname $logfile)
echo "Redirecting outputs to logfile $logfile..."

for hemi in left right; do
    if [ $hemi = "left" ]; then
        hemi_upper="L"
    elif [ $hemi = "right" ]; then
        hemi_upper="R"
    fi

    # registrations necessary to make the transformation from fsaverage to native
    # will be created if they don't exist already
    registration_from_dhcpSym_to_fsLR=$path_script/standard_registrations/dHCP_HCP-YA.MSMStrain.${hemi_upper}.sphere.reg.surf.gii
    registration_from_native_to_fsLR="$path_output_data/sub-$sub/ses-$ses/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-fsLR_dens-32k_mode-sphere_reg.surf.gii"

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
    func_native="${path_out_native}/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_mesh-native_bold.func.gii"

    # output
    func_fsLR="${path_out_native}/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_mesh-fsLR_bold.func.gii"

    # 2. Combine registration from native to dhcpSym to fsaverage
    $path_wbcommand -surface-sphere-project-unproject \
        $registration_from_native_to_dhcpSym `# sphere-in: sphere with desired output mesh: native registered to dhcpSym` \
        /data/p_02915/templates/template_corticalsurfaceneonatessym_williams2023_dhcp/dhcpSym_template/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii `# sphere-project-to: sphere that aligns with sphere-in: dhcpSym` \
        $registration_from_dhcpSym_to_fsLR`# sphere-unproject-from: sphere-project-to deformed to desired output space: dhcpSym registered to fsaverage` \
        $registration_from_native_to_fsLR `# sphere-out: output sphere`
    

    # needed here as output: native sphere registered to fsaverage
    # 4. Resample retinotopy templates onto native
    
    $path_wbcommand -metric-resample \
        $func_native `# label-in: visareas_wang_template` \
        /data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_anat_pipeline/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_hemi-${hemi}_sphere.surf.gii \
        $registration_from_native_to_fsLR \
        BARYCENTRIC \
        $func_fsLR 
    

done
>>$logfile 2>&1

