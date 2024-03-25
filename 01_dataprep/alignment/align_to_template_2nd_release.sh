#!/bin/bash

# script to align native surfaces with template space
set -x -u -e
Usage() {
    echo "align_to_template.sh <topdir> <subjid> <session> <age> <volumetric template> <volumetric template name> <surface template> <surface template name> <pre_rotation> <outdir>  <config> <script dir> < MSM bin> <wb bin>"
    echo " script to align native surfaces with template space & resample native surfaces with template topology "
    echo " input args: "
    echo " topdir: top directory where subject directories are located "
    echo " subjid : subject id "
    echo " session: subject scan session "
    echo " age: in weeks gestation - this will determine which week of the spatio-temporal template the data will first mapped to"
    echo " template volume: template T2 40 week volume "
    echo " volumetric template name: extdhcp40wk or dhcp40wk"
    echo " surface template: path to the top level directory of the dHCP surface template"
    echo " surface template name: dhcpSym or dhcpASym"
    echo " pre_rotation : txt file containing rotational transform between MNI and FS_LR space (i.e. file rotational_transforms/week40_toFS_LR_rot.%hemi%.txt  ) "
    echo " outdir : base directory where output will be sent "
    echo " config : base config file "
    echo " MSM bin: msm binary"
    echo " wb bin : workbench binary"
    echo "mirtk bin : mirtk binary "
    echo "output: surface registrations"
}

if [ "$#" -lt 11 ]; then
    echo "$#"
    Usage
    exit
fi
path_script=$(dirname $0)

topdir=$1
shift
subjid=$1
shift
session=$1
shift
age=$1
shift
templatevolume=$1
shift
templatevolumename=$1
shift
templatespherepath=$1
shift
templatespherename=$1
shift
pre_rotation=$1
shift
outdir=$1
shift
config=$1
shift
MSMBIN=$1
shift
WB_BIN=$1
shift
mirtk_BIN=$1
shift

########## DEFINE PATHS TO VARIABLES ##########
# TODO define mesh and space suffixes for all surfaces
#inputs
nativedir=${topdir}/sub-${subjid}/ses-$session/anat
native_volume=${nativedir}/sub-${subjid}_ses-${session}_desc-restore_T2w.nii.gz
native_sphere=${nativedir}/sub-${subjid}_ses-${session}_hemi-%hemi%_space-T2w_sphere.surf.gii
native_data=${nativedir}/sub-${subjid}_ses-${session}_hemi-%hemi%_space-T2w_sulc.shape.gii

# inputs (template)
template00wk_sphere=$templatespherepath/week-${age}_hemi-%hemi%_space-${templatespherename}_dens-32k_sphere.surf.gii
template00wk_data=$templatespherepath/week-${age}_hemi-%hemi%_space-${templatespherename}_dens-32k_sulc.shape.gii

#outputs
sub_templatespace_dir=${outdir}/sub-${subjid}/ses-$session/space-${templatespherename}
mkdir -p $sub_templatespace_dir/volume_dofs $sub_templatespace_dir/surface_transforms $sub_templatespace_dir/anat
native_rot_sphere=${nativedir}/sub-${subjid}_ses-${session}_hemi-%hemi%_space-fslr_sphere.rot.surf.gii
outname=$sub_templatespace_dir/surface_transforms/sub-${subjid}_ses-${session}_hemi-%hemi%_from-native_to-${templatespherename}40_dens-32k_mode-
transformed_sphere=${outname}sphere.reg40.surf.gii

for hemi in left right; do

    # capitalize and extract first letter of hemi (left --> L)
    hemi_upper_tmp=${hemi:0:1}
    hemi_upper=${hemi_upper_tmp^}

    # swap in correct hemisphere label
    pre_rotation_hemi=$(echo ${pre_rotation} | sed "s/%hemi%/$hemi_upper/g")
    native_sphere_hemi=$(echo ${native_sphere} | sed "s/%hemi%/$hemi_upper/g")
    native_data_hemi=$(echo ${native_data} | sed "s/%hemi%/$hemi_upper/g")
    template00wk_sphere_hemi=$(echo ${template00wk_sphere} | sed "s/%hemi%/$hemi/g")
    template00wk_data_hemi=$(echo ${template00wk_data} | sed "s/%hemi%/$hemi/g")
    native_rot_sphere_hemi=$(echo ${native_rot_sphere} | sed "s/%hemi%/$hemi_upper/g")
    outname_hemi=$(echo ${outname} | sed "s/%hemi%/$hemi_upper/g")
    transformed_sphere_hemi=$(echo ${transformed_sphere} | sed "s/%hemi%/$hemi_upper/g")

    if [ ! -f ${transformed_sphere_hemi} ]; then
        ########## ROTATE LEFT AND RIGHT HEMISPHERES INTO APPROXIMATE ALIGNMENT WITH MNI SPACE ##########
        $path_script/pre_rotation.sh \
            $native_volume \
            $native_sphere_hemi \
            $templatevolume \
            $pre_rotation_hemi \
            $sub_templatespace_dir/volume_dofs/${subjid}-${session}_space-$templatevolumename.dof \
            $native_rot_sphere_hemi \
            $mirtk_BIN $WB_BIN

        ########## RUN MSM NON-LINEAR ALIGNMENT TO TEMPLATE FOR LEFT AND RIGHT HEMISPHERES ##########
        indata=$native_data_hemi
        inmesh=$native_rot_sphere_hemi
        refmesh=$template00wk_sphere_hemi
        refdata=$template00wk_data_hemi

        ${MSMBIN} \
            --inmesh=${inmesh} \
            --refmesh=${refmesh} \
            --indata=${indata} \
            --refdata=${refdata} \
            -o ${outname_hemi} \
            --conf=${config} \
            --verbose

        if [ "$age" == "40" ]; then
            # rename to emphasize registration to 40 (sphere.reg40.surf.gii)
            mv ${outname_hemi}sphere.reg.surf.gii ${transformed_sphere_hemi}
        else
            # need to concatenate msm warp to local template with warp from local template to 40 week template
            ${WB_BIN} -surface-sphere-project-unproject \
                ${outname_hemi}sphere.reg.surf.gii \
                $refmesh \
                $templatespherepath/week-to-40-registrations/${hemi}.${age}-to-40/${hemi}.${age}-to-40sphere.reg.surf.gii \
                $transformed_sphere_hemi ### LZJW added hemi and changed filepath for between template ###
        fi

        # the output sphere represents the full warp from Native to 40 week template space - save this
        cp "$transformed_sphere_hemi" "$nativedir"
    else
        echo "Transformed sphere $transformed_sphere_hemi already exists, skipping..."

    fi
done
