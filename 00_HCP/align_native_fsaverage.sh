#!/bin/bash
# script for third data release
# script to align native surfaces with template space
set -x -u -e
Usage() {
    echo "align_to_template.sh <topdir> <subjid> <session> <age> <volumetric template> <volumetric template name> <surface template> <surface template name> <pre_rotation> <outdir>  <config> <script dir> < MSM bin> <wb bin>"
    echo " script to align native surfaces with template space & resample native surfaces with template topology "
    echo " input args: "
    echo " topdir: top directory where subject directories are located "
    echo " subjid : subject id "
    echo " outdir : base directory where output will be sent "
    echo " config : base config file "
    echo " MSM bin: msm binary"
    echo " wb bin : workbench binary"
    echo "mirtk bin : mirtk binary "
    echo "output: surface registrations"
}
 
if [ "$#" -lt 7 ]; then
    echo "$#"
    Usage
    exit
fi
path_script=$(dirname $0)

topdir=$1
shift
subjid=$1
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
#inputs
native_sphere=${topdir}/${subjid}.%hemi%.sphere.native.surf.gii
native_data=${topdir}/${subjid}.%hemi%.sulc.native.shape.gii

# inputs (template)fetal.week21.left.sulc.shape.gii
fsLR_sphere=/data/p_02915/templates/fs_LR_32-master/fs_LR.32k.%hemi%.sphere.surf.gii
fsLR_data=/data/p_02915/templates/fs_LR_32-master/fs_LR.32k.%hemi%.sulc.shape.gii
#outputs
sub_output_dir=${outdir}/${subjid}
mkdir -p $sub_output_dir/volume_dofs $sub_output_dir/surface_transforms
native_rot_sphere=${topdir}/${subjid}.%hemi%.sphere.rot.native.surf.gii
outname=$sub_output_dir/surface_transforms/${subjid}_hemi-%hemi%_from-native_to-fsLR_dens-32k_mode-
transformed_sphere=${outname}sphere.reg.surf.gii

for hemi in left right; do

    # capitalize and extract first letter of hemi (left --> L)
    hemi_upper_tmp=${hemi:0:1}
    hemi_upper=${hemi_upper_tmp^}

    # swap in correct hemisphere label
    native_sphere_hemi=$(echo ${native_sphere} | sed "s/%hemi%/$hemi_upper/g")
    native_data_hemi=$(echo ${native_data} | sed "s/%hemi%/$hemi_upper/g")
    fsLR_sphere_hemi=$(echo ${fsLR_sphere} | sed "s/%hemi%/$hemi_upper/g")
    fsLR_data_hemi=$(echo ${fsLR_data} | sed "s/%hemi%/$hemi_upper/g")
    native_rot_sphere_hemi=$(echo ${native_rot_sphere} | sed "s/%hemi%/$hemi_upper/g")
    outname_hemi=$(echo ${outname} | sed "s/%hemi%/$hemi_upper/g")
    transformed_sphere_hemi=$(echo ${transformed_sphere} | sed "s/%hemi%/$hemi_upper/g")
    if [ ! -f ${transformed_sphere_hemi} ]; then
        ########## RUN MSM NON-LINEAR ALIGNMENT TO TEMPLATE FOR LEFT AND RIGHT HEMISPHERES ##########
        indata=$native_data_hemi
        inmesh=$native_rot_sphere_hemi
        refmesh=$fsLR_sphere_hemi
        refdata=$fsLR_data_hemi

        ${MSMBIN} \
            --inmesh=${inmesh} \
            --refmesh=${refmesh} \
            --indata=${indata} \
            --refdata=${refdata} \
            -o ${outname_hemi} \
            --verbose
        

    else
        echo "Transformed sphere $transformed_sphere_hemi already exists, skipping..."

    fi
    
done
