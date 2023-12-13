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
    echo " script dir: path to scripts"
    echo " MSM bin: msm binary"
    echo " wb bin : workbench binary"
    echo "mirtk bin : mirtk binary "
    echo "output: 1) surface registrations; 2)  native giftis resampled with template topology "
}

/data/u_kieslinger_software/code/dHCP_template_alignment/surface_to_template_alignment/align_to_template_3rd_release.sh \
    /data/p_02495/dhcp_derivatives/dhcp_anat_pipeline \
    CC00058XX09 \
    11300 \
    41 \
    /data/p_02495/templates/template_augmentedvolumetricatlas_dhcp/atlas/T2/template-40.nii.gz \
    dhcp40wk \
    /data/p_02495/templates/template_corticalsurfaceneonatessym_williams2023_dhcp/dhcpSym_template \
    dhcpSym \
    /data/u_kieslinger_software/code/dHCP_template_alignment/rotational_transforms/week40_toFS_LR_rot.%hemi%.txt \
    /data/p_02495/dhcp_derivatives/dhcp_surface \
    /data/u_kieslinger_software/code/dHCP_template_alignment/configs/config_subject_to_40_week_template \
    /data/u_kieslinger_software/code/dHCP_template_alignment \
    /data/u_kieslinger_software/fsldevdir/bin/newmsm \
    /bin/wb_command \
    /afs/cbs.mpg.de/software/mirtk/0.20231123/debian-bullseye-amd64/bin/mirtk

#/data/u_kieslinger_software/miniconda3/envs/fsl_env/bin/msm \
