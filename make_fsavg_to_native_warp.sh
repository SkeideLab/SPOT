# “fsaverage?_std_sphere” are the standard spheres from freesurfer, converted to GIFTI format.
# “fs_LR-deformed_to-fsaverage” are the fs_LR spheres registered to fsaverage

# GOAL: registration from fsaverage to hcp fs_lr to dhcpSym40 to native

# wb_command -surface-sphere-project-unproject \
#     a `# sphere-in: sphere with desired output mesh` \
#     b `# sphere-project-to: sphere that aligns with sphere-in` \
#     c `# sphere-unproject-from: sphere-project-to deformed to desired output space` \
#     d `# sphere-out: output sphere`
sub="CC00060XX03"
ses="12501"
hemi_upper="L"
hemi="left"

# 1. Combine registration from dhcpSym to fs_lr to fsaverage
wb_command -surface-sphere-project-unproject \
    /data/u_kieslinger_software/code/CorticalAsymmetry/dHCP_HCP-YA/dHCP_HCP-YA.MSMStrain.sphere.reg.surf.gii `# sphere-in: sphere with desired output mesh: dhcpsym registered to fslr` \
    /data/u_kieslinger_software/code/HCPpipelines/global/templates/standard_mesh_atlases/${hemi_upper}.sphere.32k_fs_LR.surf.gii `# sphere-project-to: sphere that aligns with sphere-in: fs_lr` \
    /data/u_kieslinger_software/code/HCPpipelines/global/templates/standard_mesh_atlases/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemi_upper}.sphere.32k_fs_LR.surf.gii `# sphere-unproject-from: sphere-project-to deformed to desired output space: fs_lr registered to fsaverage` \
    /data/u_kieslinger_software/code/flux/registrations/hemi-${hemi_upper}_from-dhcpSym_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii`# sphere-out: output sphere`

# 2. Combine registration from native to dhcpSym to fsaverage
wb_command -surface-sphere-project-unproject \
    /data/p_02495/dhcp_derivatives/dhcp_surface/sub-$sub/ses-$ses/space-dhcpSym_32k/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-dhcpSym40_dens-32k_mode-sphere.LR.reg.surf.gii `# sphere-in: sphere with desired output mesh: native registered to dhcpSym` \
    /data/p_02495/templates/template_corticalsurfaceneonatessym_williams2023_dhcp/dhcpSym_template/week-40_hemi-${hemi}_space-dhcpSym_dens-32k_sphere.surf.gii `# sphere-project-to: sphere that aligns with sphere-in: dhcpSym` \
    /data/u_kieslinger_software/code/flux/registrations/hemi-${hemi_upper}_from-dhcpSym_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii`# sphere-unproject-from: sphere-project-to deformed to desired output space: dhcpSym registered to fsaverage` \
    /data/p_02495/dhcp_derivatives/dhcp_surface/sub-$sub/ses-$ses/space-dhcpSym_32k/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii `# sphere-out: output sphere`

# 3. Inverse registration: fsaverage -> native
wb_command -surface-sphere-project-unproject \
    a `# sphere-in: sphere with desired output mesh: fsaverage 32k ` \
    /data/p_02495/dhcp_derivatives/dhcp_surface/sub-$sub/ses-$ses/space-dhcpSym_32k/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-native_to-fsaverage32k_dens-32k_mode-sphere_reg.surf.gii `# sphere-project-to: sphere that aligns with sphere-in` \
    /data/p_02495/dhcp_derivatives/dhcp_anat_pipeline/sub-$sub/ses-$ses/anat/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_space-T2w_sphere.surf.gii `# sphere-unproject-from: sphere-project-to deformed to desired output space: native (not rotated)` \
    /data/p_02495/dhcp_derivatives/dhcp_surface/sub-$sub/ses-$ses/space-T2w/surface_transforms/sub-${sub}_ses-${ses}_hemi-${hemi_upper}_from-fsaverage32k_to-nativeT2w_dens-32k_mode-sphere_reg.surf.gii `# sphere-out: output sphere`

# /data/u_kieslinger_software/code/neuropythy/neuropythy/lib/data/fsaverage/surf/lh.wang15_mplbl.v1_0.mgz
