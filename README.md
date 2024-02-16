Population connective field modeling reveals retinotopic visual cortex organization in utero

- uses dHCP data
- more info see https://docs.google.com/document/d/1pcI7m-MnbXyYh-rkaZqW_hl4WtpSIwAdHC8ZdyhtiGs/edit?usp=sharing

# 1. Getting the anatomical visual areas template (Wang 2015)

`transform_wang_to_gifti.py`: Creates a gifti label file that is usable in the connectome workbench universe and readable in freeview as annotation. Needs a copy of neuropythy source code. Adjust input and output paths. 


# 2. Creating a surface registration from individual native anatomical surface to dhcpSym surface space and mesh

Run using `apply_template.sh`, which calls `align_to_template_2nd_release.sh`

# 3. Warping the template to individual subject native surfaces using previously computed registration

`make_fsavg_to_native_warp.sh`

# 4. Project functional data to native surface

`hcp_surface.sh` Creates a ribbon mask, chooses voxels with good signal based on covariance, projects good voxels from ribbon to surface in native func space.

#TODO
- specify versions of other repos
- specify folder structure
- instruct how to download templates