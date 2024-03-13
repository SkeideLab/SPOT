# Population connective field modeling reveals retinotopic visual cortex organization in utero

- creates registration from fsaverage to individual surfaces[^10]
- resamples template of visual areas to individual surfaces
- projects functional data to individual surfaces
- applies cortical connective field model to functional data

## Prerequisites

- dHCP second release data[^1]
- transformations from individual anatomical volume space to group anatomical volume space[^2]
- 40 weeks anatomical template in volume space (must be the one that corresponds to the transformation)[^3]
- target surface template[^4]
- newMSM installation[^5]
- workbench command installation[^6]
- MIRTK installation[^7]
- HCP templates "standard mesh atlases"[^8]
- fsaverage subject from freesurfer[^9]
- FSL installation, callable from command line

## Preparing the individual surfaces

Adjust paths in `surfaceprep/transform_template_to_individual.sh` and run.

# 4. Project functional data to native surface

`hcp_surface.sh` Creates a ribbon mask, chooses voxels with good signal based on covariance, projects good voxels from ribbon to surface in native func space.

#TODO
- specify versions of other repos
- specify folder structure

# How to align an individual anatomical surface to a group surface






## Notes 

- uses dHCP data
- more info see https://docs.google.com/document/d/1pcI7m-MnbXyYh-rkaZqW_hl4WtpSIwAdHC8ZdyhtiGs/edit?usp=sharing

[^1]: [Developing Human Connectome Project](http://www.developingconnectome.org/data-release/second-data-release/)

[^2]: included in the dHCP data

[^3]: The group volumetric atlas used in the dHCP second release is called ['dhcp-volumetric-atlas-groupwise'](https://gin.g-node.org/BioMedIA/dhcp-volumetric-atlas-groupwise). Alternatively, the same atlas can be downloaded from Sean Fitzgibbon's ['augmented version'](https://git.fmrib.ox.ac.uk/seanf/dhcp-resources/-/blob/master/docs/dhcp-augmented-volumetric-atlas.md).

[^4]: The dHCP project has published the [dhcpSym atlas](https://brain-development.org/brain-atlases/atlases-from-the-dhcp-project/cortical-surface-template/), which can be used here.

[^5]: newMSM runs on FSL and can be installed from Renato Besenczi's github repo [newMSM](https://github.com/rbesenczi/newMSM).

[^6]: [Human Connectome Project Workbench Command](https://www.humanconnectome.org/software/workbench-command) 

[^7]: [Medical Image Registration ToolKit MIRTK](http://mirtk.github.io/), needs root access for installation.

[^8]: [Human Connectome Project templates](https://github.com/Washington-University/HCPpipelines/tree/master/global/templates/standard_mesh_atlases) can be downloaded from the HCPpipelines repo

[^9]: included in a freesurfer installation at `$FREESURFER_HOME/subjects/fsaverage`

[^10]: adapted from [dHCP template alignment repo](https://github.com/ecr05/dHCP_template_alignment)