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

## Prenatal data
Run 01_dataprep/alignment/separte_nii_to_gii.py first before run the program.

## HCP data
Run 00_HCP/HCP_file_prep.sh first before run the program.

# Code avaliabilty
All data analysis was performed with custom software using Python 3.11 and various third-party packages (matplotlib 3.8.3, numpy 1.26.3, pandas 2.2.0, scipy 1.12.0, and nibabel  5.2.0). Brain surfaces were visualized with nilearn (0.10.4).

# Data availability
The data for this study are available after institutional registration through public links in the National Institute of Mental Health Data Archive (NDA). Prenatal and neonatal data are available at https://nda.nih.gov/edit_collection.html?id=3955. Adolescent and adult data are available at https://nda.nih.gov/edit_collection.html?id=2846.



## Notes 

- uses dHCP data (01_dataprep, 02_ccfanalysis, 03_retinotopyanalzsis), uses HCP data (00_HCP), group comparison (03_retinotopyanalysis) and visualization (04_visualizazion)
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
