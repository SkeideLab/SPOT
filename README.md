# Population connective field modeling reveals retinotopic visual cortex organization in utero

- creates registration from fsaverage to individual surfaces[^13]
- resamples template of visual areas to individual surfaces
- projects functional data to individual surfaces
- applies cortical connective field model to functional data

## Prerequisites

Data and templates
- dHCP fourth release data[^1]
- The Lifespan Human Connectome Project Development (HCP-D) 2.0 data [^2]
- transformations from individual anatomical volume space to group anatomical volume space[^3]
- 40 weeks anatomical template in volume space (must be the one that corresponds to the transformation)[^4]
- target surface template[^5]
- fsaverage subject from freesurfer[^6] 
- Developing Human Connectome Project spatio-temporal surface atlas of the fetal brain [^7]
- HCP templates "standard mesh atlases"[^8]

Software
- newMSM (v0.6.3-BETA) installation[^9]
- workbench command (v2.0.0) installation[^10]
- MIRTK (v0.20231123) installation[^11]
- FSL (v6.0.3) installation, callable from command line[^12]

# Code avaliabilty
All data analysis was performed with custom software using Python 3.11 and various third-party packages (matplotlib 3.8.3, numpy 1.26.3, pandas 2.2.0, scipy 1.12.0, and nibabel  5.2.0). Brain surfaces were visualized with nilearn (0.10.4).

# Data availability
The data for this study are available after institutional registration through public links in the National Institute of Mental Health Data Archive (NDA). Prenatal and neonatal data are available at https://nda.nih.gov/edit_collection.html?id=3955. Adolescent and adult data are available at https://nda.nih.gov/edit_collection.html?id=2846.



## Notes 
- Please change the path before you run the code

### dHCP dataset (prenatal and neonates)

for prenatal dataset, run 01_dataprep/alignment/separte_nii_to_gii.py first before run the program.
run `run_all.py` which contains codes in (01_dataprep, 02_ccfanalysis, 03_retinotopyanalysis)

### HCP dataset

run 00_HCP/HCP_file_prep.sh first before run the program.
run `00_HCP/HCP_run_all.py`

### Group comparison 

see scripts in `03_retinotopyanalysis`

### Visualization

see scripts in `04_visualization`

[^1]: [Developing Human Connectome Project](https://biomedia.github.io/dHCP-release-notes/)

[^2]: [Human Connectome Project Development](https://www.humanconnectome.org/study/hcp-lifespan-development)

[^3]: included in the dHCP data

[^4]: The group volumetric atlas used in the dHCP second release is called ['dhcp-volumetric-atlas-groupwise'](https://gin.g-node.org/BioMedIA/dhcp-volumetric-atlas-groupwise). Alternatively, the same atlas can be downloaded from Sean Fitzgibbon's ['augmented version'](https://git.fmrib.ox.ac.uk/seanf/dhcp-resources/-/blob/master/docs/dhcp-augmented-volumetric-atlas.md).

[^5]: The dHCP project has published the [dhcpSym atlas](https://brain-development.org/brain-atlases/atlases-from-the-dhcp-project/cortical-surface-template/), which can be used here.

[^6]: included in a freesurfer installation at `$FREESURFER_HOME/subjects/fsaverage`

[^7]: Developing Human Connectome Project spatio-temporal surface atlas of the fetal brain (https://doi.gin.g-node.org/10.12751/g-node.qj5hs7/)

[^8]: [Human Connectome Project templates](https://github.com/Washington-University/HCPpipelines/tree/master/global/templates/standard_mesh_atlases) can be downloaded from the HCPpipelines repo

[^9]: newMSM runs on FSL and can be installed from Renato Besenczi's github repo [newMSM](https://github.com/rbesenczi/newMSM).

[^10]: [Human Connectome Project Workbench Command](https://www.humanconnectome.org/software/workbench-command) 

[^11]: [Medical Image Registration ToolKit MIRTK](http://mirtk.github.io/), needs root access for installation.

[^12]: fsl installation link (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)

[^13]: adapted from [dHCP template alignment repo](https://github.com/ecr05/dHCP_template_alignment)

