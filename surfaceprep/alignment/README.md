# How to align an individual anatomical surface to a group surface

adapted from [dHCP template alignment repo](https://github.com/ecr05/dHCP_template_alignment)

## Prerequisites

- dHCP second release data[^1]
- subject label and session label of subject that is to be aligned
- subject age in months (can be extracted from .tsv file)
- transformation from individual anatomical volume space to group anatomical volume space[^2]
- 40 weeks anatomical template in volume space (must be the one that corresponds to the transformation)[^3]
- target surface template[^4]
- newMSM installation[^5]
- workbench command installation[^6]
- MIRTK installation[^7]

[^1]: [Developing Human Connectome Project](http://www.developingconnectome.org/data-release/second-data-release/)

[^2]: included in the dHCP data

[^3]: The group volumetric atlas used in the dHCP second release is called ['dhcp-volumetric-atlas-groupwise'](https://gin.g-node.org/BioMedIA/dhcp-volumetric-atlas-groupwise). Alternatively, the same atlas can be downloaded from Sean Fitzgibbon's ['augmented version'](https://git.fmrib.ox.ac.uk/seanf/dhcp-resources/-/blob/master/docs/dhcp-augmented-volumetric-atlas.md).

[^4]: The dHCP project has published the [dhcpSym atlas](https://brain-development.org/brain-atlases/atlases-from-the-dhcp-project/cortical-surface-template/), which can be used here.

[^5]: newMSM runs on FSL and can be installed from Renato Besenczi's github repo [newMSM](https://github.com/rbesenczi/newMSM).

[^6]: [Human Connectome Project Workbench Command](https://www.humanconnectome.org/software/workbench-command) 

[^7]: [Medical Image Registration ToolKit MIRTK](http://mirtk.github.io/), needs root access for installation.