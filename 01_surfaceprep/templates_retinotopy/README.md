# Maps of the visual system from Wang et al.[^1] and Benson et al.[^4]

Two retinotopy templates/atlases with information about areas in the visual system and retinotopy parameters are used here, labelled as `wang2015`, `benson2014`. 

The script `transform_templatemgz_to_gifti.py` was used to create the template files here from neuropythy `.mgz` files [^5]. Basis were the atlases downloaded from [^3].

All resulting maps are based on the 164k fsaverage mesh. 

## Contents of the directory

### Templates:

- `hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii`: `wang2015` maximum probability visual area label map 
- `hemi-{hemi}_space-fsaverage_dens-164k_desc-angleretinotbenson2014_seg.shape.gii`: `benson2014` polar angle map
- `hemi-{hemi}_space-fsaverage_dens-164k_desc-eccentretinotbenson2014_seg.shape.gii`: `benson2014` eccentricity map
- `hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii`: `benson2014` visual area segmentation

### Scripts and helper files:

- `transform_templatemgz_to_gifti.py`: Transforms Wang and Benson atlas from .mgz format to gifti file format. Uses atlases downloaded from [^3].
- `labels_visualareas_<template>.txt`[^2]: region labels as used in the atlas and corresponding region names
- `labeltable_<template>`: Freesurfer look-up tables. File with region labels, corresponding names and randomly chosen colors. Created manually from roi labels.





[^1]: Wang, L., Mruczek, R. E. B., Arcaro, M. J. & Kastner, S. Probabilistic Maps of Visual Topography in Human Cortex. Cerebral Cortex 25, 3911–3931 (2015).

[^2]: wang2015 roi labels: downloaded from https://napl.scholar.princeton.edu/resources. Last Update: 12/9/14. Probability Atlas v4. 
benson2014 roi labels: copied from https://github.com/noahbenson/neuropythy/wiki/Retinotopy#using-neuropythyvisionpredict_retinotopy

[^3]: retinotopy atlases on fsaverage: https://github.com/noahbenson/neuropythy/blob/master/neuropythy/lib/data/fsaverage/surf/ 

[^4]: Benson, Noah C., Omar H. Butt, David H. Brainard, und Geoffrey K. Aguirre. „Correction of Distortion in Flattened Representations of the Cortical Surface Allows Prediction of V1-V3 Functional Organization from Anatomy“. PLOS Computational Biology 10, Nr. 3 (27. März 2014): e1003538. https://doi.org/10.1371/journal.pcbi.1003538.

[^5]: Benson NC, Winawer J (2018) Bayesian Analysis of Retinotopic Maps. bioRxiv doi:10.1101/325597. https://github.com/noahbenson/neuropythy
