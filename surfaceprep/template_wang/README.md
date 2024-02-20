# Maps of the visual system from Wang et al.[^1]

- `hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii`: maximum probability map on fsaverage. Created by `transform_wang_to_gifti.py`. Can be used directly.
- `transform_wang_to_gifti.py`: Transforms Wang atlas from .mgz format to gifti file format. Uses atlas downloaded from [^3].
- `ROIfiles_Labeling.txt`[^2]: region labels as used in the atlas and corresponding region names
- `labeltable`: Freesurfer look-up table. File with region labels, corresponding names and randomly chosen colors. Created manually





[^1]: Wang, L., Mruczek, R. E. B., Arcaro, M. J. & Kastner, S. Probabilistic Maps of Visual Topography in Human Cortex. Cerebral Cortex 25, 3911–3931 (2015).

[^2]: downloaded from https://napl.scholar.princeton.edu/resources. Last Update: 12/9/14. Probability Atlas v4

[^3]: https://github.com/noahbenson/neuropythy/blob/master/neuropythy/lib/data/fsaverage/surf/rh.wang15_mplbl.v1_0.mgz, https://github.com/noahbenson/neuropythy/blob/master/neuropythy/lib/data/fsaverage/surf/lh.wang15_mplbl.v1_0.mgz