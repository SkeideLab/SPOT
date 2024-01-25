"""Transform the Wang 2015 retinotopy template as published by neuropythy into a usable gifti file.
The template file has the extension '.mgz', which is a freesurfer format for volume files. However, we know
that the template is on the surface, not in volumetric space. The shape of the data array when read 
in with nibabel is (1, 1, 163842), 
so the third dimension are the fsaverage 164k nodes. The array is transformed to 1 dimension and saved as 
a gifti '.gii' image.
"""
import nibabel as nib
import numpy as np

##----------------------CHANGE PATHS--------------------------------
input_wangtemplate = "/data/u_kieslinger_software/code/neuropythy/neuropythy/lib/data/fsaverage/surf/rh.wang15_mplbl.v1_0.mgz"
output_wangtemplate = "/data/p_02495/templates/template_visualtopography_wang2015/hemi-R_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.gii"
##--------------------------------------------------------------------------------

img_wang_mgz = nib.load(input_wangtemplate)
wang_data = img_wang_mgz.get_fdata()
print(wang_data.shape)

img_wang = nib.gifti.GiftiImage(
    darrays=[nib.gifti.GiftiDataArray(np.int32(np.squeeze(wang_data)))]
)
nib.save(img_wang, output_wangtemplate)
