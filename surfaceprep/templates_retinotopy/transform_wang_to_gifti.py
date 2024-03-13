"""Transform the Wang 2015 retinotopy template as published by neuropythy into a usable gifti label file.
The template file has the extension '.mgz', which is a freesurfer format for volume files. However, we know
that the template is on the surface, not in volumetric space. The shape of the data array when read 
in with nibabel is (1, 1, 163842), 
so the third dimension are the fsaverage 164k nodes. The array is transformed to 1 dimension and saved as 
a gifti '.gii' image. Then, a labeltable with the correct region labels is added. 
Needs a copy of neuropythy source code.
Left and right input_wangtemplate need to be downloaded from https://github.com/noahbenson/neuropythy/blob/master/neuropythy/lib/data/fsaverage/surf/{l,r}h.wang15_mplbl.v1_0.mgz
Needs workbench command to be callable from command line (https://www.humanconnectome.org/software/workbench-command).
"""

import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

labeltable = Path(__file__).parent.resolve() / "labeltable"
print(f"Using labeltable: {labeltable}")
templates = [
    (
        "/data/u_kieslinger_software/code/neuropythy/neuropythy/lib/data/fsaverage/surf/{hemi}h.wang15_mplbl.v1_0.mgz",
        "hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii",
        np.int32,
    ),
    (
        "/data/u_kieslinger_software/code/neuropythy/neuropythy/lib/data/fsaverage/surf/{hemi}h.benson14_eccen.v4_0.mgz",
        "hemi-{hemi}_space-fsaverage_dens-164k_desc-eccentretinotbenson2014_seg.shape.gii",
        np.float32,
    ),
    (
        "/data/u_kieslinger_software/code/neuropythy/neuropythy/lib/data/fsaverage/surf/{hemi}h.benson14_angle.v4_0.mgz",
        "hemi-{hemi}_space-fsaverage_dens-164k_desc-angleretinotbenson2014_seg.shape.gii",
        np.float32,
    ),
]

for template_paths in templates:
    for hemi in ["l", "r"]:
        ##----------------------CHANGE PATHS--------------------------------
        input_template = template_paths[0].format(hemi=hemi)
        output_template = Path(__file__).parent.resolve() / template_paths[1].format(
            hemi=hemi.upper()
        )
        dtype_templ = template_paths[2]
        ##--------------------------------------------------------------------------------

        Path(output_template).parent.mkdir(
            exist_ok=True,
            parents=True,
        )

        if not output_template.exists():

            if dtype_templ == np.int32:
                # save into intermediate file without label table
                intermediate_metric_name = str(output_template).split(".", maxsplit=1)
                intermediate_metric_name = (
                    intermediate_metric_name[0]
                    + "_metric"
                    + intermediate_metric_name[1]
                )
            else:
                intermediate_metric_name = str(output_template)

            # import data from .mgz and save as gifti
            img_mgz = nib.load(input_template)
            templ_data = img_mgz.get_fdata()
            print(f"Original template data shape: {templ_data.shape}")

            img_gifti = nib.gifti.GiftiImage(
                darrays=[nib.gifti.GiftiDataArray(dtype_templ(np.squeeze(templ_data)))]
            )
            nib.save(img_gifti, intermediate_metric_name)

            if dtype_templ == np.int32:
                # transform into label file with label table
                proc = subprocess.run(
                    [
                        "wb_command",
                        "-metric-label-import",
                        intermediate_metric_name,
                        labeltable,
                        output_template,
                    ],
                    check=False,
                    capture_output=True,
                    text=True,
                )

                assert (
                    proc.returncode == 0
                ), f"Error while trying to convert metric file to label file using wb_command. Error: {proc.stderr}"

                # remove intermediate file
                subprocess.run(
                    ["rm", intermediate_metric_name],
                    text=True,
                    check=True,
                    capture_output=True,
                )
