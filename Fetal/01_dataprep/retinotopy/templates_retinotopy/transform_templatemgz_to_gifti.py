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

parentdir = Path(__file__).parent.resolve()
labeltable_wang = parentdir / "labeltable_wang2015"
labeltable_benson = parentdir / "labeltable_benson2014"
print(f"Using labeltables: {labeltable_wang} and {labeltable_benson}.")
templates = [
    (
        "/data/u_kieslinger_software/code/neuropythy/neuropythy/lib/data/fsaverage/surf/{hemi}h.wang15_mplbl.v1_0.mgz",
        "hemi-{hemi}_space-fsaverage_dens-164k_desc-visualtopographywang2015_label-maxprob_seg.label.gii",
        np.int32,
        labeltable_wang,
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
    (
        "/data/u_kieslinger_software/code/neuropythy/neuropythy/lib/data/fsaverage/surf/{hemi}h.benson14_varea.v4_0.mgz",
        "hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii",
        np.int32,
        labeltable_benson,
    ),
]

for template_paths in templates:
    for hemi in ["l", "r"]:

        # fill paths for this hemi
        input_template = template_paths[0].format(hemi=hemi)
        output_template = Path(__file__).parent.resolve() / template_paths[1].format(
            hemi=hemi.upper()
        )

        # we need template datatype to be able to cast the array to fit gifti standard
        dtype_templ = template_paths[2]

        Path(output_template).parent.mkdir(
            exist_ok=True,
            parents=True,
        )

        if not output_template.exists():

            # import data from .mgz and save as gifti
            img_mgz = nib.load(input_template)
            templ_data = img_mgz.get_fdata()
            print(f"Original template data shape: {templ_data.shape}")

            img_gifti = nib.gifti.GiftiImage(
                darrays=[nib.gifti.GiftiDataArray(dtype_templ(np.squeeze(templ_data)))]
            )

            if dtype_templ != np.int32:
                nib.save(img_gifti, str(output_template))
            else:
                # save into intermediate file without label table
                # needs conversion into proper 'label' file with label table later
                intermediate_metric_name = str(output_template).split(".", maxsplit=1)
                intermediate_metric_name = (
                    intermediate_metric_name[0]
                    + "_metric"
                    + intermediate_metric_name[1]
                )

                nib.save(img_gifti, intermediate_metric_name)

                assert (
                    len(template_paths) == 4
                ), "ERROR: Missing labeltable (must be 4. entry)"

                # transform into label file with label table
                proc = subprocess.run(
                    [
                        "wb_command",
                        "-metric-label-import",
                        intermediate_metric_name,
                        template_paths[3],
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
                Path(intermediate_metric_name).unlink()
            print(f"Saving template {input_template} at {output_template}...")
        else:
            print(
                f"Template {input_template} already exists at {output_template}. Skipping this..."
            )
