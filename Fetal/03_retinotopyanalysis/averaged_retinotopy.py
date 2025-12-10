"""
Average ccf results
"""

import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np

def get_indices_roi(labels_area, visparc):
    """Get indices of vertices in gifti image that are located in a specific region of interest.

    Args:
        labels_area (list): labels for the region of interest as used in the parcellation
        visparc (nibabel.gifti.GiftiImage): parcellation of brain surface into regions, using labels

    Returns:
        numpy.array: n_vertices_roi. Indices of vertices that lie in the ROI.
    """
     # Ensure labels_area is a list
    if not isinstance(labels_area, list):
        labels_area = [labels_area]
    
    # Collect indices for all labels in labels_area
    indices_area = np.concatenate([
        np.nonzero(visparc.agg_data() == label)[0]
        for label in labels_area
    ])

    return indices_area

VISPARC_PATH = (
    "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k_desc-retinotbenson2014_label-visarea_seg.label.gii")
LABELS_V2 = [2, 3]
 
PREFIX_MODEL = (
    "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel_var1/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k"
)
INPUT = (
    "{prefix_model}_label-{param}_desc-{model}_roi-v2th00_metric.gii"
)
#INPUT = (
#    "{prefix_model}_desc-{model}_{param}.gii"
#)
OUTPUT = (
    "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel_var1/Averaged_younger_fetal_{hemi}_label-{param}_desc-{model}_roi-v2th00_metric.gii"
)
# OUTPUT = (
#    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_younger_fetal_{hemi}_desc-{model}_{param}.gii"
# )
# INPUT = (
#    "{prefix_model}_desc-real_{parm}.gii"
# )
# OUTPUT = (
#    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/Averaged_{hemi}_desc-real_{parm}.gii"
# )
subject_info = pd.read_csv(
    '/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_young.csv')
sub_num = len(subject_info["sub_id"])

sum_data = []  # Initialize outside the loop
for param in ["eccentricity", "polarangle"]:
#for param in ["sigma"]:
    for hemi in ["L", "R"]:
        # FORMAT PATHS FOR INPUT AND OUTPUT
        temp = nib.load(VISPARC_PATH.format(hemi=hemi))
        indices = get_indices_roi(LABELS_V2, temp)

        for model in ["real"]:
            sum_data = []
            for index, row in subject_info.iterrows():
                sub_id = subject_info.at[index, "sub_id"]
                sess_id = subject_info.at[index, "sess_id"]
                sub = sub_id.replace('sub-', '')
                ses = sess_id.replace('ses-', '')
                prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                output_path = OUTPUT.format(
                    prefix_model=prefix_model, hemi=hemi, model=model, param=param
                )
                input_path = INPUT.format(
                    prefix_model=prefix_model, model=model, param=param
                )
                try:
                    ccf = surface.load_surf_data(
                        input_path).astype(np.float64)
                    # if hemi == "L" and param == "polarangle":
                    #    # Apply transformations
                    #    ccf[(ccf > 90) & (ccf < 270)] = 0
                    #    ccf[(ccf > 0) & (ccf <= 90)] = ccf[(ccf > 0) & (ccf <= 90)] + 90
                    #    ccf[(ccf >= 270) & (ccf <= 360)] = ccf[(ccf >= 270) & (ccf <= 360)] - 270
                    #    #print(ccf[(ccf > 0) & (ccf <= 90)])
                    # elif hemi == "R" and param == "polarangle":
                    #    ccf[(ccf >= 0) & (ccf < 90)] = 0
                    #    ccf[(ccf > 270) & (ccf <= 360)] = 0
                    #    ccf[(ccf >= 90) & (ccf <= 270)] = 270 - ccf[(ccf >= 90) & (ccf <= 270)]

                    sum_data.append(ccf)

                except ValueError:
                    print(
                        f"No model results found under {INPUT.format(prefix_model=prefix_model, model=model, param = param)}."
                    )
                    print(" Skipping this...", flush=True)
                    continue
            sum_data_array = np.array(sum_data)
            #sum_data_array = np.where(sum_data_array == 0, np.nan, sum_data_array)
            averaged_data = np.nanmean(sum_data_array, axis=0)
            # save as gifti
            img_gifti = nib.gifti.GiftiImage(
                darrays=[
                    nib.gifti.GiftiDataArray(
                        np.float32(np.squeeze(averaged_data)))
                ]
            )

            nib.save(img_gifti, output_path)
