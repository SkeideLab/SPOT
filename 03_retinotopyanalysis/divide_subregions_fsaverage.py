import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
PREFIX_ALLDATA = "/data/p_02915/SPOT/01_dataprep/retinotopy/templates_retinotopy/hemi-{hemi}_space-fsaverage_dens-164k"

BENSON_ECC = "_desc-eccentretinotbenson2014_seg.shape.gii"
BENSON_ANGLE = "_desc-angleretinotbenson2014_seg.shape.gii"
BENSON_VISAREAS = "_desc-retinotbenson2014_label-visarea_seg.label.gii"

OUTSUFFIX_FOVEAL = "_desc-foveal_roi.shape.gii"
OUTSUFFIX_MIDDLE = "_desc-middle_roi.shape.gii"
OUTSUFFIX_PERIPHERAL = "_desc-peripheral_roi.shape.gii"
OUTSUFFIX_DORSAL = "_desc-dorsal_roi.shape.gii"
OUTSUFFIX_VENTRAL = "_desc-ventral_roi.shape.gii"

PREFIX_MODEL = (
    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k"
)
#sub-CC00056XX07_ses-10700_hemi-L_mesh-fsaverage_dens-164k_label-polarangle_desc-real_roi-v2th00_metric
INPUT = (
    "{prefix_model}_label-{parm}_desc-{model}_roi-v2th00_metric.gii"
)
OUTPUT = (
    "{prefix_model}_{hemi}_label-{parm}_desc-{model}_v2th00"
)
subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT.csv')
sub_num = len(subject_info["sub_id"])

for hemi in ["L", "R"]:
    boundaries = {"L": (3.225, 12.42), "R": (3.28, 15.5)}
    prefix_thisdata = PREFIX_ALLDATA.format(
        hemi=hemi,
    )
    ecc = surface.load_surf_data(prefix_thisdata + BENSON_ECC)
    pangle = surface.load_surf_data(prefix_thisdata + BENSON_ANGLE)
    visareas = surface.load_surf_data(prefix_thisdata + BENSON_VISAREAS)

    # restrict all masks to only V2 vertices
    mask_v2 = visareas == 2 # | (visareas == 3)  # 2 is code for V2
    mask_foveal = (ecc > 0) & (ecc <= boundaries[hemi][0]) & mask_v2
    mask_middle = (
        (ecc > boundaries[hemi][0]) & (ecc <= boundaries[hemi][1]) & mask_v2
    )
    mask_peripheral = (ecc > boundaries[hemi][1]) & mask_v2

    # TODO the masks have very different numbers of vertices (359, 977,11996).
    # Is this because there is no restriction to V2? change boundaries?
    # No, restriction to V2/V2+V3 does not change anything.
    # The reason must be different numbers of vertices on native vs fsaverage mesh in these regions.
    # Change the boundaries so that the numbers are similar? not sure
    for region, region_mask in [
                    ("foveal", mask_foveal),
                    ("middle", mask_middle),
                    ("peripheral", mask_peripheral),
                ]:
        print(f"Size of {region} mask is: {region_mask.sum()} vertices.")

    mask_ventral = (pangle > 90) & mask_v2
    mask_dorsal = (pangle < 90) & (pangle > 0) & mask_v2


    for param in ["eccentricity", "polarangle"]:              
        # FORMAT PATHS FOR INPUT AND OUTPUT     
           
        for model in ["real", "simulated"]:
            sum_data=[]    
            for index, row in subject_info.iterrows(): 
                sub_id = subject_info.at[index, "sub_id"]
                sess_id = subject_info.at[index, "sess_id"]
                sub = sub_id.replace('sub-','')
                ses = sess_id.replace('ses-','')
                prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                output_path = OUTPUT.format(
                    prefix_model=prefix_model, hemi=hemi, model=model, parm=param
                )
                input_path =  INPUT.format(
                    prefix_model=prefix_model, model=model, parm=param
                )
                ccf_v0 = surface.load_surf_data(input_path).astype(np.float64)
                
    # assert (
    #     np.abs(np.sum(mask_foveal) - np.sum(mask_middle)) < 10
    # ), "Middle and foveal mask have a difference of more than 10 vertices"
    # assert (
    #     np.abs(np.sum(mask_peripheral) - np.sum(mask_middle)) < 10
    # ), "Middle and peripheral mask have a difference of more than 10 vertices"
                for data, outpath in [
                (mask_foveal, OUTSUFFIX_FOVEAL),
                [mask_middle, OUTSUFFIX_MIDDLE],
                (mask_peripheral, OUTSUFFIX_PERIPHERAL),
                (mask_ventral, OUTSUFFIX_VENTRAL),
                (mask_dorsal, OUTSUFFIX_DORSAL),
            ]:
                    ccf_sub = ccf_v0[data].astype(np.float64)
                    img = nib.gifti.GiftiImage(
                        darrays=[nib.gifti.GiftiDataArray(np.int32(ccf_sub))]
                    )
                    nib.save(img, output_path + outpath)
                    # print(f"Saving mask under {prefix_thisdata+outpath}...")


            
                

                
                

