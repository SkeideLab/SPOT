import pandas as pd
import os
import numpy as np
import nibabel as nib

def get_folders_in_directory(directory):
    folders = []
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        if os.path.isdir(item_path) and not os.path.isfile(item_path):
            folders.append(item)
    return folders

file = '/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_fmri_pipeline/participants.tsv'
df = pd.read_csv(file, sep='\t')
subject_column = df[['participant_id', "birth_age"]]

base_folder = '/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_fmri_pipeline/'
seven_folders = os.listdir(base_folder)

subject_info = pd.DataFrame(columns=["sub_id", "birth_age","scan_age", "sess_id", "func_path", "xfm", "T2w", "wsr", "wsl", "psr", "psl", "msr", "msl", "isr", "isl", "vsr", "vsl", "ssr", "ssl", "mesr", "mesl"])
# Iterate over each subject ID and each of the seven folders

counter = 0
for index, row in subject_column.iterrows():
    file_name = "sub-" + row['participant_id']
    sess_dir=f"/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_fmri_pipeline/{file_name}/{file_name}_sessions.tsv"
    folder_path = os.path.join(base_folder, file_name)
    sess_if = pd.read_csv(sess_dir, sep='\t') 
    if os.path.isdir(folder_path):
        sess_ids = sess_if["session_id"]
        print(sess_ids)
        if len(sess_ids) > 1:
            for sess_id in sess_ids:
                idx = sess_if.index.get_loc(sess_if[sess_if["session_id"] == sess_id].index[0])
                print(idx)
                sess_id = "ses-" + str(sess_id)      
                #ses_row2 =sess_if[sess_if["session_id"] == int(sess_id[4:])]                    
                #print(ses_row)
                subject_info.at[index + counter, "sub_id"] = file_name
                subject_info.at[index + counter, "birth_age"] = row["birth_age"]
                subject_info.at[index + counter, "scan_age"] = sess_if.at[idx,"scan_age"]
                if subject_info.at[index + counter, "scan_age"] < 34:
                    subject_info.at[index + counter, "scan_age"] = np.nan
                subject_info.at[index + counter, "sess_id"] = sess_id
                file_path_xfm = os.path.join(folder_path, sess_id, "xfm")
                if os.path.exists(file_path_xfm):
                    for xfm_file in os.listdir(file_path_xfm):
                        file_path_xfm2 = os.path.join(file_path_xfm, xfm_file)
                        if xfm_file.endswith("bold_to-T2w_mode-image.mat"):
                            subject_info.at[index + counter, "xfm"] = file_path_xfm2
                file_path_func = os.path.join(folder_path, sess_id, "func")
                if os.path.exists(file_path_func):
                    for func_file in os.listdir(file_path_func):
                        file_path_func2 = os.path.join(file_path_func, func_file)
                        if func_file.endswith("bold.nii.gz"):
                            subject_info.at[index + counter, "func_path"] = file_path_func2
                file_path_an = "/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_anat_pipeline/"
                file_path_anat = os.path.join(file_path_an, file_name, sess_id, "anat")
                if os.path.exists(file_path_anat):
                    for anat_file in os.listdir(file_path_anat):
                        file_path_anat2 = os.path.join(file_path_anat, anat_file)
                        if anat_file.endswith("desc-restore_T2w.nii.gz"):
                            subject_info.at[index + counter, "T2w"] = file_path_anat2
                        elif anat_file.endswith("right_wm.surf.gii"):
                            subject_info.at[index + counter, "wsr"] = file_path_anat2
                        elif anat_file.endswith("left_wm.surf.gii"):
                            subject_info.at[index + counter, "wsl"] = file_path_anat2
                        elif anat_file.endswith("right_pial.surf.gii"):
                            subject_info.at[index + counter, "psr"] = file_path_anat2
                        elif anat_file.endswith("left_pial.surf.gii"):
                            subject_info.at[index + counter, "psl"] = file_path_anat2
                        elif anat_file.endswith("right_midthickness.surf.gii"):
                            subject_info.at[index + counter, "msr"] = file_path_anat2
                        elif anat_file.endswith("left_midthickness.surf.gii"):
                            subject_info.at[index + counter, "msl"] = file_path_anat2
                        elif anat_file.endswith("right_inflated.surf.gii"):
                            subject_info.at[index + counter, "isr"] = file_path_anat2
                        elif anat_file.endswith("left_inflated.surf.gii"):
                            subject_info.at[index + counter, "isl"] = file_path_anat2
                        elif anat_file.endswith("right_vinflated.surf.gii"):
                            subject_info.at[index + counter, "vsr"] = file_path_anat2
                        elif anat_file.endswith("left_vinflated.surf.gii"):
                            subject_info.at[index + counter, "vsl"] = file_path_anat2
                        elif anat_file.endswith("right_sphere.surf.gii"):
                            subject_info.at[index + counter, "ssr"] = file_path_anat2
                        elif anat_file.endswith("left_sphere.surf.gii"):
                            subject_info.at[index + counter, "ssl"] = file_path_anat2
                        elif anat_file.endswith("right_desc-medialwall_mask.shape.gii"):
                            subject_info.at[index + counter, "mesr"] = file_path_anat2
                        elif anat_file.endswith("left_desc-medialwall_mask.shape.gii"):
                            subject_info.at[index + counter, "mesl"] = file_path_anat2
                counter += 1

        elif len(sess_ids) == 1:
            ses_row = sess_if.index[0]                                  
            #print(ses_row)
            subject_info.at[index + counter, "sub_id"] = file_name
            subject_info.at[index + counter, "birth_age"] = row["birth_age"]
            subject_info.at[index + counter, "scan_age"] = sess_if.at[ses_row,"scan_age"]
            if subject_info.at[index + counter, "scan_age"] < 34:
                subject_info.at[index + counter, "scan_age"] = np.nan
            sess_id = "ses-" + str(sess_ids[0])
            subject_info.at[index + counter, "sess_id"] = sess_id
            file_path_xfm = os.path.join(folder_path, sess_id, "xfm")
            if os.path.exists(file_path_xfm):
                for xfm_file in os.listdir(file_path_xfm):
                    file_path_xfm2 = os.path.join(file_path_xfm, xfm_file)
                    if xfm_file.endswith("bold_to-T2w_mode-image.mat"):
                        subject_info.at[index + counter, "xfm"] = file_path_xfm2
            file_path_func = os.path.join(folder_path, sess_id, "func")
            if os.path.exists(file_path_func):
                for func_file in os.listdir(file_path_func):
                    file_path_func2 = os.path.join(file_path_func, func_file)
                    if func_file.endswith("bold.nii.gz"):
                        subject_info.at[index + counter, "func_path"] = file_path_func2
            file_path_an = "/data/pt_02880/Package_1225541/fmriresults01/rel3_derivatives/rel3_dhcp_anat_pipeline/"
            file_path_anat = os.path.join(file_path_an, file_name, sess_id, "anat")
            if os.path.exists(file_path_anat):
                for anat_file in os.listdir(file_path_anat):
                    file_path_anat2 = os.path.join(file_path_anat, anat_file)
                    if anat_file.endswith("desc-restore_T2w.nii.gz"):
                        subject_info.at[index + counter, "T2w"] = file_path_anat2
                    elif anat_file.endswith("right_wm.surf.gii"):
                        subject_info.at[index + counter, "wsr"] = file_path_anat2
                    elif anat_file.endswith("left_wm.surf.gii"):
                        subject_info.at[index + counter, "wsl"] = file_path_anat2
                    elif anat_file.endswith("right_pial.surf.gii"):
                        subject_info.at[index + counter, "psr"] = file_path_anat2
                    elif anat_file.endswith("left_pial.surf.gii"):
                        subject_info.at[index + counter, "psl"] = file_path_anat2
                    elif anat_file.endswith("right_midthickness.surf.gii"):
                        subject_info.at[index + counter, "msr"] = file_path_anat2
                    elif anat_file.endswith("left_midthickness.surf.gii"):
                        subject_info.at[index + counter, "msl"] = file_path_anat2
                    elif anat_file.endswith("right_inflated.surf.gii"):
                        subject_info.at[index + counter, "isr"] = file_path_anat2
                    elif anat_file.endswith("left_inflated.surf.gii"):
                        subject_info.at[index + counter, "isl"] = file_path_anat2
                    elif anat_file.endswith("right_vinflated.surf.gii"):
                        subject_info.at[index + counter, "vsr"] = file_path_anat2
                    elif anat_file.endswith("left_vinflated.surf.gii"):
                        subject_info.at[index + counter, "vsl"] = file_path_anat2
                    elif anat_file.endswith("right_sphere.surf.gii"):
                        subject_info.at[index + counter, "ssr"] = file_path_anat2
                    elif anat_file.endswith("left_sphere.surf.gii"):
                        subject_info.at[index + counter, "ssl"] = file_path_anat2
                    elif anat_file.endswith("right_desc-medialwall_mask.shape.gii"):
                        subject_info.at[index + counter, "mesr"] = file_path_anat2
                    elif anat_file.endswith("left_desc-medialwall_mask.shape.gii"):
                        subject_info.at[index + counter, "mesl"] = file_path_anat2   
           
subject_info.dropna(subset=['func_path'], inplace=True)
subject_info.dropna(subset=['T2w'], inplace=True)
subject_info.dropna(subset=['vsr'], inplace=True)
subject_info.dropna(subset=['xfm'], inplace=True)
subject_info.dropna(subset=['wsr'], inplace=True)
subject_info.dropna(subset=['scan_age'], inplace=True)
subject_info.reset_index(drop=True, inplace=True)
nan_locations = np.where(subject_info.isna())
print("Row indices with NaN values:", nan_locations[0])
print("Column indices with NaN values:", nan_locations[1])
#check the number of vertices are same in pial and wm
for index, row in subject_info.iterrows():
    white_surf_right=nib.load(subject_info.at[index, "wsr"])
    white_surf_left=nib.load(subject_info.at[index, "wsl"])
    pial_surf_right=nib.load(subject_info.at[index, "psr"])
    pial_surf_left=nib.load(subject_info.at[index, "psl"])
    co_wr, tri_wr=white_surf_right.agg_data()
    co_wl, tri_wl=white_surf_left.agg_data()
    co_pr, tri_pr=pial_surf_right.agg_data()
    co_pl, tri_pl=pial_surf_left.agg_data()
    if not np.array_equal(tri_pr, tri_wr):
        subject_info.at[index, "wsr"] = "INVALID"
        #print(subject_info.at[index, "wsr"])
    elif not np.array_equal(tri_pl, tri_wl):
        subject_info.at[index, "wsl"] = "INVALID"
        #print(subject_info.at[index, "wsl"])

subject_info = subject_info[subject_info["wsr"] != "INVALID"]
subject_info = subject_info[subject_info["wsl"] != "INVALID"]
subject_info.reset_index(drop=True, inplace=True)

print(subject_info)            # Perform actions with the folder as needed

subject_info.to_csv("/data/p_02915/SPOT/dhcp_subj_path_SPOT.csv")