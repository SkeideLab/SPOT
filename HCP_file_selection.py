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



base_folder = '/data/pt_02880/HCP_D/fmriresults01/'

seven_folders = os.listdir(base_folder)
seven_folders.remove("manifests")
subject_completeness = pd.read_csv("/data/pt_02880/HCP_D/HCD_LS_2.0_subject_completeness.csv", skiprows=[1])
vision1 = pd.read_csv("/data/p_02915/dhcp_derivatives_SPOT/HCP-D/asr01.csv", skiprows=[1])
vision2 = pd.read_csv("/data/p_02915/dhcp_derivatives_SPOT/HCP-D/ysr01.csv", skiprows=[1])
vision3 = pd.read_csv("/data/p_02915/dhcp_derivatives_SPOT/HCP-D/cbcl01.csv", skiprows=[1])
vision4 = pd.read_csv("/data/p_02915/dhcp_derivatives_SPOT/HCP-D/cbcl1_501.csv", skiprows=[1])
subject_info = pd.DataFrame(columns=["sub_id","scan_age", "sex","wsr", "wsl", "psr", "psl", "msr", "msl", "isr", "isl", "vsr", "vsl", "ssr", "ssl", "fwsr", "fwsl", "fpsr", "fpsl", "fmsr", "fmsl", "fisr", "fisl", "fvsr", "fvsl", "fssr", "fssl"])
# Iterate over each subject ID and each of the seven folders
counter = 0
visual = 0
failed_QC=0
for index, row in enumerate(seven_folders):
    file_name = row    
    folder_path = os.path.join(base_folder, file_name)
    if os.path.isdir(folder_path):
        print(file_name[:10])
        idx = subject_completeness.index.get_loc(subject_completeness[subject_completeness["src_subject_id"] == file_name[:10]].index[0])
        if file_name[:10] in vision1["src_subject_id"].values:
            vision = vision1
            idx2 = vision.index.get_loc(vision[vision["src_subject_id"] == file_name[:10]].index[0])
            value = "asr10_6"
            parameter = 1
        elif file_name[:10] in vision2["src_subject_id"].values:
            vision = vision2
            idx2 = vision.index.get_loc(vision[vision["src_subject_id"] == file_name[:10]].index[0])
            value = "ysr_59"
            parameter = 0
        elif file_name[:10] in vision3["src_subject_id"].values:
            vision = vision3
            idx2 = vision.index.get_loc(vision[vision["src_subject_id"] == file_name[:10]].index[0])
            value = "cbcl56d"
            parameter = 0
        elif file_name[:10] in vision4["src_subject_id"].values:
            vision = vision4
            idx2 = vision.index.get_loc(vision[vision["src_subject_id"] == file_name[:10]].index[0])
            value = "cbcl56d"
            parameter = 0

        if (subject_completeness.at[idx,"RS-fMRI_Count"] >= 2 and
            pd.isna(subject_completeness.at[idx,"QC_Issue_Codes"]) and
            vision.at[idx2,value] == parameter):            
            subject_info.at[index + counter, "sub_id"] = file_name
            subject_info.at[index + counter, "sex"] = subject_completeness.at[idx,"sex"]          
            subject_info.at[index + counter, "scan_age"] = subject_completeness.at[idx,"interview_age"]
            if subject_info.at[index + counter, "scan_age"] < 18*12:
                    subject_info.at[index + counter, "scan_age"] = np.nan          
            file_path_anat_native = os.path.join(base_folder, file_name, "MNINonLinear","Native")
            file_path_anat_fsaverage = os.path.join(base_folder, file_name, "MNINonLinear","fsaverage_LR32k")
            for anat_file in os.listdir(file_path_anat_native):
                file_path_anat2 = os.path.join(file_path_anat_native, anat_file)
                if anat_file.endswith("R.white.native.surf.gii"):
                    subject_info.at[index + counter, "wsr"] = file_path_anat2
                elif anat_file.endswith("L.white.native.surf.gii"):
                    subject_info.at[index + counter, "wsl"] = file_path_anat2
                elif anat_file.endswith("R.pial.native.surf.gii"):
                    subject_info.at[index + counter, "psr"] = file_path_anat2
                elif anat_file.endswith("L.pial.native.surf.gii"):
                    subject_info.at[index + counter, "psl"] = file_path_anat2
                elif anat_file.endswith("R.midthickness.native.surf.gii"):
                    subject_info.at[index + counter, "msr"] = file_path_anat2
                elif anat_file.endswith("L.midthickness.native.surf.gii"):
                    subject_info.at[index + counter, "msl"] = file_path_anat2
                elif anat_file.endswith("R.inflated.native.surf.gii"):
                    subject_info.at[index + counter, "isr"] = file_path_anat2
                elif anat_file.endswith("L.inflated.native.surf.gii"):
                    subject_info.at[index + counter, "isl"] = file_path_anat2
                elif anat_file.endswith("R.very_inflated.native.surf.gii"):
                    subject_info.at[index + counter, "vsr"] = file_path_anat2
                elif anat_file.endswith("L.very_inflated.native.surf.gii"):
                    subject_info.at[index + counter, "vsl"] = file_path_anat2
                elif anat_file.endswith("R.sphere.native.surf.gii"):
                    subject_info.at[index + counter, "ssr"] = file_path_anat2
                elif anat_file.endswith("L.sphere.native.surf.gii"):
                    subject_info.at[index + counter, "ssl"] = file_path_anat2
            for anat_file in os.listdir(file_path_anat_fsaverage):
                file_path_anat2 = os.path.join(file_path_anat_fsaverage, anat_file)
                if anat_file.endswith("R.white_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fwsr"] = file_path_anat2
                elif anat_file.endswith("L.white_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fwsl"] = file_path_anat2
                elif anat_file.endswith("R.pial_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fpsr"] = file_path_anat2
                elif anat_file.endswith("L.pial_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fpsl"] = file_path_anat2
                elif anat_file.endswith("R.midthickness_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fmsr"] = file_path_anat2
                elif anat_file.endswith("L.midthickness_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fmsl"] = file_path_anat2
                elif anat_file.endswith("R.inflated_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fisr"] = file_path_anat2
                elif anat_file.endswith("L.inflated_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fisl"] = file_path_anat2
                elif anat_file.endswith("R.very_inflated_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fvsr"] = file_path_anat2
                elif anat_file.endswith("L.very_inflated_MSMAll.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fvsl"] = file_path_anat2
                elif anat_file.endswith("R.sphere.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fssr"] = file_path_anat2
                elif anat_file.endswith("L.sphere.32k_fs_LR.surf.gii"):
                    subject_info.at[index + counter, "fssl"] = file_path_anat2
            counter += 1
        elif (subject_completeness.at[idx,"RS-fMRI_Count"] < 2 or
            not pd.isna(subject_completeness.at[idx,"QC_Issue_Codes"]) and
            subject_completeness.at[idx,"interview_age"] >= 18*12):
            failed_QC = failed_QC +1
             
        elif (vision.at[idx2,value] != parameter and subject_completeness.at[idx,"interview_age"] >= 17*12):
            visual = visual + 1

print(f"failed QC: {failed_QC}")
print(f"visual impairment: {visual}")

subject_info.dropna(subset=['scan_age'], inplace=True)
subject_info.reset_index(drop=True, inplace=True)
nan_locations = np.where(subject_info.isna())
print("Row indices with NaN values:", nan_locations[0])
print("Column indices with NaN values:", nan_locations[1])

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

for index, row in subject_info.iterrows():
    white_surf_right=nib.load(subject_info.at[index, "fwsr"])
    white_surf_left=nib.load(subject_info.at[index, "fwsl"])
    pial_surf_right=nib.load(subject_info.at[index, "fpsr"])
    pial_surf_left=nib.load(subject_info.at[index, "fpsl"])
    co_wr, tri_wr=white_surf_right.agg_data()
    co_wl, tri_wl=white_surf_left.agg_data()
    co_pr, tri_pr=pial_surf_right.agg_data()
    co_pl, tri_pl=pial_surf_left.agg_data()
    if not np.array_equal(tri_pr, tri_wr):
        subject_info.at[index, "fwsr"] = "INVALID"
        #print(subject_info.at[index, "wsr"])
    elif not np.array_equal(tri_pl, tri_wl):
        subject_info.at[index, "fwsl"] = "INVALID"
        #print(subject_info.at[index, "wsl"])

subject_info = subject_info[subject_info["fwsr"] != "INVALID"]
subject_info = subject_info[subject_info["fwsl"] != "INVALID"]
subject_info.reset_index(drop=True, inplace=True)
subject_info = subject_info.sort_values(by="scan_age").reset_index(drop=True)

print(subject_info)            # Perform actions with the folder as needed

#subject_info.to_csv("/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv")