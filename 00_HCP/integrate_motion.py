import pandas as pd
import os

subject_info = pd.read_csv('/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_subj_path_SPOT_old.csv')
#sub_num = int(sys.argv[1])

# Function to extract last six columns from a line
def extract_last_six_columns(line):
    columns = line.split()
    last_six_columns = columns[6:12]
    return " ".join(last_six_columns)

for index, row in subject_info.iterrows():
    # List to store extracted data
    extracted_data = []

    sub_id = subject_info.at[index, "sub_id"]
    sub = sub_id.replace('sub-','')
    input_file1=(f"/data/pt_02880/HCP_D/fmriresults01/{sub}/MNINonLinear/Results/rfMRI_REST1_AP/Movement_Regressors_hp0_clean.txt")
    input_file2=(f"/data/pt_02880/HCP_D/fmriresults01/{sub}/MNINonLinear/Results/rfMRI_REST1_PA/Movement_Regressors_hp0_clean.txt")
    input_file3=(f"/data/pt_02880/HCP_D/fmriresults01/{sub}/MNINonLinear/Results/rfMRI_REST2_AP/Movement_Regressors_hp0_clean.txt")
    input_file4=(f"/data/pt_02880/HCP_D/fmriresults01/{sub}/MNINonLinear/Results/rfMRI_REST3_PA/Movement_Regressors_hp0_clean.txt")

    for file_path in [input_file1, input_file2, input_file3, input_file4]:
        if os.path.isfile(file_path):  # Check if the file exists
            with open(file_path, 'r') as file:
                # Read each line from the file
                for line in file:
                    # Extract last six columns from the line
                    extracted_data.append(extract_last_six_columns(line))

    output_file_path=(f"/data/p_02915/dhcp_derivatives_SPOT/HCP-D/hcp_surface/{sub}/func/{sub}_motion_derivatives.txt")
    # Write extracted data to a new file
    
    try:
        with open(output_file_path, 'w') as output_file:
            for data in extracted_data:
                output_file.write(data + "\n")
        print(f"File successfully saved at {output_file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")
