import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import chisquare
from scipy.stats import bootstrap   
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

def chi_squared_test(dist1, dist2):
    dist1 = np.histogram(dist1,bins=range(20))[0] +0.5   
    dist2 = np.histogram(dist2,bins=range(20))[0] +0.5
    chi2, p, _, _ = chi2_contingency([dist1,dist2])
    return chi2, p

PREFIX_MODEL = (
    "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_"
)
PREFIX_MODEL_1 = (
    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_"
)

PARAM = "desc-real_{param}.gii"
#PARAM = "label-{param}_desc-real_roi-v2th00_metric.gii"

subject_34w = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_2.csv')
subject_30w = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
sub_num_30w = len(subject_30w["sub_id"])
sub_num_34w = len(subject_34w["sub_id"])
print("fetal (younger) vs. neonates over 37")
for hemi in ["L", "R"]:              
    for model in ["34w", "neonates"]:
        sum_data_r=[]
        parameters_r = []   
        sum_data_sig=[]
        parameters_sig = []           
        if model =="34w":
            subject_info = subject_34w
            sub_num = sub_num_34w
        #elif model =="30w":        
        elif model == "neonates":
            subject_info = subject_30w
            sub_num = sub_num_30w
        for param in ["r", "sigma"]:                   
        #for param in ["eccentricity", "polarangle"]:
        # FORMAT PATHS FOR INPUT AND OUTPUT 
            for index, row in subject_info.iterrows(): 
                sub_id = subject_info.at[index, "sub_id"]
                sess_id = subject_info.at[index, "sess_id"]
                sub = sub_id.replace('sub-','')
                ses = sess_id.replace('ses-','')
                if model == "34w":
                    prefix_model = PREFIX_MODEL.format(
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                    )
                elif model == "neonates":
                    prefix_model = PREFIX_MODEL_1.format(
                        sub=sub,
                        ses=ses,
                        hemi=hemi,
                    )
                pram_model = PARAM.format(
                    param=param,
                )
                input_path = prefix_model + pram_model
                
                ccf_v0 = surface.load_surf_data(input_path).astype(np.float64)
                ccf_v0 = ccf_v0[ccf_v0 != 0]               
                
                if param == "r":
                #if param =="eccentricity":
                    bins = np.linspace(0, 8, 11)
                    frequency,r_bin =np.histogram(ccf_v0, bins=bins)                                
                    sum_data_r.append(frequency)
                    parameters_r.append(ccf_v0)
                elif param == "sigma":
                #elif param == "polarangle":
                    bins = np.linspace(0, 4.5, 11)
                    frequency,sig_bin =np.histogram(ccf_v0, bins=bins)                                
                    sum_data_sig.append(frequency)
                    parameters_sig.append(ccf_v0)

            if model =="34w" and param =="r":#"eccentricity":
                print(hemi, model, param)
                r_array_34w = np.array(sum_data_r)
            elif model =="34w" and param =="sigma": # "polarangle":
                print(hemi, model, param)
                sig_array_34w = np.array(sum_data_sig)
            elif model =="neonates" and param =="r": #"eccentricity":
                print(hemi, model, param)
                r_array_30w = np.array(sum_data_r)
            elif model =="neonates" and param =="sigma": # "polarangle":
                print(hemi, model, param)
                sig_array_30w = np.array(sum_data_sig)

    group1_data = r_array_34w
    group2_data = r_array_30w
    group3_data = sig_array_34w
    group4_data = sig_array_30w
    # Perform bootstrapping by subjects
    num_bootstraps = 10000
    chi_squared_values = []
    chi_squared_values1 = []
    print(hemi)
    for _ in range(num_bootstraps):
        # Resample subjects with replacement
        resampled_indices = np.random.choice(range(group1_data.shape[0]), size=group1_data.shape[0], replace=True)
        resampled_indices1 = np.random.choice(range(group2_data.shape[0]), size=group2_data.shape[0], replace=True)
        resampled_indices2 = np.random.choice(range(group3_data.shape[0]), size=group3_data.shape[0], replace=True)
        resampled_indices3 = np.random.choice(range(group4_data.shape[0]), size=group4_data.shape[0], replace=True)
        resampled_group1 = group1_data[resampled_indices]
        resampled_group2 = group2_data[resampled_indices1]
        resampled_group3 = group3_data[resampled_indices2]
        resampled_group4 = group4_data[resampled_indices3]
        
        # Calculate chi-squared statistic for the resampled data
        chi_squared,_ = chi_squared_test(resampled_group1, resampled_group2)
        chi_squared_values.append(chi_squared)
        # Calculate chi-squared statistic for the resampled data
        chi_squared1,_ = chi_squared_test(resampled_group3, resampled_group4)
        chi_squared_values1.append(chi_squared1)

    # Calculate p-value
    observed_chi_squared, p = chi_squared_test(group1_data, group2_data)
    print(p)
    p_value = np.mean(np.array(chi_squared_values) >= observed_chi_squared)
    print("Observed chi-squared statistic for correlation:", observed_chi_squared)
    print("Bootstrapped p-value:", p_value)

    observed_chi_squared1, p1 = chi_squared_test(group3_data, group3_data)
    print(p1)
    p_value = np.mean(np.array(chi_squared_values1) >= observed_chi_squared1)

    print("Observed chi-squared statistic for sigma:", observed_chi_squared)
    print("Bootstrapped p-value:", p_value)


                

