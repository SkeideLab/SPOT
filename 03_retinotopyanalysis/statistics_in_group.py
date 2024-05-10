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
    dist1 = np.histogram(dist1,bins=range(20))[0]   
    dist2 = np.histogram(dist2,bins=range(20))[0]
    chi2, p, _, _ = chi2_contingency([dist1,dist2])
    return chi2, p

PREFIX_MODEL = (
    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_{hemi}_label"
)

BENSON_ECC = "-eccentricity_desc-real_v2th00"
BENSON_ANGLE = "-eccentricity_desc-real_v2th00"
INPUTSUFFIX_FOVEAL = "_desc-foveal_roi.shape.gii"
INPUTSUFFIX_MIDDLE = "_desc-middle_roi.shape.gii"
INPUTSUFFIX_PERIPHERAL = "_desc-peripheral_roi.shape.gii"
INPUTSUFFIX_DORSAL = "_desc-dorsal_roi.shape.gii"
INPUTSUFFIX_VENTRAL = "_desc-ventral_roi.shape.gii"

subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT.csv')
sub_num = len(subject_info["sub_id"])

            
for param, inpath in [
                (BENSON_ECC, INPUTSUFFIX_FOVEAL),
                [BENSON_ECC, INPUTSUFFIX_MIDDLE],
                (BENSON_ECC, INPUTSUFFIX_PERIPHERAL),
                (BENSON_ANGLE, INPUTSUFFIX_VENTRAL),
                (BENSON_ANGLE, INPUTSUFFIX_DORSAL),
            ]:
    sum_data=[]
    parameters = []   
    for hemi in ["L"]:        
        # FORMAT PATHS FOR INPUT AND OUTPUT                  
        for model in ["real"]:
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
                input_path = prefix_model + param + inpath
                
                ccf_v0 = surface.load_surf_data(input_path).astype(np.float64)               
                
                if inpath == INPUTSUFFIX_FOVEAL:
                    frequency,f_bin =np.histogram(ccf_v0, bins=range(20))                                
                    sum_data.append(frequency)
                    parameters.append(ccf_v0)
                elif inpath == INPUTSUFFIX_MIDDLE:
                    frequency,m_bin =np.histogram(ccf_v0, bins=range(20))                                
                    sum_data.append(frequency)
                    parameters.append(ccf_v0)
                elif inpath == INPUTSUFFIX_PERIPHERAL:
                    frequency,p_bin =np.histogram(ccf_v0, bins=range(20))                                
                    sum_data.append(frequency)
                    parameters.append(ccf_v0)
                elif inpath == INPUTSUFFIX_DORSAL:
                    frequency,d_bin =np.histogram(ccf_v0, bins=range(0,180,10))                                
                    sum_data.append(frequency)
                    parameters.append(ccf_v0)
                elif inpath == INPUTSUFFIX_VENTRAL:
                    frequency,v_bin =np.histogram(ccf_v0, bins=range(0,180,10))                                
                    sum_data.append(frequency)
                    parameters.append(ccf_v0)

        

    if inpath == INPUTSUFFIX_FOVEAL:        
        foveal_array = np.array(sum_data)
        fig, ax = plt.subplots()
        plt.stairs(np.sum(sum_data,0),f_bin, fill=True)
        ax.set_xlabel('Eccentricity')
        ax.set_ylabel('Number of vertex')
        ax.set_title(f"Histogram of {model} data on hemisphere {hemi} for foveal")
        #plt.show()
        foveal_observed_frequencies = np.array(parameters)
    elif inpath == INPUTSUFFIX_MIDDLE:        
        middle_array = np.array(sum_data)
        fig, ax = plt.subplots()
        plt.stairs(np.sum(sum_data,0),m_bin, fill=True)
        ax.set_xlabel('Eccentricity')
        ax.set_ylabel('Number of vertex')
        ax.set_title(f"Histogram of {model} data on hemisphere {hemi} for middle")
        #plt.show()
        middle_observed_frequencies = np.array(parameters)
    elif inpath == INPUTSUFFIX_PERIPHERAL:        
        peripheral_array = np.array(sum_data)
        fig, ax = plt.subplots()
        plt.stairs(np.sum(sum_data,0),p_bin, fill=True)
        ax.set_xlabel('Eccentricity')
        ax.set_ylabel('Number of vertex')
        ax.set_title(f"Histogram of {model} data on hemisphere {hemi} for peripheral")
        #plt.show()
        peripheral_observed_frequencies = np.array(parameters)
    elif inpath == INPUTSUFFIX_DORSAL:        
        dorsal_array = np.array(sum_data)
        fig, ax = plt.subplots()
        plt.stairs(np.sum(sum_data,0),d_bin, fill=True)
        ax.set_xlabel('Polarangle')
        ax.set_ylabel('Number of vertex')
        ax.set_title(f"Histogram of {model} data on hemisphere {hemi} for dorsal")
        #plt.show()
        dorsal_observed_frequencies = np.array(parameters)
    elif inpath == INPUTSUFFIX_VENTRAL:        
        ventral_array = np.array(sum_data)
        ffig, ax = plt.subplots()
        plt.stairs(np.sum(sum_data,0),v_bin, fill=True)
        ax.set_xlabel('Polarangle')
        ax.set_ylabel('Number of vertex')
        ax.set_title(f"Histogram of {model} data on hemisphere {hemi} for ventral")
        #plt.show()
        ventral_observed_frequencies = np.array(parameters)

group1_data = dorsal_observed_frequencies
group2_data = ventral_observed_frequencies
# Perform bootstrapping by subjects
num_bootstraps = 10000
chi_squared_values = []
print(group1_data.shape)
print(group2_data.shape)
for _ in range(num_bootstraps):
    # Resample subjects with replacement
    resampled_indices = np.random.choice(range(193), size=193, replace=True)
    resampled_group1 = group1_data[resampled_indices]
    resampled_group2 = group2_data[resampled_indices]
    
    # Calculate chi-squared statistic for the resampled data
    chi_squared,_ = chi_squared_test(resampled_group1, resampled_group2)
    chi_squared_values.append(chi_squared)

# Calculate p-value
observed_chi_squared, p = chi_squared_test(group1_data, group2_data)
print(p)
p_value = np.mean(np.array(chi_squared_values) >= observed_chi_squared)

print("Observed chi-squared statistic:", observed_chi_squared)
print("Bootstrapped p-value:", p_value)

                

