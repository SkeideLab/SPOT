import pandas as pd
from nilearn import surface
import nibabel as nib
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import chisquare
from scipy.stats import bootstrap   
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

PREFIX_MODEL = (
    "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_{hemi}_label"
)

BENSON_ECC = "-eccentricity_desc-real_v2th00"
BENSON_ANGLE = "-polarangle_desc-real_v2th00"
INPUTSUFFIX_FOVEAL = "_desc-foveal_roi.shape.gii"
INPUTSUFFIX_MIDDLE = "_desc-middle_roi.shape.gii"
INPUTSUFFIX_PERIPHERAL = "_desc-peripheral_roi.shape.gii"
INPUTSUFFIX_DORSAL = "_desc-dorsal_roi.shape.gii"
INPUTSUFFIX_VENTRAL = "_desc-ventral_roi.shape.gii"

subject_info = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal.csv')
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
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.spines['bottom'].set_visible(False)
        ax1.tick_params(axis='x',which='both',bottom=False)
        ax2.spines['top'].set_visible(False)

        ax2.set_ylim(0,0.3)
        ax1.set_ylim(0.5,0.6)
        plot_data = np.sum(sum_data,0)
        normalized_data = plot_data / np.sum(plot_data)
        print(normalized_data)
        bin_centers = (f_bin[:-1] + f_bin[1:]) / 2
        bars1 = ax1.bar(bin_centers, normalized_data, width=(f_bin[1] - f_bin[0]))
        bars2 = ax2.bar(bin_centers, normalized_data, width=(f_bin[1] - f_bin[0]))
        #fig, ax = plt.subplots()
        #plt.stairs(np.sum(sum_data,0),f_bin, fill=True)
        ax2.set_xlabel('Eccentricity')
        f.text(0.04, 0.5, 'Normalized number of vertex', va='center', rotation='vertical')
        ax1.set_title(f"Histogram of {model} data on hemisphere {hemi} for foveal")

        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        plt.show()
        foveal_observed_frequencies = np.array(parameters)
    elif inpath == INPUTSUFFIX_MIDDLE:        
        middle_array = np.array(sum_data)
        fig, ax = plt.subplots()
        plot_data = np.sum(sum_data,0)
        normalized_data = plot_data / np.sum(plot_data)
        print(normalized_data)
        bin_centers = (m_bin[:-1] + m_bin[1:]) / 2
        bars = ax.bar(bin_centers, normalized_data, width=(m_bin[1] - m_bin[0]))
        ax.set_xlabel('Eccentricity')
        fig.text(0.04, 0.5, 'Normalized number of vertex', va='center', rotation='vertical')
        ax.set_title(f"Histogram of {model} data on hemisphere {hemi} for middle")
        plt.show()
        middle_observed_frequencies = np.array(parameters)
    elif inpath == INPUTSUFFIX_PERIPHERAL:        
        peripheral_array = np.array(sum_data)
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.spines['bottom'].set_visible(False)
        ax1.tick_params(axis='x',which='both',bottom=False)
        ax2.spines['top'].set_visible(False)

        ax2.set_ylim(0,0.05)
        ax1.set_ylim(0.8,1)
        plot_data = np.sum(sum_data,0)
        normalized_data = plot_data / np.sum(plot_data)
        print(normalized_data)
        bin_centers = (p_bin[:-1] + p_bin[1:]) / 2
        bars1 = ax1.bar(bin_centers, normalized_data, width=(p_bin[1] - p_bin[0]))
        bars2 = ax2.bar(bin_centers, normalized_data, width=(p_bin[1] - p_bin[0]))
        #fig, ax = plt.subplots()
        #plt.stairs(np.sum(sum_data,0),f_bin, fill=True)
        ax2.set_xlabel('Eccentricity')
        f.text(0.04, 0.5, 'Normalized number of vertex', va='center', rotation='vertical')
        ax1.set_title(f"Histogram of {model} data on hemisphere {hemi} for peripheral")

        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
        plt.show()
        peripheral_observed_frequencies = np.array(parameters)
    elif inpath == INPUTSUFFIX_DORSAL:        
        dorsal_array = np.array(sum_data)
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.spines['bottom'].set_visible(False)
        ax1.tick_params(axis='x',which='both',bottom=False)
        ax2.spines['top'].set_visible(False)

        ax2.set_ylim(0,0.2)
        ax1.set_ylim(0.5,0.7)
        plot_data = np.sum(sum_data,0)
        normalized_data = plot_data / np.sum(plot_data)
        print(normalized_data)
        bin_centers = (d_bin[:-1] + d_bin[1:]) / 2
        bars1 = ax1.bar(bin_centers, normalized_data, width=(d_bin[1] - d_bin[0]))
        bars2 = ax2.bar(bin_centers, normalized_data, width=(d_bin[1] - d_bin[0]))
        #fig, ax = plt.subplots()
        #plt.stairs(np.sum(sum_data,0),f_bin, fill=True)
        ax2.set_xlabel('Polar angle')
        f.text(0.04, 0.5, 'Normalized number of vertex', va='center', rotation='vertical')
        ax1.set_title(f"Histogram of {model} data on hemisphere {hemi} for dorsal")

        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
        plt.show()
        dorsal_observed_frequencies = np.array(parameters)
    elif inpath == INPUTSUFFIX_VENTRAL:        
        ventral_array = np.array(sum_data)
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.spines['bottom'].set_visible(False)
        ax1.tick_params(axis='x',which='both',bottom=False)
        ax2.spines['top'].set_visible(False)

        ax2.set_ylim(0,0.2)
        ax1.set_ylim(0.5,0.7)
        plot_data = np.sum(sum_data,0)
        normalized_data = plot_data / np.sum(plot_data)
        print(normalized_data)
        bin_centers = (v_bin[:-1] + v_bin[1:]) / 2
        bars1 = ax1.bar(bin_centers, normalized_data, width=(v_bin[1] - v_bin[0]))
        bars2 = ax2.bar(bin_centers, normalized_data, width=(v_bin[1] - v_bin[0]))
        #fig, ax = plt.subplots()
        #plt.stairs(np.sum(sum_data,0),f_bin, fill=True)
        ax2.set_xlabel('Polar angle')
        f.text(0.04, 0.5, 'Normalized number of vertex', va='center', rotation='vertical')
        ax1.set_title(f"Histogram of {model} data on hemisphere {hemi} for ventral")

        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
        ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
        plt.show()
        ventral_observed_frequencies = np.array(parameters)
