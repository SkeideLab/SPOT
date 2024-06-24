import pandas as pd
from nilearn import surface
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# Get the screen dimensions using tkinter
import tkinter as tk
root = tk.Tk()
screen_width_px = root.winfo_screenwidth()
screen_height_px = root.winfo_screenheight()
root.destroy()

# Convert pixels to inches (assuming 100 dpi for conversion)
dpi = 100
screen_width_in = screen_width_px / dpi
screen_height_in = screen_height_px / dpi

# Set desired figure size to half of screen width and full screen height
fig_width_in = screen_width_in / 2
fig_height_in = screen_height_in

PREFIX_MODEL = (
    "/data/p_02915/dhcp_derivatives_SPOT/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_"
)
PREFIX_MODEL_fetal = (
    "/data/p_02915/dhcp_derivatives_SPOT/fetal/ccfmodel/sub-{sub}/ses-{ses}/"
    "sub-{sub}_ses-{ses}_hemi-{hemi}_mesh-fsaverage_dens-164k_"
)

#PARAM = "desc-real_{param}.gii"
PARAM = "label-{param}_desc-real_roi-v2th00_metric.gii"

subject_fetal = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_3.csv')
subject_fetal_young = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_fetal_4.csv')
subject_l37 = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_less_37.csv')
subject_o37 = pd.read_csv('/data/p_02915/SPOT/dhcp_subj_path_SPOT_over_37.csv')
sub_num_fetal = len(subject_fetal["sub_id"])
sub_num_fetal_young = len(subject_fetal_young["sub_id"])
sub_num_l37 = len(subject_l37["sub_id"])
sub_num_o37 = len(subject_o37["sub_id"])
 
for hemi in ["L", "R"]:              
    for model in ["l37", "o37", "fetal", "fetal_young"]:
        sum_data_ecc=[]
        parameters_ecc = []   
        sum_data_png=[]
        parameters_png = []           
        if model =="l37":
            subject_info = subject_l37
            sub_num = sub_num_l37   
        elif model == "o37":
            subject_info = subject_o37
            sub_num = sub_num_o37
        elif model == "fetal":
            subject_info = subject_fetal
            sub_num = sub_num_fetal
        elif model == "fetal_young":
            subject_info = subject_fetal_young
            sub_num = sub_num_fetal_young
        #for param in ["r", "sigma"]:                   
        for param in ["eccentricity", "polarangle"]:
        # FORMAT PATHS FOR INPUT AND OUTPUT 
            for index, row in subject_info.iterrows(): 
                sub_id = subject_info.at[index, "sub_id"]
                sess_id = subject_info.at[index, "sess_id"]
                sub = sub_id.replace('sub-','')
                ses = sess_id.replace('ses-','')
                if model =="l37":
                    prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                elif model == "o37":
                    prefix_model = PREFIX_MODEL.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                elif model == "fetal":
                    prefix_model = PREFIX_MODEL_fetal.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                elif model == "fetal_young":
                    prefix_model = PREFIX_MODEL_fetal.format(
                    sub=sub,
                    ses=ses,
                    hemi=hemi,
                )
                
                pram_model = PARAM.format(
                    param=param,
                )
                input_path = prefix_model + pram_model
                
                ccf_v0 = surface.load_surf_data(input_path).astype(np.float64)
                #ccf_v0 = ccf_v0[ccf_v0 != 0]               
                
                #if param == "r":
                if param =="eccentricity":
                    e_bins = np.linspace(0, 23, 20)
                    frequency,r_bin =np.histogram(ccf_v0[ccf_v0 != 0], bins=e_bins)                                
                    sum_data_ecc.append(frequency)
                    parameters_ecc.append(ccf_v0)
                #elif param == "sigma":
                elif param == "polarangle":
                    p_bins = np.linspace(0, 180, 30)
                    frequency,sig_bin =np.histogram(ccf_v0[ccf_v0 != 0], bins=p_bins)                                
                    sum_data_png.append(frequency)
                    parameters_png.append(ccf_v0)

            if model =="l37" and param =="eccentricity": 
                print(hemi, model, param)
                ecc_array_l37 = np.array(sum_data_ecc)
                ecc_data_l37=np.array(parameters_ecc).reshape(-1)
                ecc_data_l37=ecc_data_l37[ecc_data_l37 !=0]
            elif model =="l37" and param =="polarangle": 
                print(hemi, model, param)
                png_array_l37 = np.array(sum_data_png)
                png_data_l37=np.array(parameters_png).reshape(-1)
                png_data_l37=png_data_l37[png_data_l37 !=0]
            elif model =="o37" and param =="eccentricity": #"r": #
                print(hemi, model, param)
                ecc_array_o37 = np.array(sum_data_ecc)
                ecc_data_o37=np.array(parameters_ecc).reshape(-1)
                ecc_data_o37=ecc_data_o37[ecc_data_o37 !=0]
            elif model =="o37" and param == "polarangle": #"sigma": #
                print(hemi, model, param)
                png_array_o37 = np.array(sum_data_png)
                png_data_o37=np.array(parameters_png).reshape(-1)
                png_data_o37=png_data_o37[png_data_o37 !=0]
            elif model =="fetal" and param =="eccentricity": #"r": #
                print(hemi, model, param)
                ecc_array_fetal = np.array(sum_data_ecc)
                ecc_data_fetal=np.array(parameters_ecc).reshape(-1)
                ecc_data_fetal=ecc_data_fetal[ecc_data_fetal !=0]
            elif model =="fetal" and param == "polarangle": #"sigma": #
                print(hemi, model, param) 
                png_array_fetal = np.array(sum_data_png)
                png_data_fetal=np.array(parameters_png).reshape(-1)
                png_data_fetal=png_data_fetal[png_data_fetal !=0]
            elif model =="fetal_young" and param =="eccentricity": #"r": #
                print(hemi, model, param)
                ecc_array_fetaly = np.array(sum_data_ecc)
                ecc_data_fetaly=np.array(parameters_ecc).reshape(-1)
                ecc_data_fetaly=ecc_data_fetaly[ecc_data_fetaly !=0]
            elif model =="fetal_young" and param == "polarangle": #"sigma": #
                print(hemi, model, param) 
                png_array_fetaly = np.array(sum_data_png)
                png_data_fetaly=np.array(parameters_png).reshape(-1)
                png_data_fetaly=png_data_fetaly[png_data_fetaly !=0]

    
    #group1_data = np.sum(ecc_array_l37,0)/np.sum(np.sum(ecc_array_l37,0))    
    #group2_data = np.sum(ecc_array_o37,0)/np.sum(np.sum(ecc_array_o37,0))
    #group3_data = np.sum(ecc_array_fetal,0)/np.sum(np.sum(ecc_array_fetal,0))
    #group4_data = np.sum(png_array_l37,0)/np.sum(np.sum(png_array_l37,0))
    #group5_data = np.sum(png_array_o37,0)/np.sum(np.sum(png_array_o37,0))
    #group6_data = np.sum(png_array_fetal,0)/np.sum(np.sum(png_array_fetal,0))
    # Flatten the arrays to create a 1D array for each group
    group1_data = ecc_data_fetal.flatten()
    group2_data = ecc_data_l37.flatten()
    group3_data = ecc_data_o37.flatten()
    group4_data = ecc_data_fetaly.flatten()
    group5_data = png_data_fetal.flatten()
    group6_data = png_data_l37.flatten()
    group7_data = png_data_o37.flatten()
    group8_data = png_data_fetaly.flatten()

    # Prepare data for plotting
    data_l37 = pd.DataFrame({'Eccentricity': group2_data, 'Group': 'Neonates < 37w'})
    data_o37 = pd.DataFrame({'Eccentricity': group3_data, 'Group': 'Neonates > 37w'})
    data_fetal = pd.DataFrame({'Eccentricity': group1_data, 'Group': 'Fetal > 29w'})
    data_fetaly = pd.DataFrame({'Eccentricity': group4_data, 'Group': 'Fetal < 29w'})

    # Combine all data into a single DataFrame
    df = pd.concat([data_fetaly, data_fetal, data_l37, data_o37], ignore_index=True)

     # Define colors for each group
    colors = {'Neonates < 37w': 'limegreen', 'Neonates > 37w': 'mediumaquamarine', 'Fetal < 29w': 'gold', 'Fetal > 29w': 'yellowgreen'}

    # Create subplots
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(fig_width_in, fig_height_in))

    # Plot normalized histograms for each group
    for ax, group in zip(axes, df['Group'].unique()):
        sns.histplot(data=df[df['Group'] == group], x='Eccentricity', ax=ax, color=colors[group], bins=e_bins, stat='density')
        ax.text(-0.3, 0.5, f'{group}', transform=ax.transAxes, 
            fontsize=10, verticalalignment='center', rotation='horizontal')
        ax.set_ylabel('')
        ax.set_xlim(0, 23)
        ax.set_ylim(0, 0.18)

    # Set common x-axis label
    axes[-1].set_xlabel('Eccentricity')

    fig.suptitle(f'Eccentricity on {hemi} hemisphere', fontsize=16)
    
    # Adjust layout and spacing between subplots
    plt.tight_layout(h_pad=0.5)

    # Show plot
    plt.show()

    # Prepare data for plotting
    data_l37 = pd.DataFrame({'Polarangle': group6_data, 'Group': 'Neonates < 37w'})
    data_o37 = pd.DataFrame({'Polarangle': group7_data, 'Group': 'Neonates > 37w'})
    data_fetal = pd.DataFrame({'Polarangle': group5_data, 'Group': 'Fetal > 29w'})
    data_fetaly = pd.DataFrame({'Polarangle': group5_data, 'Group': 'Fetal < 29w'})

    # Combine all data into a single DataFrame
    df = pd.concat([data_fetaly, data_fetal, data_l37, data_o37], ignore_index=True)

    # Create subplots
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(fig_width_in, fig_height_in))

    # Plot normalized histograms for each group
    for ax, group in zip(axes, df['Group'].unique()):
        sns.histplot(data=df[df['Group'] == group], x='Polarangle', ax=ax, color=colors[group], bins=p_bins, stat='density')
        ax.text(-0.3, 0.5, f'{group}', transform=ax.transAxes, 
            fontsize=10, verticalalignment='center', rotation='horizontal')
        ax.set_ylabel('')
        ax.set_xlim(0, 180)
        ax.set_ylim(0, 0.03)

    # Set common x-axis label
    axes[-1].set_xlabel('Polarangle')

    fig.suptitle(f'Polarangle on {hemi} hemisphere', fontsize=16)

    # Adjust layout and spacing between subplots
    plt.tight_layout(h_pad=0.5)

    # Show plot
    plt.show()