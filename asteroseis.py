#!/usr/bin/env python
# coding: utf-8


from astropy.io import fits
# from PSD_PATH import * # WHEN THE KIC FILE IS DOWNLOADED, GET RID OF THIS
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import pickle as pickle 
import bokeh 
#from echelle import plot_echelle
#from echelle import interact_echelle
#import csv
from skimage.transform import rescale, resize


# 400
rep_hundreds = "/Users/kanahsmith/Dropbox (Simons Foundation)/INTERNSHIP2022/DATA/DATA_Yu2018/"
path_hundreds = PSD_PATH_IN_REPERTORY(rep_hundreds)


# collecting the data for all 40 fits files
def kepler_data(num_files, path): # where num_files = 40, path = path_data
    stellar = []
    for num in range(num_files):
        ii = num
        hdul = fits.open(path[ii])
        hdu = hdul[0]
        
        KIC = int(path[ii][path[ii].find("_") + 1: path[ii].find("_COR")])
        data = np.array(hdu.data)
        hdul.close()
        freq = data[:,0]*1e6
        power = data[:,1]
        
        stellar.append((freq, power, KIC))
    return stellar

all_kepler = kepler_data(len(path_hundreds), path_hundreds)

CORRECTED_PATH_TO_FILE = "/Users/kanahsmith/Dropbox (Simons Foundation)/INTERNSHIP2022/DATA/"
with open(CORRECTED_PATH_TO_FILE+"eps_dnu_with_espilon_corrected.pkl", "rb") as g: # epsilon is misspelt on purpose
    epsilon, del_nu, KIC, numax = pickle.load(g) 

eps_kic = []
delnu_kic = []
numax_kic = []

for num in range(len(epsilon)):
    eps_kic.append((epsilon[num], KIC[num]))
    delnu_kic.append((del_nu[num], KIC[num]))
    numax_kic.append((numax[num], KIC[num]))

spec_kepler = []
for num in range(len(all_kepler)):
    
    f = all_kepler[num][0]
    p = all_kepler[num][1]
    k = all_kepler[num][2]

    for i in range(len(delnu_kic)):
        if delnu_kic[i][1] == k and eps_kic[i][1] == k and numax_kic[i][1] == k:
            d = delnu_kic[i][0]
            e = eps_kic[i][0]
            n = numax_kic[i][0]
            
            spec_kepler.append([k, f, p, n, d, e]) # freq = x and power = y
# --------------------------------------------------------------

for num in range(len(spec_kepler)):
    kic_num = spec_kepler[num][0]
    nu_max = spec_kepler[num][3]
    delta_nu = spec_kepler[num][4]
    eps = spec_kepler[num][5]
    
    x = spec_kepler[num][1]
    y = spec_kepler[num][2]
    
    slice_interval = 5
    x_lim = (x > nu_max - slice_interval*delta_nu) & (x < nu_max + slice_interval*delta_nu)
    y_max = np.max(spec_kepler[num][2][x_lim])
    
    spec_y = y[(x > nu_max - slice_interval*delta_nu) & (x < nu_max + slice_interval*delta_nu)]
    spec_x = x[(x > nu_max - slice_interval*delta_nu) & (x < nu_max + slice_interval*delta_nu)]
    
    spec_kepler[num].append(spec_x)
    spec_kepler[num].append(spec_y)
# ech_kepler is now [KIC, frequency, PSD, nu max, delta nu, epsilon, freq slice, PSD slice] 
# where both slices are done with 5*delta nu
            

#np.shape(spec_kepler)
# SPEC_KEPLER NEEDS TO BE SAVED TO MY COMPUTER ------ DONE

------------------------------------------------- WANT TO START FROM HERE -------------------------------------------

def create_echelle_diagram(power, freq, dnu, numax, epsilon, min_interval, max_interval):
    """
    Create the echelle diagram corresponding to the given PSD contained in 'power'
    arguments:
        power = PSD associated with a star
        freq = corresponding frequencies
        dnu =  large frequency seperation
        numax = frequency of max aplitude
        epsilon = epsilon parameter
        min_interval = lower boundary (in number of dnu) of the desired frequency interval (default = 3.)
        max_interval = upper boundary (in number of dnu) of the desired frequency interval (default = 3.)
    """
    freq_l0 = int(numax/dnu)*dnu+(epsilon*dnu)     #frequency of the mode ell = 0, the closest to numax (and larger than numax)
    fbin = np.mean(np.diff(freq))     #size of a frequency bin (interval between 2 consecutive frequencies)
    
    #select the part of the PSD with the oscillations
    power_osc = power[np.where((freq>=(freq_l0-min_interval*dnu)) & (freq <= (freq_l0+max_interval*dnu)))]
    freq_osc = freq[np.where((freq>=(freq_l0-min_interval*dnu)) & (freq <= (freq_l0+max_interval*dnu)))]

    #number of frequency bins in dnu
    number_of_bin_for_dnu=int(dnu/fbin)
    #number of slices of the spectra
    sizeh=int(len(freq_osc)/number_of_bin_for_dnu)

    # makes sure that freqs_osc and power_osc have a length that is an integer number of slices
    freq_osc=freq_osc[:sizeh*number_of_bin_for_dnu]
    power_osc=power_osc[:sizeh*number_of_bin_for_dnu]

    # Echelle diagram
    # Stack the slices in a 2D array
    Z = np.reshape(np.transpose(power_osc), (sizeh, number_of_bin_for_dnu))

    return Z


for num in range(len(spec_kepler)):
    fig = plt.figure(figsize = (15,10))
    
    KIC_ech = spec_kepler[num][0]
    delnu_ech = spec_kepler[num][4]
    numax_ech = spec_kepler[num][3]
    freq_ech = spec_kepler[num][6] # spec_x , change to ech_kepler[num][1] for full PSD, ech_kepler[num][6] for sliced by 5
    power_ech = spec_kepler[num][7] # spec_y , change to ech_kepler[num][2] for full PSD, ech_kepler[num][7] for sliced by 5
    epsilon_ech = spec_kepler[num][5]
    
    ech = create_echelle_diagram(power_ech, freq_ech, delnu_ech, numax_ech, epsilon_ech, min_interval = 3.3, max_interval = 2.7)
    rs = resize(ech, (128, 128))
    
    plt.pcolor(rs, cmap='plasma', vmin=np.min(rs), vmax=np.max(rs))
    plt.title(f"KIC {KIC_ech}, f0 = {round((numax_ech/delnu_ech)*delnu_ech+(delnu_ech*epsilon_ech), 3)}", size = 45)
    plt.ylabel('Frequency ($\mu$Hz)', size = 35)
    plt.xlabel('Frequency [$\Delta \\nu$] ($\mu$Hz)', size = 35)

    #path_save = "/Users/kanahsmith/Documents/astero/echelle_plasma/" # THIS PATH NEEDS TO CHANGE TO SUNNYVALE CLUSTER
    img_name = "KIC_{}.png".format(KIC_ech)
    plt.savefig(path_save+img_name, dpi = 800)
    plt.close(fig)
    
