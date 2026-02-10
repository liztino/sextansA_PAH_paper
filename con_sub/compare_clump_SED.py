#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plots all of the clump SEDs and compares with original

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import seaborn as sns

from get_pivot_wave import get_pivot_wave

# open original SED
reg_name = 'SexA_pah_small_box'
textdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/obs/'
file_name = f'sexa_SED_region_{reg_name}.txt'

reg_sed = ascii.read(textdir + file_name)

# open clump table
# using the sum, not max or median (since it pretty much looks the same)
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_name = 'filter_flux_dendro_clump_mask_corrected.txt'

clump_sed = ascii.read(clump_path + clump_name)

plt.figure(1, figsize = (10,5))
plt.clf()
plt.hlines(reg_sed['med'], reg_sed['mu_low'], reg_sed['mu_high'], lw = 3.0, colors = 'k')
# plt.title('Region ' + str(reg_name))
# plt.xlabel('Wavelength ($\mu$m)', size = 'large')

# get filters then their central wavelengths
filt_list = ['F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M', 'F560W','F770W', 'F1000W', 'F1130W', 'F1500W']

wave = np.zeros(len(filt_list))

for i, filt in enumerate(filt_list):
    wave[i] = get_pivot_wave(filt)
    

norm_filt = 'F1500W'

sed_norm = reg_sed[reg_sed['filt'] == norm_filt]['med']

plot_clumps = [1, 6, 7, 8, 10, 11, 14]
colors = sns.color_palette('hls', len(plot_clumps))

# symbol array 
symb = ['o', '^', 's', 'p', '*', 'x', '8']

c_i = 0
# get value for each clump
for i in range(len(clump_sed)):
    if clump_sed['clump_num'][i] in plot_clumps:

        flux_arr = np.zeros(len(filt_list))
        err_arr = np.zeros(len(filt_list))
        for j, filt in enumerate(filt_list):
            flux_arr[j] = clump_sed[filt + '_flux'][i]
            err_arr[j] = clump_sed[filt + '_err'][i]
                
        # normalize at 7 micron
        norm = sed_norm/clump_sed[norm_filt + '_flux'][i]
    
        plt.scatter(wave, flux_arr*norm, alpha = 0.5, color = colors[c_i], label = clump_sed['clump_num'][i], marker = symb[c_i], s = 100)
        # plt.errorbar(wave, flux_arr*norm, yerr = err_arr*1e6, alpha = 0.5, fmt = 'none', capsize=5, capthick=2, color = colors[c_i])

        c_i = c_i + 1

plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Normalized Flux ($\mu$Jy)')
plt.legend(loc='upper right')       
plt.xlim(0, 20)
        
savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/clump_SEDs/'
savename = 'clump_SED_compare_norm_{:s}.png'.format(norm_filt)
plt.savefig(savedir + savename, dpi = 300)
        
        
        
        
        
        