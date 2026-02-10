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
plt.title('Region ' + str(reg_name))
plt.xlabel('Wavelength ($\mu$m)', size = 'large')

# get filters then their central wavelengths
filt_list = ['F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M', 'F560W','F770W', 'F1000W', 'F1130W', 'F1500W']

wave = get_pivot_wave(filt_list)