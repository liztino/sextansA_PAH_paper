#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Uses the slope from the continuum subtraction to select out galaxies and stars that are not well subtracted

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np

def load_PAH(filt_mid, k):
    datapath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
    dataname = 'SextansA_{:s}_k_method_pah_k_{:3.2f}.fits'.format(filt_mid, k)
    pah_hdr = fits.open(datapath + dataname)[0]
    pah_data = pah_hdr.data
    pah_head = pah_hdr.header
    
    errname = 'SextansA_{:s}_k_method_err_k_{:3.2f}.fits'.format(filt_mid, k)
    err_hdr = fits.open(datapath + errname)[0]
    err_data = err_hdr.data
    err_head = err_hdr.header
    
    conname = 'SextansA_{:s}_k_method_con_k_{:3.2f}.fits'.format(filt_mid, k)
    con_hdr = fits.open(datapath + conname)[0]
    con_data = con_hdr.data
    con_head = con_hdr.header    
    
    slopename = 'SextansA_{:s}_k_method_slope_k_{:3.2f}.fits'.format(filt_mid, k)
    slope_hdr = fits.open(datapath + slopename)[0]
    slope_data = slope_hdr.data
    slope_head = slope_hdr.header    
    
    return pah_data, pah_head, err_data, con_data, slope_data


pah_11, head_11, err_11, con_11, slope_11 = load_PAH('F1130W', 7.21)

mask = np.zeros_like(pah_11, dtype = bool)

upper_slope = 0.13
lower_slope = -0.01

# Create boolean mask for values to keep
mask = (slope_11 <= upper_slope) & (slope_11 >= lower_slope)

# Apply mask but retain 2D shape
masked_data = np.where(mask, pah_11, np.nan)

plt.figure(1)
plt.clf()
plt.imshow(masked_data, origin = 'lower', vmin = -0.02, vmax = 0.2)