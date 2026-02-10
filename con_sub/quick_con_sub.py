#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 15:20:00 2023

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
plt.ion()

miri_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
filter_path = '/astro/dust_kg/etarantino/JWST/filts/'

# load the middle filter we will be continuum subtracting from
filt_mid = 'F335M'
# mid_name = miri_path + f'miri_{filt_mid}_stage3_asn_fixed_wcs_skysub_i2d.fits'
mid_name = miri_path + f'{filt_mid}_reproject_to_F1500W_rot'
hdu_mid = fits.open(mid_name + '.fits')
# header_mid = hdu_mid['SCI'].header
# data_mid = hdu_mid['SCI'].data
header_mid = hdu_mid[0].header
data_mid = hdu_mid[0].data
# data_mid_err = fits.open(mid_name + '_ERR.fits')['SCI'].data

# load the upper filter     that we will use to match to the other filters
filt_up = 'F360M'
up_name = miri_path + f'{filt_up}_reproject_to_F1500W_rot'
hdu_up = fits.open(up_name  + '.fits')
# header_up = hdu_up['SCI'].header
# data_up = hdu_up['SCI'].data
header_up = hdu_up[0].header
data_up = hdu_up[0].data
s = np.shape(data_up)
# data_up_err = fits.open(up_name + '_ERR.fits')['SCI'].data

# load the lowest wavelength filter
filt_low = 'F300M'
low_name = miri_path + f'{filt_low}_reproject_to_F1500W_rot'
hdu_low = fits.open(low_name + '.fits')
# header_low = hdu_low['SCI'].header
# data_low = hdu_low['SCI'].data
header_low = hdu_low[0].header
data_low = hdu_low[0].data
s = np.shape(data_low)
# data_low_err = fits.open(low_name + '_ERR.fits')['SCI'].data


miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

# gets the median wavelengths for a given filter
def get_wave(filt):
    # open filter file 
    filter_path = '/astro/dust_kg/etarantino/JWST/filts/'
    
    if filt in miri_filts:
        filt_file = 'JWST_MIRI.{:s}.dat'.format(filt)
    else:
        filt_file = 'JWST_NIRCam.{:s}.dat'.format(filt)
    
    # wavelength is in Angstroms, convert to microns
    filt_curve = ascii.read(filter_path + filt_file, names = ['wave', 'trans'])
    filt_curve['wave'] = filt_curve['wave']/1e4
    
    # get position of filter where it is half it's max transmission to plot
    cutoff = np.nanmax(filt_curve['trans'])*0.5
    ind = np.where(filt_curve['trans'] >  cutoff)[0]
    ind1 = ind[0]
    ind2 = ind[-1]
    micron_low = filt_curve['wave'][ind1]
    micron_high = filt_curve['wave'][ind2]
    micron = micron_low + (micron_high-micron_low)/2
    
    return micron

micron_mid = get_wave(filt_mid)
micron_up = get_wave(filt_up)
micron_low = get_wave(filt_low)

# loop through the pixels and fit a simple linear model to the continuum 
s = np.shape(data_mid)
xlen = s[1]
ylen = s[0]

mid_sub = np.zeros(s)
cont = np.zeros(s)

micron_arr = np.array([micron_low, micron_up])

for x in range(xlen-1):
    for y in range(ylen-1):
        # print(y,x)
        # fit the continuum  
        isnan = np.isnan(data_mid[y,x]) or np.isnan(data_low[y,x]) or np.isnan(data_up[y,x])
        if isnan:
            mid_sub[y,x] = np.nan
            cont[y,x] = np.nan
        else:
            flux_arr = np.array([data_low[y,x], data_up[y,x]])
            fit = np.polyfit(micron_arr, flux_arr, 1)
            p = np.poly1d(fit)
            
            cont[y,x] = p(micron_mid)
            
            mid_sub[y,x] = data_mid[y,x] - p(micron_mid)
            
        
# code to plot the data for a check
vmin = 0
vmax = 0.03

plt.figure(1)
plt.clf()
plt.imshow(mid_sub, origin = 'lower', vmin = vmin, vmax = vmax)

plt.figure(2)
plt.clf()
plt.imshow(cont, origin = 'lower', vmin = vmin, vmax = vmax)

plt.figure(3)
plt.clf()
plt.imshow(data_mid, origin = 'lower', vmin = vmin, vmax = vmax)


# save as fits files
save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/interp_method/'

save_name = f'SextansA_{filt_mid}_interp_method_pah.fits'
fits.writeto(save_path + save_name, mid_sub, header = header_up, overwrite = True)

save_name = f'SextansA_{filt_mid}_interp_method_con.fits'
fits.writeto(save_path + save_name, cont, header = header_up, overwrite = True)
