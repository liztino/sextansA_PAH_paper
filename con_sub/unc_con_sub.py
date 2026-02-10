#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Continuum subtraction of filters with ODR

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
from scipy.odr import Model, RealData, ODR
plt.ion()

miri_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject/'
filter_path = '/astro/dust_kg/etarantino/JWST/filts/'

# load the middle filter we will be continuum subtracting from
filt_mid = 'F1130W'
# mid_name = miri_path + f'miri_{filt_mid}_stage3_asn_fixed_wcs_skysub_i2d.fits'
mid_name = miri_path + f'{filt_mid}_reproject_to_F1500W'
hdu_mid = fits.open(mid_name + '.fits')
header_mid = hdu_mid['SCI'].header
data_mid = hdu_mid['SCI'].data
data_mid_err = fits.open(mid_name + '_ERR.fits')['SCI'].data

# load the upper filter     that we will use to match to the other filters
filt_up = 'F1500W'
up_name = miri_path + f'{filt_up}_reproject_to_F1500W'
hdu_up = fits.open(up_name  + '.fits')
header_up = hdu_up['SCI'].header
data_up = hdu_up['SCI'].data
s = np.shape(data_up)
data_up_err = fits.open(up_name + '_ERR.fits')['SCI'].data

# load the lowest wavelength filter
filt_low = 'F1000W'
low_name = miri_path + f'{filt_low}_reproject_to_F1500W'
hdu_low = fits.open(low_name + '.fits')
header_low = hdu_low['SCI'].header
data_low = hdu_low['SCI'].data
s = np.shape(data_low)
data_low_err = fits.open(low_name + '_ERR.fits')['SCI'].data


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
mid_sub_err = np.zeros(s)

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
            err_arr = np.array([data_low_err[y,x], data_up_err[y,x]])
            
            # fit, cov = np.polyfit(micron_arr, flux_arr, 1,  w=1/err_arr, cov='unscaled')
            fit= np.polyfit(micron_arr, flux_arr, 1)

            p = np.poly1d(fit)
        
            cont[y,x] = p(micron_mid)
            
            mid_sub[y,x] = data_mid[y,x] - p(micron_mid)
            
            mid_sub_err[y,x] = np.sqrt(data_low_err[y,x]**2 + data_up_err[y,x]**2 + data_mid_err[y,x]**2)
            
        
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
save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images/'

save_name = f'SextansA_{filt_mid}_sub_match_pah_v2.fits'
fits.writeto(save_path + save_name, mid_sub, header = header_up, overwrite = True)

save_name = f'SextansA_{filt_mid}_sub_match_con_v2.fits'
fits.writeto(save_path + save_name, cont, header = header_up, overwrite = True)

save_name = f'SextansA_{filt_mid}_sub_match_err_v2.fits'
fits.writeto(save_path + save_name, mid_sub_err, header = header_up, overwrite = True)

save_name = f'SextansA_{filt_mid}_sub_match_snr_v2.fits'
mid_sub_snr = mid_sub/mid_sub_err
fits.writeto(save_path + save_name, mid_sub_snr, header = header_up, overwrite = True)