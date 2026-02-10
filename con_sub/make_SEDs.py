#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates SEDs of points with a full PAH detection in each band

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
plt.ion()

con_sub_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images/'

F335M_file = 'SextansA_F335M_sub_match_pah.fits'
F335M_snr_file = 'SextansA_F335M_sub_match_snr.fits'
F335M_err_file = 'SextansA_F335M_sub_match_err.fits'

F770W_file = 'SextansA_F770W_sub_match_pah_v2.fits'
F770W_snr_file = 'SextansA_F770W_sub_match_snr_v2.fits'
F770W_err_file = 'SextansA_F770W_sub_match_err_v2.fits'

F1130W_file = 'SextansA_F1130W_sub_match_pah_v2.fits'
F1130W_snr_file = 'SextansA_F1130W_sub_match_snr_v2.fits'
F1130W_err_file = 'SextansA_F1130W_sub_match_err_v2.fits'

F335M_hdu = fits.open(con_sub_dir + F335M_file)
F335M_data = F335M_hdu[0].data
F335M_head = F335M_hdu[0].header
F335M_snr = fits.open(con_sub_dir + F335M_snr_file)[0].data
F335M_err = fits.open(con_sub_dir + F335M_err_file)[0].data


F770W_hdu = fits.open(con_sub_dir + F770W_file)
F770W_data = F770W_hdu[0].data
F770W_head = F770W_hdu[0].header
F770W_snr = fits.open(con_sub_dir + F770W_snr_file)[0].data
F770W_err = fits.open(con_sub_dir + F770W_err_file)[0].data

F1130W_hdu = fits.open(con_sub_dir + F1130W_file)
F1130W_data = F1130W_hdu[0].data
F1130W_head = F1130W_hdu[0].header
F1130W_snr = fits.open(con_sub_dir + F1130W_snr_file)[0].data
F1130W_err = fits.open(con_sub_dir + F1130W_err_file)[0].data

# grab the area that has PAHs
x1 = 665
x2 = 746
y1 = 495
y2 = 625

# grab the area that has PAHs
# x1 = 665
# x2 = 690
# y1 = 510
# y2 = 580

reg_ind = np.zeros(np.shape(F1130W_data), dtype = bool)
F335M_ind = np.zeros(np.shape(F1130W_data), dtype = bool)
F770W_ind = np.zeros(np.shape(F1130W_data), dtype = bool)
F1130W_ind = np.zeros(np.shape(F1130W_data), dtype = bool)

reg_ind[y1:y2, x1:x2] = True

F335M_ind[F335M_snr > 3] = True
F770W_ind[F770W_snr > 3] = True
F1130W_ind[F1130W_snr > 3] = True

master_ind = np.logical_and(np.logical_and(np.logical_and(reg_ind, F335M_ind), F770W_ind), F1130W_ind)

data_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject/'
filter_path = '/astro/dust_kg/etarantino/JWST/filts/'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/SEDs/pixel_SED/'

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

F300M = {} 
F335M = {} 
F360M = {} 
F560W = {} 
F770W = {} 
F1000W = {}
F1130W = {} 
F1500W = {}
filt_dict = [F300M, F335M, F360M,F560W,  F770W, F1000W, F1130W, F1500W ]
filt_name = ['F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W' ]
color_array =['k', 'm', 'k', 'k', 'c', 'k', 'orange', 'k']

# populate the filter dictionary
for i, filt in enumerate(filt_dict):
    # put the name of the filter
    filt['name'] = filt_name[i]
    
    # load the reprojected file
    file_name = data_path + '{:s}_reproject_to_F1500W'.format(filt['name'])
    hdu = fits.open(file_name + '.fits')
    header = hdu['SCI'].header
    data = hdu['SCI'].data
    
    err = fits.open(file_name + '_ERR.fits')['SCI'].data

    #  asign the data and header
    filt['data'] = data
    filt['header'] = header
    filt['err'] = err
    
    # find the transmission of the filter
    if filt['name'] in nircam_filts:
        filt_file = 'JWST_NIRCam.{:s}.dat'.format(filt['name'])
    else:
        filt_file = 'JWST_MIRI.{:s}.dat'.format(filt['name'])
        
    # open filter file 
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
    micron_med = micron_low + (micron_high - micron_low )/2

    filt['mu_low'] = micron_low
    filt['mu_high'] = micron_high
    filt['mu_med'] = micron_med
    
# now make the SEDs
plt.figure(1, figsize = (10,5))
plt.clf()

flux = np.zeros(len(filt_dict))
err = np.zeros(len(filt_dict))
mu_low = np.zeros(len(filt_dict))
mu_high = np.zeros(len(filt_dict))
name = []


points  = 0 
for x in np.arange(x1, x2, 1):
    for y in np.arange(y1, y2, 1):
        if master_ind[y,x]:
            if F770W['data'][y,x] > 2.0:
                points = points + 1
                
                for i, filt in enumerate(filt_dict):
                    flux[i] = filt['data'][y,x]
                    err[i] = filt['err'][y,x]
                    mu_low[i] = filt['mu_low']
                    mu_high[i] = filt['mu_high']
                    name.append(filt['name'])
                    
                plt.clf()
                plt.hlines(flux, mu_low, mu_high, lw = 3.0, colors = color_array, zorder = 10000)
                # for i in range(len(filt_dict)):
                #     print('flux', flux[i])
                #     print('ymin', flux[i] - err[i])
                #     print('ymax',flux[i] + err[i])
                #     plt.axvspan(xmin = mu_low[i], xmax = mu_high[i], ymin = flux[i] - err[i], ymax = flux[i] + err[i], color = color_array[i], alpha = 0.3)
                
                plt.xlabel('Wavelength ($\mu$m)', size = 'large')
                plt.ylabel('Flux (MJy/sr)', size = 'large')
                plt.title('Position x={:d} y={:d}'.format(x,y))
                plt.minorticks_on()
                
                savename = 'SextansA_SED_SNR_mask_x_{:d}_y{:d}'.format(x,y)
                plt.savefig(savepath + '/bright_pixels/' + savename + '.png', dpi = 300)
                
                savedata = np.array([name, flux, err])
                ascii.write(savepath + '/bright_pixels/' + savename + '.txt')
print(points)

# make integrated SED
flux_sum =  np.zeros(len(filt_dict))
for i, filt in enumerate(filt_dict):
    flux_sum[i] = np.nansum(filt['data'][master_ind])
    mu_low[i] = filt['mu_low']
    mu_high[i] = filt['mu_high']

plt.clf()
plt.hlines(flux_sum/len(filt['data'][master_ind]), mu_low, mu_high, lw = 3.0, colors = color_array, zorder = 10000)

plt.xlabel('Wavelength ($\mu$m)', size = 'large')
plt.ylabel('Flux (MJy/sr)', size = 'large')
plt.minorticks_on()


