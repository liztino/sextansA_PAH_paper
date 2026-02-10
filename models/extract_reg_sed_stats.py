#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:14:29 2023

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from astropy.visualization import make_lupton_rgb
from regions.core import PixCoord
from regions.shapes.rectangle import RectanglePixelRegion
from astropy import wcs
from regions import Regions
import matplotlib as mpl
from astropy.modeling import models, fitting
from astropy import units as u
import seaborn as sns


plt.ion()

# constants
c = 3e14        # in microns/sec
h = 6.62e-27    # in erg s 

# data paths
nircam_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/nircam/'
miri_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/miri/'

# analysis paths
reg_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/reg/'
savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/init/'
textdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/obs/'


reg_name = 'SexA_pah_small_box'

filt_list = ['F1500W', 'F1130W', 'F1000W', 'F770W', 'F560W', 'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M']

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

reg = Regions.read(reg_path + reg_name + '.reg')[0]

micron_low = np.zeros(len(filt_list))
micron_high = np.zeros(len(filt_list))
flux = np.zeros(len(filt_list))
avg = np.zeros(len(filt_list))
med = np.zeros(len(filt_list))
std = np.zeros(len(filt_list))
max_val = np.zeros(len(filt_list))
max_err = np.zeros(len(filt_list))

c_list = []

for i, filt in enumerate(filt_list):
    # check if it's a miri or nircam filter
    if filt in miri_filts:
        infile = glob.glob(miri_path + f'*{filt}*')[0]
            
        filt_file = '{:s}_mean_system_throughput.txt'.format(filt)
        filter_path = '/Users/etarantino/Documents/JWST/filts/miri/'
        spitzer = False


    elif filt in nircam_filts:
        filt_name = filt.lower()
        infile = nircam_path + f'jw02391007001_nircam_{filt_name}_asn_modb_1fcor_bkgrsub_i2d.fits'
        filt_file = '{:s}_mean_system_throughput.txt'.format(filt)
        filter_path = '/Users/etarantino/Documents/JWST/filts/nircam/'
        spitzer = False
    
    
    # # original data path not PSF matched
    # filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
    # filename = f'{filt}_reproject_to_F1500W_rot.fits'
    
    # # load error file
    # err_file = 'err_reg_reproject_to_F1500W_rot.txt'
    # err_filt = ascii.read(filepath + err_file)
    
    # # load the filter and info
    # hdu = fits.open(filepath + filename)
    
    # # load the filter and info
    # hdu = fits.open(filepath + filename)
    # header = hdu[0].header
    # data = hdu[0].data

    # open fits file
    hdu = fits.open(infile)

    header = hdu['SCI'].header
    data = hdu['SCI'].data
    err = hdu['ERR'].data
    
    # open filter file 
    # wavelength is in Angstroms, convert to microns
    filt_curve = ascii.read(filter_path + filt_file, names = ['wave', 'trans'])
    filt_curve['wave'] = filt_curve['wave']
    
    # # get the middle of the filter for plotting purposes
    # micron[i] = ((filt_curve['wave'][-1] - filt_curve['wave'][0])/2) + filt_curve['wave'][0]
    # print(micron[i])
    
    # get position of filter where it is half it's max transmission to plot
    cutoff = np.nanmax(filt_curve['trans'])*0.5
    ind = np.where(filt_curve['trans'] >  cutoff)[0]
    ind1 = ind[0]
    ind2 = ind[-1]
    micron_low[i] = filt_curve['wave'][ind1]
    micron_high[i] = filt_curve['wave'][ind2]

    # convert from pixel/degree to pixel/sr
    w = wcs.WCS(header)
    pix_area = wcs.utils.proj_plane_pixel_area(w) * 3600**2 / (4.25e10)

    #converts region to frame of data and creates a mask 
    pixel_reg = reg.to_pixel(w)
    pix_mask = pixel_reg.to_mask()
    shape = np.shape(data)
    mask = pix_mask.to_image(shape)
    final_mask = np.array(mask, dtype=bool)
    
    # add up the flux and multiply by 
    flux[i] = np.nansum(data[final_mask]) * pix_area
    avg[i] = np.nanmean(data[final_mask])
    med[i] = np.nanmedian(data[final_mask])
    std[i] = np.nanstd(data[final_mask])
    max_val[i] = np.nanmax(data[final_mask])
    max_err[i] = err[final_mask][np.argmax(data[final_mask])]
    
    print(flux[i], avg[i], std[i])
    
    # plt.figure()
    # plt.clf()
    # plt.hist(data[final_mask], bins = 50, histtype = 'step', lw = 2, color = 'k')
    # plt.xlabel('Flux (MJy/sr)')
    # plt.ylabel('N Pixels')
    # plt.title('Filter {:s}'.format(filt))    
        
    # plt.axvline(med[i], c = 'r')
    # plt.axvline(avg[i], c = 'g')
    
    # # test to make sure region is correct 
    # test_data = np.zeros(np.shape(data))
    # test_data[final_mask] = data[final_mask]
    # plt.figure()
    # plt.imshow(test_data, origin = 'lower')
    # plt.show()    
        
            
    if filt == 'F335M':
        c_list.append('magenta')
    elif filt == 'F770W':
        c_list.append('cyan')
    elif filt == 'F1130W':
        c_list.append('orange')
    else:
        c_list.append('black')
        

plt.figure(1, figsize = (10,5))
plt.clf()
plt.hlines(med, micron_low, micron_high, lw = 3.0, colors = c_list, zorder = 10000)
plt.title('Region ' + str(reg_name))
plt.xlabel('Wavelength ($\mu$m)', size = 'large')
plt.ylabel('Flux Avg (MJy/sr)', size = 'large')
plt.minorticks_on()
plt.xlim(0.25,20)
plt.ylim(0,4.0)
plt.show()

mic_mid = micron_low + (micron_high - micron_low)/2

savename = f'sexa_SED_region_{reg_name}.pdf'
# plt.savefig(savedir + savename, bbox_inches='tight',pad_inches = 0.1)

# get position of filter where it is half it's max transmission to plot
cutoff = np.nanmax(filt_curve['trans'])*0.5
ind = np.where(filt_curve['trans'] >  cutoff)[0]
ind1 = ind[0]
ind2 = ind[-1]
micron_low[i] = filt_curve['wave'][ind1]
micron_high[i] = filt_curve['wave'][ind2]

# save text file with all information
vals = [filt_list, micron_low, micron_high, mic_mid, flux, avg, med, std]
names = ['filt', 'mu_low', 'mu_high', 'mu_mid', 'tot_flux', 'avg', 'med', 'std', 'max', 'err']

file_name = f'sexa_SED_region_{reg_name}_max_val.txt'

# ascii.write(vals, textdir + file_name, names = names, delimiter = '\t', overwrite = True)

    
    
    
        
