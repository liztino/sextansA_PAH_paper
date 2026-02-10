#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Tests the uncertainty given from pipeline and PSF matched/regridded

@author: etarantino
"""

from astropy.io import fits, ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian1DKernel, Box1DKernel
import scipy.signal
import numpy as n
import scipy.interpolate
import scipy.ndimage
import glob
from astropy.stats import sigma_clip
from astropy.table import Table

plt.ion()

filepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
filt_list = ['F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W' ]
# filt_list = ['F300M']

# box defined by err_test_box.reg
x1 = 607
x2 = 722
y1 = 426
y2 = 521

rms_arr = np.zeros(len(filt_list))
err_arr = np.zeros(len(filt_list))

for i,filt in enumerate(filt_list):
    # errname =  f'{filt}_reproject_to_F1500W_rot_ERR_kernelsq_norm.fits'
    # hdu_err = fits.open(filepath + errname)
    # err = hdu_err[0].data
    
    origname =  f'{filt}_reproject_to_F1500W_rot.fits'
    hdu_orig = fits.open(filepath + origname)
    orig = hdu_orig[0].data
    
    # grab the box
    # err_box = err[y1:y2, x1:x2]
    orig_box = orig[y1:y2, x1:x2]
    
    # sigma clip stars from the original data 
    clip_orig_box = sigma_clip(orig_box)
    
    # plt.figure()
    # plt.clf()
    # plt.imshow(clip_orig_box, origin = 'lower')
    # plt.title(filt)
    
    orig_rms = np.sqrt(np.nanmean(np.square(clip_orig_box)))
    
    print(filt)
    print('RMS of data', orig_rms)
    
    # med_err = np.nanmedian(err_box)
    # print('Median of ERR', med_err)
    
    rms_arr[i] = orig_rms
    # err_arr[i] = med_err
    
savename = 'err_reg_reproject_to_F1500W_rot.txt'

tab = Table()
tab['filt'] = filt_list
tab['err'] = rms_arr

ascii.write(tab, filepath + savename, overwrite = True)

    
# plt.figure(1)
# plt.clf()
# plt.scatter(rms_arr, err_arr)
# # plt.plot(rms_arr, rms_arr)
# for i in range(len(rms_arr)):
#     plt.text(rms_arr[i], err_arr[i], s = filt_list[i])
    
# lim1 = 0
# lim2 = 0.025
# plt.xlim(lim1, lim2)
# plt.ylim(lim1, lim2)

# # xx = np.linspace(np.min(rms_arr), np.max(rms_arr))
# xx = np.linspace(0, 0.025)
# plt.plot(xx, xx, c = 'k')

# plt.xlabel('RMS of data')
# plt.ylabel('Median of ERR array')
# plt.title('PSF Matched and Regridded Images')

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/plots/'
# savename = 'RMS_Median_err_compare_kernelsq_norm.pdf'

# plt.savefig(savepath + savename)

