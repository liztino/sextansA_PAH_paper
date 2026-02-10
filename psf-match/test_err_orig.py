#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Tests the uncertainty given from pipeline for the original images

@author: etarantino
"""

from astropy.io import fits, ascii
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian1DKernel, Box1DKernel
import scipy.signal
import numpy as n
import scipy.interpolate
import scipy.ndimage
import glob
from astropy.stats import sigma_clip
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5 
import astropy.units as u
from astropy.table import Table

plt.ion()

filepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'
filt_list = ['F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W']
# filt_list = ['F300M']

corner1 = SkyCoord(['10:11:08.5969 -04:42:47.085', '10:11:07.7728 -04:42:36.525'], frame = FK5, unit = (u.hourangle, u.deg))
# corner1 = SkyCoord(['10:11:07.7728 -04:42:36.525'], frame = FK5, unit = (u.hourangle, u.deg))
rms_arr = np.zeros(len(filt_list))
err_arr = np.zeros(len(filt_list))


miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
lw_filts = {'F300M', 'F335M', 'F360M'}



for i,filt in enumerate(filt_list):
    # errname =  f'jw2391_{filt}_err_align_north.fits'
    # hdu_err = fits.open(filepath + errname)
    # err = hdu_err[0].data
    
    # origname =  f'jw2391_{filt}_align_north.fits'
    # hdu_orig = fits.open(filepath + origname)
    # orig = hdu_orig[0].data
    # head = hdu_orig[0].header
    
    # original data path not PSF matched
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/'
    
    if filt in miri_filts:
        filepath = filepath + 'miri/'
        
        filename = glob.glob(filepath + f'*{filt}*i2d.fits')[0]
        
    elif filt in nircam_filts:
        filepath = filepath + 'nircam/'
        
        if filt in lw_filts:
            filename = glob.glob(filepath + '*{:s}*jhat_i2d.fits'.format(filt.lower()))[0]
        else:
            filename = glob.glob(filepath + '*{:s}*i2d.fits'.format(filt.lower()))[0]
    
    # load the filter and info
    hdu = fits.open(filename)['SCI']
    orig = hdu.data
    head = hdu.header
    err = fits.open(filename)['ERR'].data
    
    w = WCS(head)
    
    xy = w.world_to_pixel(corner1)
    if xy[0][0] < xy[1][0]:
        x1 = int(xy[0][0])
        y1 = int(xy[0][1])
        x2 = int(xy[1][0])
        y2 = int(xy[1][1])
    else:
        x2 = int(xy[0][0])
        y2 = int(xy[0][1])
        x1 = int(xy[1][0])
        y1 = int(xy[1][1])
    
    # # grab the box
    err_box = err[y1:y2, x1:x2]
    orig_box = orig[y1:y2, x1:x2]
    
    # sigma clip stars from the original data 
    clip_orig_box = sigma_clip(orig_box)
    
    orig_rms = np.sqrt(np.nanmean(np.square(clip_orig_box)))
    
    print(filt)
    print('RMS of data', orig_rms)
    
    med_err = np.nanmedian(err_box)
    print('Median of ERR', med_err)
    
    rms_arr[i] = orig_rms
    err_arr[i] = med_err
    
    
    plt.figure()
    plt.clf()
    plt.imshow(clip_orig_box, origin = 'lower')
    plt.title(filt)
    ax = plt.gca()
    plt.text(0.95, 0.95, 'RMS of data = {:3.5f} \nMed of err = {:3.5f}'.format(orig_rms, med_err), 
             transform=ax.transAxes, va = 'top', ha = 'right')
    
    
plt.figure(1)
plt.clf()
plt.scatter(rms_arr, err_arr)
for i in range(len(rms_arr)):
    plt.text(rms_arr[i], err_arr[i], s = filt_list[i])
    
plt.xlabel('RMS of data')
plt.ylabel('Median of ERR array')
plt.title('Original Images')

lim1 = 0.001
lim2 = 0.035
plt.xlim(lim1, lim2)
plt.ylim(lim1, lim2)

xx = np.linspace(lim1, lim2)
plt.plot(xx, xx, c = 'k')

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/plots/'
# savename = 'Orig_image_RMS_Median_err_compare_aug24_data_orig_noalign.pdf'
# plt.savefig(savepath + savename)

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/err_vals/'
# savename = 'SextansA_err_vals_aug24_data_orig_noalign.txt'

# tab = Table()
# tab['filter'] = filt_list
# tab['rms_data'] = rms_arr
# tab['med_err'] = err_arr

# ascii.write(tab, savepath + savename)



