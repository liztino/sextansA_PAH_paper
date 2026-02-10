#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculates uncertainty of a given filter by finding the RMS in an empty region box

@author: etarantino
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian1DKernel, Box1DKernel
import scipy.signal
import numpy as n
import scipy.interpolate
import scipy.ndimage
import glob
from astropy.stats import sigma_clip

plt.ion()

def get_err(filt):


    filepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
    
    # box defined by err_test_box.reg
    x1 = 607
    x2 = 722
    y1 = 426
    y2 = 521
    
    origname =  f'{filt}_reproject_to_F1500W_rot.fits'
    hdu_orig = fits.open(filepath + origname)
    orig = hdu_orig[0].data
    
    # grab the box
    orig_box = orig[y1:y2, x1:x2]
    
    # sigma clip stars from the original data 
    clip_orig_box = sigma_clip(orig_box)
    
    plt.figure()
    plt.clf()
    plt.imshow(clip_orig_box, origin = 'lower')
    plt.title(filt)
    
    orig_rms = np.sqrt(np.nanmean(np.square(clip_orig_box)))
    
    return orig_rms

filt_list = ['F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W' ]
for filt in filt_list:
    print(filt, get_err(filt))