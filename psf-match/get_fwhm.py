#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 14:39:36 2024

@author: etarantino
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from fwhm import fwhm

filt_list = ['F300M', 'F335M', 'F360M', 'F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W']

for filt in filt_list:
    
    psf1_file = f'{filt}_webbpsfv121_rot'
    
    psf_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/psfs/'

    # open psf files
    psf1_hdu = fits.open(psf_path + psf1_file + '.fits')
    psf1 = psf1_hdu[0].data
    psf1_head = psf1_hdu[0].header
    
    width1 = fwhm(psf1, psf1_head)
    
    print('FWHM of {:s}: x = {:5.4f}as y = {:5.4f}as'.format(filt, width1['xarcsec'], width1['yarcsec']))
    