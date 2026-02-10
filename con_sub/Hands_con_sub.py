#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Testing Lindsey Hand's continuum subtraction 

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np

filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'

filt_mid = 'F770W'
filt_low = 'F1000W'
filt_up = 'F1130W'

# load the middle filter we will be continuum subtracting from
mid_name = filepath + f'{filt_mid}_reproject_to_F1500W_rot'
hdu_mid = fits.open(mid_name + '.fits')
header_mid = hdu_mid[0].header
data_mid = hdu_mid[0].data
# pivot_mid = get_pivot_wave(filt_mid)

# load the upper filter     that we will use to match to the other filters
up_name = filepath + f'{filt_up}_reproject_to_F1500W_rot'
hdu_up = fits.open(up_name  + '.fits')
header_up = hdu_up[0].header
data_up = hdu_up[0].data
# pivot_up = get_pivot_wave(filt_up)

# load the lowest wavelength filter
low_name = filepath + f'{filt_low}_reproject_to_F1500W_rot'
hdu_low = fits.open(low_name + '.fits')
header_low = hdu_low[0].header
data_low = hdu_low[0].data
# pivot_low = get_pivot_wave(filt_low)