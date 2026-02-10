#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 09:37:54 2023

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
plt.ion()

miri_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject/'
filter_path = '/astro/dust_kg/etarantino/JWST/filts/'

# load the middle filter we will be continuum subtracting from
filt_mid = 'F1130W'
# mid_name = miri_path + f'miri_{filt_mid}_stage3_asn_fixed_wcs_skysub_i2d.fits'
mid_name = miri_path + f'{filt_mid}_reproject_to_F1500W.fits'
hdu_mid = fits.open(mid_name)
header_mid = hdu_mid['SCI'].header
data_mid = hdu_mid['SCI'].data

# load the upper filter     that we will use to match to the other filters
filt_up = 'F1500W'
up_name = miri_path + f'{filt_up}_reproject_to_F1500W.fits'
hdu_up = fits.open(up_name)
header_up = hdu_up['SCI'].header
data_up = hdu_up['SCI'].data
s = np.shape(data_up)

# load the lowest wavelength filter
filt_low = 'F1000W'
low_name = miri_path + f'{filt_low}_reproject_to_F1500W.fits'
hdu_low = fits.open(low_name)
header_low = hdu_low['SCI'].header
data_low = hdu_low['SCI'].data
s = np.shape(data_low)

x1 = 665
x2 = 736
y1 = 495
y2 = 592

# con_ratio = (data_up/data_low).flatten()
# line_ratio = (data_mid/data_low).flatten()

#F770W
con_ratio = (data_low[y1:y2, x1:x2]/data_up[y1:y2, x1:x2]).flatten()
line_ratio = (data_mid[y1:y2, x1:x2]/data_up[y1:y2, x1:x2]).flatten()

#F1130W
con_ratio = (data_up[y1:y2, x1:x2]/data_low[y1:y2, x1:x2]).flatten()
line_ratio = (data_mid[y1:y2, x1:x2]/data_low[y1:y2, x1:x2]).flatten()


def con_line(x, A, B):
    return A + B*x

plt.figure(1)
plt.clf()
plt.scatter(con_ratio, line_ratio, c = 'k', s = 3)
plt.xlim(0,7)
plt.ylim(0,7)

# x_vals = np.linspace(0, 5, 100)
# lai_vals = con_line(x_vals, 0.35, 0.65)
# plt.plot(x_vals, lai_vals, label = 'L20')

# s23_vals = con_line(x_vals, -0.2, 1.6)
# plt.plot(x_vals, s23_vals, label = 'S23')

# plt.xlabel('F560W/F1000W', size = 'large')
# plt.ylabel('F770W/F1000W', size = 'large')

plt.xlabel('F1500W/F1000W', size = 'large')
plt.ylabel('F1130W/F1000W', size = 'large')

# plt.legend(loc = 'best')

plt.title('PAH Region')

save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_compare/'
save_name = 'F1130W_con_ratio_pah_reg'
plt.savefig(save_path + save_name + '.pdf')