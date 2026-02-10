#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates three color images of the zoom in area around the PAH region 

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
import numpy.ma as ma
import astropy.visualization.wcsaxes as vis
import astropy.units as u 
import matplotlib
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
import seaborn as sns

from astropy.visualization import make_lupton_rgb

plt.ion()

from get_pivot_wave import get_pivot_wave


def load_filter(filt):
    
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_mcf/'
    
    # load the middle filter we will be continuum subtracting from
    filename = filepath + f'SextansA_jw2391_{filt}_reproject_mcf_north'
    hdu = fits.open(filename + '.fits')
    header = hdu[0].header
    data = hdu[0].data
    pivot = get_pivot_wave(filt)
    
    filt_dict = {'name': filt, 'data': data, 'header': header, 'wave': pivot}
    
    return filt_dict

F300M = load_filter('F300M')
F335M = load_filter('F335M')
F360M = load_filter('F360M')
F560W = load_filter('F560W')
F770W = load_filter('F770W')
F1000W = load_filter('F1000W')
F1130W = load_filter('F1130W')
F1500W = load_filter('F1500W')

filt_list = [F300M, F335M, F360M, F560W, F770W, F1000W, F1130W, F1500W]

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/RGB_images/cutouts/'

box = ''
# region defined by PAH box in the reprojection space
if box == 'big':
    # big PAH box
    x1 = 3864
    x2 = 5362
    y1 = 4740
    y2 = 6147
else:
    x1 = 4164
    x2 = 4973
    y1 = 4908
    y2 = 5671

min_vals = [0,0,0]

### 3.3 RGB ####
image = make_lupton_rgb(F360M['data'], F335M['data'], F300M['data'], stretch=0.2, Q = 5, minimum = min_vals)

w = WCS(F300M['header'])
fig1 = plt.figure(1)
plt.clf()
ax = fig1.add_subplot(projection = w)

im = ax.imshow(image, origin = 'lower')

plt.xlim(x1, x2)
plt.ylim(y1, y2)

ax.set_xticks([])
ax.set_yticks([])

vis.add_scalebar(ax = ax, length = 1 * u.arcsec, corner = 'bottom right', label = '1" = 7 pc', color = 'white')


savename = f'SextansA_lupton_RGB_F360M_F335M_F300M_pah_box_{box}_scalebar.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)


### 7.7 RGB ####
image = make_lupton_rgb(F1000W['data'], F770W['data'], F560W['data'], stretch=0.3, Q = 5, minimum = min_vals)

w = WCS(F560W['header'])
fig2 = plt.figure(2)
plt.clf()
ax = fig2.add_subplot(projection = w)

im = ax.imshow(image, origin = 'lower')

plt.xlim(x1, x2)
plt.ylim(y1, y2)

ax.set_xticks([])
ax.set_yticks([])

vis.add_scalebar(ax = ax, length = 1 * u.arcsec, corner = 'bottom right', label = '1" = 7 pc', color = 'white')


savename = f'SextansA_lupton_RGB_F1000W_F770W_F560W_pah_box_{box}_scalebar.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)

### 11.3 RGB ####
min_vals = [0.01, 0.01, 0.01]
image = make_lupton_rgb(F1500W['data'], F1130W['data'], F1000W['data'], stretch=0.3, Q = 5, minimum = min_vals)

w = WCS(F1000W['header'])
fig3 = plt.figure(3)
plt.clf()
ax = fig3.add_subplot(projection = w)

im = ax.imshow(image, origin = 'lower')

plt.xlim(x1, x2)
plt.ylim(y1, y2)

ax.set_xticks([])
ax.set_yticks([])

vis.add_scalebar(ax = ax, length = 1 * u.arcsec, corner = 'bottom right', label = '1" = 7 pc', color = 'white')


savename = f'SextansA_lupton_RGB_F1500W_F1130W_F1000W_pah_box_{box}_scalebar.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)