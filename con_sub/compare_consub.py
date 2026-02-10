#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates plots that compare the PAH continuum subtraction methods

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
plt.ion()

from get_pivot_wave import get_pivot_wave

filt = 'F1130W'
cutoff = 0.002
consub_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/'

k_method_file = consub_path + f'k_method/SextansA_{filt}_k_method_pah.fits'
k_method_hdu = fits.open(k_method_file)
k_method_data = k_method_hdu[0].data

interp_method_file = consub_path + f'interp_method/SextansA_{filt}_interp_method_pah.fits'
interp_method_hdu = fits.open(interp_method_file)
interp_method_data = interp_method_hdu[0].data

if filt == 'F335M':
    S23_method_file = consub_path + f'S23_method/SextansA_{filt}_S23_method_pah.fits'
    S23_method_hdu = fits.open(S23_method_file)
    S23_method_data = S23_method_hdu[0].data
    
    L20_method_file = consub_path + f'S23_method/SextansA_{filt}_lai20_method_pah.fits'
    L20_method_hdu = fits.open(L20_method_file)
    L20_method_data = L20_method_hdu[0].data
    

# define the pixel region to use 
# see the box defined in the region file 
x1 = 793
x2 = 948
y1 = 651
y2 = 798

mask = np.zeros_like(k_method_data, dtype = bool)
mask[y1:y2, x1:x2] = True

ind = np.where(k_method_data[mask] > cutoff)[0]

print('{:d} pixels are kept out of {:d} pixels masked, {:5.4f} fraction'.format(len(ind), len(k_method_data[mask]),len(ind)/len(k_method_data[mask]) ))

k_mask = k_method_data[mask][ind]
interp_mask = interp_method_data[mask][ind]

if filt == 'F335M':
    S23_mask = S23_method_data[mask][ind]
    
    L20_mask = L20_method_data[mask][ind]

plt.figure(1)
plt.clf()
plt.scatter(k_mask, interp_mask, c = 'cornflowerblue', s = 30, alpha = 0.7, label = 'Interp Method')

if filt == 'F335M':
    plt.scatter(k_mask, S23_mask, c = 'tomato', s = 30, alpha = 0.7, label = 'S23 Method')
    plt.scatter(interp_mask, S23_mask, c = 'indigo', s = 30, alpha = 0.7, label = 'Interp and S23 Method')
    
    plt.scatter(k_mask, L20_mask, c = 'orange', s = 30, alpha = 0.7, label = 'L20 Method')
    plt.scatter(interp_mask, L20_mask, c = 'pink', s = 30, alpha = 0.7, label = 'Interp and L20 Method')
    plt.scatter(S23_mask, L20_mask, c = 'yellow', s = 30, alpha = 0.7, label = 'S23 and L20 Method')

plt.legend(loc = 'best')

# draw y = x line
xvals = np.linspace(0, 1.5, 100)
plt.plot(xvals, xvals, c = 'k', lw = 0.7)

plt.xlabel(filt + ' PAH from k method', size = 'large')
plt.ylabel(filt + ' PAH from other method (see legend)', size = 'large')

plt.xlim(0, 1.4)
plt.ylim(0, 1.4)


savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/compare/'
savename = f'SextansA_consub_compare_methods_k_method_{filt}'
# plt.savefig(savepath + savename + '.pdf')    

plt.figure(2)
plt.clf()
plt.scatter(k_mask, k_mask/interp_mask, c = 'cornflowerblue', s = 30, alpha = 0.7, label = 'k/interp')

plt.legend(loc = 'best')

# draw y = x line
# xvals = np.linspace(0, 0.2, 100)
# plt.plot(xvals, xvals, c = 'k', lw = 0.7)

plt.xlabel(filt + ' PAH from k method', size = 'large')
plt.ylabel(filt + ' Ratio', size = 'large')

plt.xlim(0, 0.2)
plt.ylim(0, 0.2)


    
    
    
