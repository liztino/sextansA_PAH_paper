#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Performs continuum subtraction with the new method that accounts for PAH contamination in the continuum filters

@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
plt.ion()

filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'

k = 0.32

F560W_file =  'F560W_reproject_to_F1500W_rot'
F560W_hdu = fits.open(filepath + F560W_file + '.fits')
F560W_head = F560W_hdu[0].header
F560W_data = F560W_hdu[0].data
F560W_pivot = 5.635

F770W_file =  'F770W_reproject_to_F1500W_rot'
F770W_hdu = fits.open(filepath + F770W_file + '.fits')
F770W_head = F770W_hdu[0].header
F770W_data = F770W_hdu[0].data
F770W_pivot = 7.639

F1000W_file =  'F1000W_reproject_to_F1500W_rot'
F1000W_hdu = fits.open(filepath + F1000W_file + '.fits')
F1000W_head = F1000W_hdu[0].header
F1000W_data = F1000W_hdu[0].data
F1000W_pivot = 9.953


# pah = (F1000W_pivot - F560W_pivot)/ (F1000W_pivot - F560W_pivot - k*F1000W_pivot - k*F770W_pivot) * (
#     (F770W_data - (((F1000W_data * (F770W_pivot - F560W_pivot)))/(F1000W_pivot - F560W_pivot)) + (((F560W_data * (F770W_pivot - F560W_pivot)))/(F1000W_pivot - F560W_pivot))))


def fp2_v2(f1, f2, f3, lam1, lam2, lam3, k):
    # equation calculated through Julia's math
    fp1 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(1 - k - ((lam2 - lam1)/(lam3 - lam1)))
    fp2 = fp1*k
    
    print('fp1', fp1)
    print('fp2', fp2)
    
    return fp2

def fc2_v2(f1, f2, f3, lam1, lam2, lam3, k):
    # equation calculated through Julia's math
    fp1 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(1 - k - ((lam2 - lam1)/(lam3 - lam1)))
    fc1 = f1 - fp1    
    fc2 = (((f3 - fc1)/(lam3 - lam1)) * (lam2 - lam1)) + fc1
    
    print('fc1', fc1)
    print('fc2', fc2)
    print('slope', (((f3 - fc1)/(lam3 - lam1))))
    # print('delta_lam', lam2 - lam1)
    
    return fc2

k = 5.505

pah = fp2_v2(F560W_data, F770W_data, F1000W_data, F560W_pivot, F770W_pivot, F1000W_pivot, k)
con = fc2_v2(F560W_data, F770W_data, F1000W_data, F560W_pivot, F770W_pivot, F1000W_pivot, k)


plt.figure(2)
plt.clf()
plt.imshow(pah, origin = 'lower', vmin = 0, vmax = 0.01)


save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images_rot/'
save_name = 'SextansA_F770W_pah_newmethod.fits'
fits.writeto(save_path + save_name, pah, header = F770W_head, overwrite = True)

save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images_rot/'
save_name = 'SextansA_F770W_con_newmethod.fits'
fits.writeto(save_path + save_name, con, header = F770W_head, overwrite = True)

# open up the original subtraction, the interpolation method
inpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images_rot/'
inname = 'SextansA_F770W_sub_match_rot_pah.fits'
data = fits.open(inpath + inname)[0].data

ratio = pah/data

save_name = 'SextansA_F770W_pah_ratio_methods.fits'
fits.writeto(save_path + save_name, ratio, header = F770W_head, overwrite = True)


