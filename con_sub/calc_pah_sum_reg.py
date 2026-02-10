#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Takes arbitrarily defined sub regions and integrates the PAH flux.
Saves in a text file to use later 

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
from astropy import wcs

plt.ion()

# custom functions 
from get_pivot_wave import get_pivot_wave
import k_eq

def load_filter(filt):
    
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
    
    # load the middle filter we will be continuum subtracting from
    filename = filepath + f'{filt}_reproject_to_F1500W_rot'
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


# load the continuum subtracted PAH map to make a mask
# using the F1130W to include the most flux
filt = 'F770W'
cutoff = 0.00
consub_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/maps/'

k_method_file = consub_path + f'k_method/SextansA_{filt}_k_method_pah.fits'
k_method_hdu = fits.open(k_method_file)
k_method_data = k_method_hdu[0].data
k_method_head = k_method_hdu[0].header

# define the pixel region to use 
# see the box defined in the region file 

reg_list = np.arange(1, 6)

x1_arr = np.array([810, 857, 900, 882, 850])
x2_arr = np.array([872, 888, 930, 947, 946])
y1_arr = np.array([670, 723, 715, 653, 763])
y2_arr = np.array([720, 760, 746, 715, 795])

pah_3_arr = []
pah_7_arr = []
pah_11_arr = []

# #vals for reg1
# zx1 = 810
# zx2 = 872
# zy1 = 670
# zy2 = 720
# reg_name = 'reg1'

# #vals for reg2
# zx1 = 857
# zx2 = 888
# zy1 = 723
# zy2 = 760
# reg_name = 'reg2'


# #vals for reg3
# zx1 = 900
# zx2 = 930
# zy1 = 715
# zy2 = 746
# reg_name = 'reg3'


# #vals for reg4
# zx1 = 882
# zx2 = 947
# zy1 = 653
# zy2 = 715
# reg_name = 'reg4'


# #vals for reg5
# zx1 = 850
# zx2 = 946
# zy1 = 763
# zy2 = 795
# reg_name = 'reg5'

# orig vals used for big region
# x1 = 793
# x2 = 948
# y1 = 651
# y2 = 798



for i in range(len(reg_list)):
    
    x1 = x1_arr[i]
    x2 = x2_arr[i]
    y1 = y1_arr[i]
    y2 = y2_arr[i]
    
    
    # also grab the PAH continuum subtracted data
    pah_3 = {'name': 'F335M'}
    pah_7 = {'name': 'F770W'}
    pah_11 = {'name': 'F1130W'}
    
    pah_filt_dict = [pah_3 , pah_7, pah_11]
    
    w = WCS(k_method_head)
    
    pix_area = wcs.utils.proj_plane_pixel_area(w) * 3600**2 / (4.25e10)
    
    # define cutoff values for each filter
    pah_3['cutoff'] = 0.005
    pah_7['cutoff'] = 0.005
    pah_11['cutoff'] = 0.005
    
    
    for pah_dict in pah_filt_dict:
        filt = pah_dict['name']
        
        k_method_file = consub_path + f'k_method/SextansA_{filt}_k_method_pah.fits'
        k_method_hdu = fits.open(k_method_file)
        k_method_data_filt = k_method_hdu[0].data
        
        mask1 = np.zeros_like(k_method_data, dtype = bool)
        mask2 = np.zeros_like(k_method_data, dtype = bool)
        
        mask1[y1:y2, x1:x2] = True
        
        ind = np.where(k_method_data > pah_dict['cutoff'])
        mask2[ind] = True
        
        mask = np.logical_and(mask1, mask2)
        
        pah_dict['data'] = k_method_data_filt
        pah_dict['mask_arr'] = ma.masked_array(pah_dict['data'], mask = ~mask, fill_value = np.nan)
        pah_dict['mask_arr'][pah_dict['mask_arr'] < 0] = np.nan
        pah_dict['mask'] = pah_dict['mask_arr'].flatten()
        
        pah_dict['reg_sum'] = np.nansum(pah_dict['mask'] * pix_area)
        
        plt.figure()
        plt.clf()
        plt.imshow(pah_dict['mask_arr'],origin = 'lower')
        plt.xlim(x1, x2)
        plt.ylim(y1,y2)
        plt.title('Reg {:d} PAH {:s}'.format(i, pah_dict['name']))

    
    pah_3_arr.append(pah_3['reg_sum'])
    pah_7_arr.append(pah_7['reg_sum'])
    pah_11_arr.append(pah_11['reg_sum'])
    
ascii_name = ['reg', 'pah_3', 'pah_7', 'pah_11']
ascii_table = [reg_list, pah_3_arr, pah_7_arr, pah_11_arr]  

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/'
savename = 'SextansA_reg_sum_real_cutoff.txt'
ascii.write(ascii_table, savepath + savename, names = ascii_name, overwrite = True)



# plt.figure(1)
# plt.clf()
# plt.imshow(pah_3['mask_arr'], origin = 'lower')

# plt.xlim(zx1, zx2)
# plt.ylim(zy1, zy2)


# for pah_dict in pah_filt_dict:

#     shi_mask = shi_reg_pix.to_mask(mode = 'center')
    
#     mask_data = shi_mask.cutout(pah_dict['data'])
#     weigh_data = shi_mask.multiply(pah_dict['data'])
    
    
#     # remove the galaxy by hand
#     x1 = 185
#     x2 = 235
#     y1 = 13
#     y2 = 64
    
#     weigh_data[y1:y2, x1:x2] = np.nan
    
#     # create a cutoff to sum up on 
#     cutoff_mask = np.zeros_like(weigh_data, dtype = bool)
#     cutoff = 0.1
#     ind = np.where(weigh_data > pah_dict['cutoff'])
#     cutoff_mask[ind] = True
    
#     pah_dict['reg_sum_arr'] = ma.masked_array(weigh_data, mask = ~cutoff_mask, fill_value = np.nan)
    
#     # remove the sr from pixels so I can add them later
#     pah_dict['reg_sum_arr'] = pah_dict['reg_sum_arr'] * pix_area
    
#     # plt.figure()
#     # plt.clf()
#     # plt.imshow(pah_dict['reg_sum'], origin = 'lower')
    
#     pah_dict['reg_sum'] = np.nansum(pah_dict['reg_sum_arr'])


