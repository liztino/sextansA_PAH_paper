#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculates Sigma PAH in the Shi+ region to create Sigma PAH/TIR ratio

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import scipy.interpolate as inter

from astropy.wcs import WCS
import numpy.ma as ma
import astropy.visualization.wcsaxes as vis
import matplotlib
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm

from regions.core import PixCoord
from regions.shapes.rectangle import RectanglePixelRegion
from astropy import wcs
from regions import Regions



# custom functions 
from get_pivot_wave import get_pivot_wave
import k_eq

# open up the continuum subtracted PAH maps 

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
cutoff = -10
consub_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/plots/'

k_method_file = consub_path + f'SextansA_F335M_k_method_pah_k_4.45.fits'
k_method_hdu = fits.open(k_method_file)
k_method_data = k_method_hdu[0].data
k_method_head = k_method_hdu[0].header

# define the pixel region to use 
# see the box defined in the region file 
x1 = 793
x2 = 948
y1 = 651
y2 = 798

mask1 = np.zeros_like(k_method_data, dtype = bool)
mask2 = np.zeros_like(k_method_data, dtype = bool)

mask1[y1:y2, x1:x2] = True

ind = np.where(k_method_data > cutoff)
mask2[ind] = True

# mask = np.logical_and(mask1, mask2)
mask = mask2

# get WCS info 
w = WCS(k_method_head)
pix_area = wcs.utils.proj_plane_pixel_area(w) * 3600**2 / (4.25e10)

# region file info 
reg_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/reg/'
reg_name = 'Shi2014_dust_sf_reg3'

shi_reg = Regions.read(reg_path + reg_name + '.reg')[0]
shi_reg_pix = shi_reg.to_pixel(w)

shi_mask = shi_reg_pix.to_mask(mode = 'center')



for filt_dict in filt_list:  
    
    mask_data = shi_mask.cutout(filt_dict['data'])
    weigh_data = shi_mask.multiply(filt_dict['data'])
    
    cutoff_mask = np.zeros_like(weigh_data, dtype = bool)
    cutoff = 0.020477*3
    ind = np.where(weigh_data > cutoff)
    cutoff_mask[ind] = True
    
    # remove the galaxy by hand
    x1 = 185
    x2 = 270
    y1 = 13
    y2 = 85
    
    weigh_data[y1:y2, x1:x2] = np.nan
    
    
    # filt_dict['mask_arr'] = ma.masked_array(filt_dict['data'], mask = ~mask, fill_value = np.nan)
    # filt_dict['mask'] = filt_dict['mask_arr'].flatten()
    
    # # create a cutoff to sum up on 
    # # cutoff_mask = np.zeros_like(weigh_data, dtype = bool)
    # # cutoff = 0
    # # ind = np.where(weigh_data > pah_dict['cutoff'])
    # # cutoff_mask[ind] = True
    
    # weigh_data = filt_dict['mask_arr']

    
    # # create a cutoff to sum up on 
    # cutoff_mask = np.zeros_like(weigh_data, dtype = bool)
    # cutoff = 0.020477
    # ind = np.where(weigh_data > cutoff)
    # cutoff_mask[ind] = True
    
    filt_dict['reg_sum_arr'] = ma.masked_array(weigh_data, mask = ~cutoff_mask, fill_value = np.nan)
    
    filt_dict['reg_sum_arr'] = filt_dict['reg_sum_arr']  * pix_area

    filt_dict['reg_sum'] = np.nansum(filt_dict['reg_sum_arr'])

    # filt_dict['reg_sum_arr'] = ma.masked_array(weigh_data, fill_value = np.nan)
    
    # remove the sr from pixels so I can add them later
    # filt_dict['reg_sum_arr'] = filt_dict['reg_sum_arr'] * pix_area
    
    
    # filt_dict['reg_sum'] = np.nansum(filt_dict['reg_sum_arr'])
    
    if filt_dict['name'] == 'F1500W':
    
        plt.figure()
        plt.clf()
        plt.imshow(filt_dict['reg_sum_arr'], origin = 'lower', vmin = 0, vmax = 0.2* pix_area)
        plt.title('F1500W')
        
        # plt.xlim(x1, x2)
        # plt.ylim(y1, y2)
        
        figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
        figname = 'SextansA_calc_Sigma_PAH_F1500W.png'
        plt.savefig(figpath + figname, dpi = 300)
    


# also grab the PAH continuum subtracted data
pah_3 = {'name': 'F335M', 'k': 2.07}
pah_7 = {'name': 'F770W', 'k': 4.33}
pah_11 = {'name': 'F1130W', 'k': 7.21}

pah_filt_dict = [pah_3 , pah_7, pah_11]

for pah_dict in pah_filt_dict:
    filt = pah_dict['name']
    
    k_method_file = consub_path + 'SextansA_{:s}_k_method_pah_k_{:3.2f}.fits'.format(filt, pah_dict['k'])
    k_method_hdu = fits.open(k_method_file)
    k_method_data_filt = k_method_hdu[0].data
    
    pah_dict['data'] = k_method_data_filt
    pah_dict['mask_arr'] = ma.masked_array(pah_dict['data'], mask = ~mask, fill_value = np.nan)
    pah_dict['mask_arr'][pah_dict['mask_arr'] < 0] = np.nan
    pah_dict['mask'] = pah_dict['mask_arr'].flatten()
    
    
# get WCS info 
w = WCS(k_method_head)

# uncertainties
pah3_sigma = 0.0036
pah7_sigma =  0.0066
pah11_sigma = 0.0138

pah_3['cutoff'] = pah3_sigma*3
pah_7['cutoff'] = pah7_sigma*3
pah_11['cutoff'] = pah11_sigma*3

pix_area = wcs.utils.proj_plane_pixel_area(w) * 3600**2 / (4.25e10)

for pah_dict in pah_filt_dict:

    shi_mask = shi_reg_pix.to_mask(mode = 'center')
    
    pah_dict['mask_arr'] = ma.masked_array(pah_dict['data'], mask = ~mask, fill_value = np.nan)
    
    mask_data = shi_mask.cutout(pah_dict['mask_arr'] )
    weigh_data = shi_mask.multiply(pah_dict['mask_arr'] )
    
    print(weigh_data)
    
    # pah_dict['mask_arr'] = ma.masked_array(pah_dict['data'], mask = ~mask, fill_value = np.nan)
    # pah_dict['mask'] = pah_dict['mask_arr'].flatten()
    
    # weigh_data = pah_dict['mask_arr']
    
    # remove the galaxy by hand
    gx1 = 185
    gx2 = 235
    gy1 = 13
    gy2 = 64
    
    weigh_data[gy1:gy2, gx1:gx2] = np.nan
    
    # create a cutoff to sum up on 
    cutoff_mask = np.zeros_like(weigh_data, dtype = bool)
    cutoff = pah_dict['cutoff']
    # cutoff = -10
    ind = np.where(weigh_data > cutoff)
    cutoff_mask[ind] = True
    
    pah_dict['reg_sum_arr'] = ma.masked_array(weigh_data, mask = ~cutoff_mask, fill_value = np.nan)
    
    # remove the sr from pixels so I can add them later
    pah_dict['reg_sum_arr'] = pah_dict['reg_sum_arr'] * pix_area
    
    plt.figure()
    plt.clf()
    plt.imshow(pah_dict['reg_sum_arr'], origin = 'lower', vmin = 0, vmax = 0.2* pix_area)
    plt.title(pah_dict['name'])
    # plt.xlim(x1, x2)
    # plt.ylim(y1, y2)
    
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    figname = 'SextansA_calc_Sigma_PAH_{:s}.png'.format(pah_dict['name'])
    # plt.savefig(figpath + figname, dpi = 300)
    
    pah_dict['reg_sum'] = np.nansum(pah_dict['reg_sum_arr'])
    
# convert the sum of PAH from MJy/sr to erg/s/cm2/Hz

# first take care of the sr
# in arcseconds
reg_r = 16
# arcsec^2 in a sr 
sr = 4.25e10

# in arcsecond^2
area = np.pi * reg_r**2

area_sr = area/sr

# distance to Sextans A
# in Mpc, from McQuinn+ 2017 using TRGB
d = 1.46

# 1 pc to cm 
pc = 3.086e18

d_cm = d * 1e6 * pc

# speed of light in microns
c = 3e14

# luminosity of the sun in erg/s
L_sol = 3.846e33

# convert to erg/s/Hz
def MJy_to_L(MJy):
    return 1e-17 * MJy * 4 * np.pi * d_cm**2

# get the wavelengths in Hz
# function to convert wavelength in microns to frequency in Hz
def um_to_hz(um):
    return c/um

# using wavelength weights from D+L2007
hz_3 = um_to_hz(F335M['wave'])
hz_7 = um_to_hz(F770W['wave'])
hz_11 = um_to_hz(F1130W['wave'])

# add up the PAH sum, taking into account the pivot wavelengths
pah_sum = hz_3*MJy_to_L(pah_3['reg_sum']) + hz_7 * MJy_to_L(pah_7['reg_sum']) + hz_11 * MJy_to_L(pah_11['reg_sum'])
print(pah_sum)
print(pah_sum/L_sol)

# from clump work
sigma_PAH_3 = 6.020359491068104e-06 * 1e-6
sigma_PAH_7 = 1.8531946634664654e-05 * 1e-6
sigma_PAH_11 = 5.0678663443125454e-05 * 1e-6

pah_sum = hz_3*MJy_to_L(sigma_PAH_3) + hz_7 * MJy_to_L(sigma_PAH_7) + hz_11 * MJy_to_L(sigma_PAH_11)
print(pah_sum)
print(pah_sum/L_sol)


tot_pah = pah_3['reg_sum'] + pah_7['reg_sum'] + pah_11['reg_sum']


print(tot_pah)
print(F1500W['reg_sum'])
print(tot_pah/F1500W['reg_sum'])

