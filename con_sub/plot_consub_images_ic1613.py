#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates band ratio plots and images from the data


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
from scipy import ndimage
import astropy.convolution.convolve
from astropy.convolution import Gaussian2DKernel


plt.ion()

# custom functions 
from get_pivot_wave import get_pivot_wave
import k_eq


# def load_filter(filt):
    
#     filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
    
#     # load the middle filter we will be continuum subtracting from
#     filename = filepath + f'{filt}_reproject_to_F1500W_rot'
#     hdu = fits.open(filename + '.fits')
#     header = hdu[0].header
#     data = hdu[0].data
#     pivot = get_pivot_wave(filt)
    
#     filt_dict = {'name': filt, 'data': data, 'header': header, 'wave': pivot}
    
#     return filt_dict

# F300M = load_filter('F300M')
# F335M = load_filter('F335M')
# F360M = load_filter('F360M')
# F560W = load_filter('F560W')
# F770W = load_filter('F770W')
# F1000W = load_filter('F1000W')
# F1130W = load_filter('F1130W')
# F1500W = load_filter('F1500W')

# filt_list = [F300M, F335M, F360M, F560W, F770W, F1000W, F1130W, F1500W]


# load the continuum subtracted PAH map to make a mask
# using the F1130W to include the most flux
filt = 'F770W'
cutoff = 0.1
consub_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613/aug23_reduction/quick_k_consub/'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613/aug23_reduction/plots/'

k_method_file = consub_path + f'IC1613_{filt}_k_method_pah.fits'
k_method_hdu = fits.open(k_method_file)
k_method_data = k_method_hdu[0].data
k_method_head = k_method_hdu[0].header

# define the pixel region to use 
# see the box defined in the region file 
# for the original box
x1 = 482
x2 = 667
y1 = 646
y2 = 900

# # slightly larger box
# x1 = 737
# x2 = 1020
# y1 = 618
# y2 = 890

mask1 = np.zeros_like(k_method_data, dtype = bool)
mask2 = np.zeros_like(k_method_data, dtype = bool)

mask1[y1:y2, x1:x2] = True

ind = np.where(k_method_data > cutoff)
mask2[ind] = True

mask = np.logical_and(mask1, mask2)


# for filt_dict in filt_list:  

#     filt_dict['mask_arr'] = ma.masked_array(filt_dict['data'], mask = ~mask, fill_value = np.nan)
#     filt_dict['mask'] = filt_dict['mask_arr'].flatten()

# also grab the PAH continuum subtracted data
pah_3 = {'name': 'F335M'}
pah_7 = {'name': 'F770W'}
pah_11 = {'name': 'F1130W'}

pah_filt_dict = [ pah_7]

for pah_dict in pah_filt_dict:
    filt = pah_dict['name']
    
    k_method_file = consub_path + f'IC1613_{filt}_k_method_pah.fits'
    k_method_hdu = fits.open(k_method_file)
    k_method_data_filt = k_method_hdu[0].data
    
    pah_dict['data'] = k_method_data_filt
    pah_dict['mask_arr'] = ma.masked_array(pah_dict['data'], mask = ~mask, fill_value = np.nan)
    pah_dict['mask_arr'][pah_dict['mask_arr'] < 0] = np.nan
    pah_dict['mask'] = pah_dict['mask_arr'].flatten()


zoom = True

#vals for reg1
zx1 = 330
zx2 = 761
zy1 = 317
zy2 = 532
reg_name = 'reg1'

# adding ancillary data onto the images
anc_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
Halpha_file = 'MUSE_SexA_pipeline_Halpha.fits'

Halpha_hdu = fits.open(anc_path + Halpha_file)
Halpha_data = Halpha_hdu[0].data
Halpha_header = Halpha_hdu[0].header
Halpha_w = WCS(Halpha_header)

cat_file = 'Lorenzo_2022_SexA_spec.csv'
cat = ascii.read(anc_path + cat_file)
c = SkyCoord(cat['RA (J2000.0)'], cat['Dec. (J2000.0)'], frame = 'fk5', unit = (u.hourangle, u.deg))

#  YSOs
star_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/dolphot/reg/'
yso_file = 'YSO_line_text.txt'
yso = ascii.read(star_path + yso_file)

anc_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
HI_file = 'DDO75_R_X0_P_R_FLATTEN.FITS'
HI_hdu = fits.open(anc_path + HI_file)[0]
HI_data = HI_hdu.data
HI_header = HI_hdu.header
HI_w = WCS(HI_header)

if not zoom: 
    reg_name = 'full'

# PSF/beam info 
# in arcseconds
F1500W_fwhm = (0.488 * u.arcsec).to(u.deg)
beam_angle = 0 * u.deg

w = WCS(k_method_head)

# first make band maps
fig1 = plt.figure(1, figsize = (9,6))
plt.clf()
ax = fig1.add_subplot(projection = w)
# plt.clf()

custom_map  = sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, reverse=False, as_cmap=True)

norm = simple_norm(pah_7['mask_arr'], 'asinh', min_cut = 0.002, max_cut = 0.4)

im = ax.imshow(pah_7['data'], origin = 'lower', cmap = custom_map,  norm = norm)
cb = plt.colorbar(im)
cb.set_label(r'MJy sr$^{-1}$', size = 'x-large')
plt.xlabel('RA')
plt.ylabel('DEC')


# smooth data for contouors 
# use `nearest` in stead of others
# resample your data grid by a factor of 4 using cubic spline interpolation
k = Gaussian2DKernel(2)
pah_sm = astropy.convolution.convolve(pah_7['data'], k, boundary='extend')

# contours for the PAH data
levels = [0.02, 0.08, 0.15, 0.2,  0.5,1]
ax.contour(pah_sm , levels = levels, colors = 'grey', transform = ax.get_transform(w), linewidths =0.7)

# # contours for the Halpha
# levels = [500, 1500, 2000, 2500]
# ax.contour(Halpha_data, levels = levels, colors = 'k', transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10)

# # scatter plot with the massive stars from Lorenzo+ 2022
# ax.scatter(c.ra.deg, c.dec.deg, marker = 'x', c = 'm', transform = ax.get_transform('world'), s = 30, zorder = 30)

vis.add_beam(ax = ax, header = None, major = F1500W_fwhm,  minor = F1500W_fwhm, angle = beam_angle, color = 'k', frame = True)
vis.add_scalebar(ax = ax, length = 1 * u.arcsec, corner = 'bottom right', label = '1" = 3.5 pc', color = 'k')

# YSOs

# plt.scatter(yso['RA'], yso['DEC'], marker = 'x', c = 'green', transform = ax.get_transform('world'), s = 30, zorder = 30)

# HI contours
# levels = [100, 150, 200]
# ax.contour(HI_data, levels = levels, colors = 'k', transform = ax.get_transform(HI_w), origin = 'lower', zorder = 10)


if zoom:
    plt.xlim(zx1, zx2)
    plt.ylim(zy1, zy2)
else:
    plt.xlim(x1, x2)
    plt.ylim(y1, y2)
    
savename = 'IC1613_pah_7_quick_k_method_zoom.pdf'.format(reg_name)
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)

