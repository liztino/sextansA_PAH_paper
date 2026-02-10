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
cutoff = 0.003
consub_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/maps/'

k_method_file = consub_path + f'k_method/SextansA_{filt}_k_method_pah.fits'
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

mask = np.logical_and(mask1, mask2)


for filt_dict in filt_list:  

    filt_dict['mask_arr'] = ma.masked_array(filt_dict['data'], mask = ~mask, fill_value = np.nan)
    filt_dict['mask'] = filt_dict['mask_arr'].flatten()

# also grab the PAH continuum subtracted data
pah_3 = {'name': 'F335M'}
pah_7 = {'name': 'F770W'}
pah_11 = {'name': 'F1130W'}

pah_filt_dict = [pah_3 , pah_7, pah_11]

for pah_dict in pah_filt_dict:
    filt = pah_dict['name']
    
    k_method_file = consub_path + f'k_method/SextansA_{filt}_k_method_pah.fits'
    k_method_hdu = fits.open(k_method_file)
    k_method_data_filt = k_method_hdu[0].data
    
    pah_dict['data'] = k_method_data_filt
    pah_dict['mask_arr'] = ma.masked_array(pah_dict['data'], mask = ~mask, fill_value = np.nan)
    pah_dict['mask_arr'][pah_dict['mask_arr'] < 0] = np.nan
    pah_dict['mask'] = pah_dict['mask_arr'].flatten()


zoom = False

#vals for reg1
zx1 = 810
zx2 = 872
zy1 = 670
zy2 = 720
reg_name = 'reg1'

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
im = ax.imshow(pah_11['mask_arr']/pah_3['mask_arr'], origin = 'lower', cmap = 'Spectral_r', vmax = 15, vmin = 5)
cb = plt.colorbar(im)
cb.set_label('11.3/3.3', size = 'x-large')
plt.xlabel('RA')
plt.ylabel('DEC')

# contours for the PAH data
levels = [0.005, 0.02, 0.05, 0.1, 0.15]
ax.contour(pah_3['data'], levels = levels, colors = 'grey', transform = ax.get_transform(w))

# # contours for the Halpha
# levels = [500, 1500, 2000, 2500]
# ax.contour(Halpha_data, levels = levels, colors = 'k', transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10)

# # scatter plot with the massive stars from Lorenzo+ 2022
# ax.scatter(c.ra.deg, c.dec.deg, marker = 'x', c = 'm', transform = ax.get_transform('world'), s = 30, zorder = 30)

vis.add_beam(ax = ax, header = None, major = F1500W_fwhm,  minor = F1500W_fwhm, angle = beam_angle, color = 'k', frame = True)

if zoom:
    plt.xlim(zx1, zx2)
    plt.ylim(zy1, zy2)
else:
    plt.xlim(x1, x2)
    plt.ylim(y1, y2)
    
savename = 'SextansA_11_3_PAH_band_ratio_{:s}_with_Halpha_stars.pdf'.format(reg_name)
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)

# first make band maps
fig2 = plt.figure(2, figsize = (9,6))
plt.clf()
ax = fig2.add_subplot(projection = w)
im = ax.imshow(pah_11['mask_arr']/pah_7['mask_arr'], origin = 'lower', cmap = 'Spectral_r', vmax = 5, vmin = 2)
cb = plt.colorbar(im)
cb.set_label('11.3/7.7', size = 'x-large')

levels = [0.005, 0.02, 0.05, 0.1, 0.2]
ax.contour(pah_7['data'], levels = levels, colors = 'k', transform = ax.get_transform(w))

# # contours for the Halpha
# levels = [500, 1500, 2000, 2500]
# ax.contour(Halpha_data, levels = levels, colors = 'k', transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10)

# scatter plot with the massive stars from Lorenzo+ 2022
# ax.scatter(c.ra.deg, c.dec.deg, marker = 'x', c = 'm', transform = ax.get_transform('world'), s = 30, zorder = 30)


vis.add_beam(ax = ax, header = None, major = F1500W_fwhm,  minor = F1500W_fwhm, angle = beam_angle, color = 'k', frame = True)

if zoom:
    plt.xlim(zx1, zx2)
    plt.ylim(zy1, zy2)
else:
    plt.xlim(x1, x2)
    plt.ylim(y1, y2)

savename = 'SextansA_11_7_PAH_band_ratio_{:s}_no_contours.pdf'.format(reg_name)
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)

# first make band maps
fig3 = plt.figure(3, figsize = (9,6))
plt.clf()
ax = fig3.add_subplot(projection = w)
im = ax.imshow(pah_7['mask_arr']/pah_3['mask_arr'], origin = 'lower', cmap = 'Spectral_r', vmax = 5, vmin = 1)
cb = plt.colorbar(im)
cb.set_label('7.7/3.3', size = 'x-large')

levels = [0.005, 0.02, 0.05, 0.1, 0.2]
ax.contour(pah_3['data'], levels = levels, colors = 'k', transform = ax.get_transform(w))

# # contours for the Halpha
# levels = [500, 1500, 2000, 2500]
# ax.contour(Halpha_data, levels = levels, colors = 'k', transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10)

# # scatter plot with the massive stars from Lorenzo+ 2022
# ax.scatter(c.ra.deg, c.dec.deg, marker = 'x', c = 'm', transform = ax.get_transform('world'), s = 30, zorder = 30)


vis.add_beam(ax = ax, header = None, major = F1500W_fwhm,  minor = F1500W_fwhm, angle = beam_angle, color = 'k', frame = True)

if zoom:
    plt.xlim(zx1, zx2)
    plt.ylim(zy1, zy2)
else:
    plt.xlim(x1, x2)
    plt.ylim(y1, y2)

savename = 'SextansA_7_3_PAH_band_ratio_{:s}_with_Halpha_stars.pdf'.format(reg_name)
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)

# now make a scatter plot of all the points
norm_val = matplotlib.colors.LogNorm(vmin = 0.002, vmax = 0.15)

plt.figure(4)
plt.clf()

im = plt.scatter(pah_3['mask']/pah_7['mask'], pah_3['mask']/pah_11['mask'], c = pah_3['mask'], cmap = 'magma_r', norm = norm_val, s = 3, alpha = 0.5)
cax = plt.colorbar(im)
cax.set_label('3.3 PAH Flux', size = 'x-large')

plt.xlabel('3.3/7.7 PAH', size = 'x-large')
plt.ylabel('3.3/11.3 PAH', size = 'x-large')
plt.minorticks_on()

plt.semilogx()
plt.semilogy()

plt.xlim(0.1, 4)
plt.ylim(0.03, 1.2)

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/'
savename = 'SextansA_PAH_band_ratio_scatter.pdf'
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)


# band ratio of non-continuum subtracted data to test
F1130W_F335M = F1130W['data']/F335M['data']

fig5 = plt.figure(5, figsize = (9,6))
plt.clf()
ax5 = fig5.add_subplot(projection = w)
# plt.clf()
im = ax5.imshow(F1130W_F335M, origin = 'lower', cmap = 'Spectral_r', vmax = 10, vmin = 1)
cb = plt.colorbar(im)
cb.set_label('F1130W/F335M', size = 'x-large')
ax5.set_xlabel('RA')
ax5.set_ylabel('DEC')

# contours for the PAH data
levels = [0.005, 0.02, 0.05, 0.1, 0.15]
ax5.contour(pah_3['data'], levels = levels, colors = 'grey', transform = ax5.get_transform(w))

plt.xlim(x1, x2)
plt.ylim(y1, y2)

F1130W_F1130W = F1130W['data']/F770W['data']

savename = 'SextansA_F1130W_F770W_ratio.fits'
fits.writeto(savepath + savename, F1130W_F1130W, header = F770W['header'], overwrite = True)
