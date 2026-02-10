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
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

pah3_sigma = 0.0036
pah7_sigma =  0.0066
pah11_sigma = 0.0138

# dictionary for all props
pah_3 = {'name': 'F335M', 'sigma': pah3_sigma, 'k': 2.07}
pah_7 = {'name': 'F770W', 'sigma': pah7_sigma, 'k': 4.33}
pah_11 = {'name': 'F1130W', 'sigma': pah11_sigma, 'k': 7.21}

# load the continuum subtracted PAH map to make a mask
# using the F1130W to include the most flux
filt = 'F770W'
cutoff = pah7_sigma
F770W_k = 4.33
consub_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/plots/'

k_method_file = consub_path + 'k_method_new/SextansA_F770W_k_method_pah_k_4.33.fits'
k_method_hdu = fits.open(k_method_file)
k_method_data = k_method_hdu[0].data
k_method_head = k_method_hdu[0].header

# define the pixel region to use 
# see the box defined in the region file 
pix_diff = 145
x1 = 800
x2 = x1 + pix_diff
y1 = 651
y2 = y1 + pix_diff


# # original-- do not delete!!!
# x1 = 793
# x2 = 948
# y1 = 651
# y2 = 794

# # gets clump to the right
# x1 = 813
# x2 = 990
# y1 = 651
# y2 = 794

# for the ALMA proposal box
# x1 = 798
# x2 = 948
# y1 = 651
# y2 = 763

# slightly larger box
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


for filt_dict in filt_list:  

    filt_dict['mask_arr'] = ma.masked_array(filt_dict['data'], mask = ~mask, fill_value = np.nan)
    filt_dict['mask'] = filt_dict['mask_arr'].flatten()


pah_filt_dict = [pah_3 , pah_7, pah_11]

for pah_dict in pah_filt_dict:
    filt = pah_dict['name']
    k = pah_dict['k']
    
    k_method_file = consub_path + 'k_method_new/SextansA_{:s}_k_method_pah_k_{:3.2f}.fits'.format(filt, k)
    k_method_hdu = fits.open(k_method_file)
    k_method_data_filt = k_method_hdu[0].data
    
    pah_dict['data'] = k_method_data_filt
    pah_dict['mask_arr'] = ma.masked_array(pah_dict['data'], mask = ~mask, fill_value = np.nan)
    pah_dict['mask_arr'][pah_dict['mask_arr'] < 0] = np.nan
    pah_dict['mask'] = pah_dict['mask_arr'].flatten()


# do mask_arr just for F1500W to plot it later
F1500W['mask_arr'] = ma.masked_array(F1500W['data'], mask = ~mask, fill_value = np.nan)
F1500W['mask_arr'][F1500W['mask_arr'] < 0] = np.nan

zoom = False

#vals for reg1
zx1 = 812
zx2 = 870
zy1 = 670
zy2 = 720
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

plt.style.use('default')

norm = simple_norm(pah_3['mask_arr'], 'asinh', min_cut = 0.002, max_cut = np.nanmax(pah_3['mask_arr']))
im = ax.imshow(pah_3['data'], origin = 'lower', cmap = 'Greys', norm = norm)

cb = plt.colorbar(im, orientation='vertical', fraction=0.046, pad=0.02)
cb.set_label(r'MJy sr$^{-1}$', size = 'x-large')
ax.set_xlabel('RA')
ax.set_ylabel('DEC')

# contours for the PAH data
unc_3 = pah_3['sigma']
levels = [unc_3*3, unc_3*5, unc_3*10, unc_3*20]
# levels = [0.005, 0.01, 0.02, 0.05, 0.1, 0.15]
ax.contour(pah_3['data'], levels = levels, colors = 'grey', transform = ax.get_transform(w), linewidths =0.7)

# # contours for the Halpha
# levels = [500, 1500, 2000, 2500]
# ax.contour(Halpha_data, levels = levels, colors = 'k', transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10)

# # scatter plot with the massive stars from Lorenzo+ 2022
# ax.scatter(c.ra.deg, c.dec.deg, marker = 'x', c = 'm', transform = ax.get_transform('world'), s = 30, zorder = 30)

vis.add_beam(ax = ax, header = None, major = F1500W_fwhm,  minor = F1500W_fwhm, angle = beam_angle, color = 'k', frame = True)
vis.add_scalebar(ax = ax, length = 1 * u.arcsec, corner = 'bottom right', label = '1\" = 7 pc')

# YSOs

plt.scatter(yso['RA'], yso['DEC'], marker = 'D', c = 'purple', transform = ax.get_transform('world'), s = 30, zorder = 30)

# HI contours
# levels = [100, 150, 200]
# ax.contour(HI_data, levels = levels, colors = 'k', transform = ax.get_transform(HI_w), origin = 'lower', zorder = 10)


if zoom:
    ax.set_xlim(zx1, zx2)
    ax.set_ylim(zy1, zy2)
else:
    ax.set_xlim(x1, x2)
    ax.set_ylim(y1, y2)
    
# savename = 'SextansA_PAH_new_k_method_3_PAH_k_{:3.2f}.pdf'.format(pah_3['k'])
savename = 'SextansA_PAH_new_k_method_3_PAH_k_{:3.2f}_yso_purple.png'.format(pah_3['k'])
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 500)

crest = sns.color_palette("magma_r", as_cmap=True)
custom_map  = sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, reverse=False, as_cmap=True)

# first make band maps
fig2 = plt.figure(2, figsize = (9,6))
plt.clf()
ax = fig2.add_subplot(projection = w)
norm = simple_norm(pah_7['mask_arr'], 'asinh', min_cut = 0.002, max_cut = np.nanmax(pah_7['mask_arr']))
im = ax.imshow(pah_7['data'], origin = 'lower', cmap =  'magma_r', norm= norm)
cb = plt.colorbar(im, orientation='vertical', fraction=0.046, pad=0.02)
cb.set_label(r'MJy sr$^{-1}$', size = 'x-large')
plt.xlabel('RA')
plt.ylabel('DEC')

sexa_norm = norm 

unc_7 = pah_7['sigma']
levels = [unc_7*3, unc_7*5, unc_7*10, unc_7*20]
# levels = [0.006, 0.02, 0.05, 0.1, 0.2]
ax.contour(pah_7['data'], levels = levels, colors = 'gray', transform = ax.get_transform(w), linewidths =0.7)

# # contours for the Halpha
levels = [500, 1500, 2000, 2500, 5000]
# ax.contourf(Halpha_data, levels = levels, colors = 'k', transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10)

cmap = sns.color_palette("bone_r", as_cmap=True)
# ax.contourf(Halpha_data, levels = levels, cmap = cmap, transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10, alpha = 0.3)
# ax.contour(Halpha_data, levels = levels, colors = 'k', linewidths = 0.5, transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10, alpha = 0.3)

# scatter plot with the massive stars from Lorenzo+ 2022
# ax.scatter(c.ra.deg, c.dec.deg, marker = 'x', c = 'b', transform = ax.get_transform('world'), s = 30, zorder = 30)


vis.add_beam(ax = ax, header = None, major = F1500W_fwhm,  minor = F1500W_fwhm, angle = beam_angle, color = 'k', frame = True)
vis.add_scalebar(ax = ax, length = 1 * u.arcsec, corner = 'bottom right', label = '1\" = 7 pc')

# YSOs
# plt.scatter(yso['RA'], yso['DEC'], marker = 'x', c = 'red', transform = ax.get_transform('world'), s = 30, zorder = 30)

# # HI contours
# levels = [100, 150, 200]
# ax.contour(HI_data, levels = levels, colors = 'k', transform = ax.get_transform(HI_w), origin = 'lower', zorder = 10)

if zoom:
    plt.xlim(zx1, zx2)
    plt.ylim(zy1, zy2)
else:
    plt.xlim(x1, x2)
    plt.ylim(y1, y2)

# savename = 'SextansA_PAH_new_k_method_7_PAH_k_{:3.2f}.pdf'.format( pah_7['k'])
savename = 'SextansA_PAH_new_k_method_7_PAH_k_{:3.2f}.png'.format( pah_7['k'])
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 500)

# first make band maps
fig3 = plt.figure(3, figsize = (9,6))
plt.clf()
# plt.style.use('dark_background')
ax = fig3.add_subplot(projection = w)
norm = simple_norm(pah_11['mask_arr'], 'asinh', min_cut = 0.02, max_cut = np.nanmax(pah_11['mask_arr']))
im = ax.imshow(pah_11['data'], origin = 'lower', cmap = 'Greys', norm = norm)

cb = plt.colorbar(im, orientation='vertical', fraction=0.046, pad=0.02)
cb.set_label(r'MJy sr$^{-1}$', size = 'x-large')
plt.xlabel('RA')
plt.ylabel('DEC')

unc_11 = pah_11['sigma']
levels = [unc_11*3, unc_11*5, unc_11*10, unc_11*20]
# levels = [0.03, 0.06, 0.1, 0.3, 0.6]
ax.contour(pah_11['data'], levels = levels, colors = 'gray', transform = ax.get_transform(w), linewidths =0.7)

# # contours for the Halpha
# levels = [500, 1500, 2000, 2500]
# ax.contour(Halpha_data, levels = levels, colors = 'k', transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10)

# # scatter plot with the massive stars from Lorenzo+ 2022
# ax.scatter(c.ra.deg, c.dec.deg, marker = 'x', c = 'm', transform = ax.get_transform('world'), s = 30, zorder = 30)


# vis.add_beam(ax = ax, header = None, major = F1500W_fwhm,  minor = F1500W_fwhm, angle = beam_angle, color = 'k', frame = True)
# vis.add_scalebar(ax = ax, length = 1 * u.arcsec, corner = 'bottom right', label = '1\" = 7 pc')

plt.scatter(yso['RA'], yso['DEC'], marker = 'D', c = 'purple', transform = ax.get_transform('world'), s = 30, zorder = 30)


if zoom:
    plt.xlim(zx1, zx2)
    plt.ylim(zy1, zy2)
else:
    plt.xlim(x1, x2)
    plt.ylim(y1, y2)

# savename = 'SextansA_PAH_new_k_method_11_PAH_k_{:3.2f}.pdf'.format(pah_11['k'])
savename = 'SextansA_PAH_new_k_method_11_PAH_k_{:3.2f}_yso_purple.png'.format(pah_11['k'])
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 500)

# add the F1500W filter as a measure of the hot dust
fig4 = plt.figure(4, figsize = (9,6))
plt.clf()
ax = fig4.add_subplot(projection = w)
norm = simple_norm(F1500W['mask_arr'], 'asinh', min_cut = 0.02, max_cut = np.nanmax(F1500W['mask_arr']))
im = ax.imshow(F1500W['data'], origin = 'lower', cmap = 'magma_r', norm = norm)

cb = plt.colorbar(im, orientation='vertical', fraction=0.046, pad=0.02)
cb.set_label(r'MJy sr$^{-1}$', size = 'x-large')
plt.xlabel('RA')
plt.ylabel('DEC')

unc_F1500W = 0.020477
levels = [unc_F1500W*3, unc_F1500W*5, unc_F1500W*10, unc_F1500W*20]
# levels = [0.03, 0.06, 0.1, 0.3, 0.6]
ax.contour(F1500W['data'], levels = levels, colors = 'gray', transform = ax.get_transform(w), linewidths =0.7)

# # contours for the Halpha
# levels = [500, 1500, 2000, 2500]
# ax.contour(Halpha_data, levels = levels, colors = 'k', transform = ax.get_transform(Halpha_w), origin = 'lower', zorder = 10)

# # scatter plot with the massive stars from Lorenzo+ 2022
# ax.scatter(c.ra.deg, c.dec.deg, marker = 'x', c = 'm', transform = ax.get_transform('world'), s = 30, zorder = 30)


vis.add_beam(ax = ax, header = None, major = F1500W_fwhm,  minor = F1500W_fwhm, angle = beam_angle, color = 'k', frame = True)
vis.add_scalebar(ax = ax, length = 1 * u.arcsec, corner = 'bottom right', label = '1\" = 7 pc')


if zoom:
    plt.xlim(zx1, zx2)
    plt.ylim(zy1, zy2)
else:
    plt.xlim(x1, x2)
    plt.ylim(y1, y2)

# savename = 'SextansA_PAH_new_k_method_F1500W_PAH_k.pdf'
savename = 'SextansA_PAH_new_k_method_F1500W_PAH_k.png'
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 500)


