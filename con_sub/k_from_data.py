#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates ratio plots of PAH and continuum filters to test the derived k empirically

@author: etarantino
"""


from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
import numpy.ma as ma
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
# using the F335M one to be most conservative
filt = 'F335M'
cutoff = -20
consub_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/compare/'

k_method_file = consub_path + f'k_method/SextansA_{filt}_k_method_pah.fits'
k_method_hdu = fits.open(k_method_file)
k_method_data = k_method_hdu[0].data

# define the pixel region to use 
# see the box defined in the region file 
# x1 = 793
# x2 = 948
# y1 = 651
# y2 = 798

x1 = 1007
x2 = 1040
y1 = 420
y2 = 470

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
    pah_dict['mask'] = pah_dict['mask_arr'].flatten()



# make ratio plots for each filter combination

####################### F335M #######################
plt.figure(1)
plt.clf()
# im = plt.scatter(F360M['mask']/F300M['mask'], F335M['mask']/F300M['mask'], c = pah_3['mask'], alpha = 0.5, cmap = 'rainbow', vmin = 0.001, vmax = 0.1, s =5)
im = plt.scatter(F360M['mask']/F300M['mask'], F335M['mask']/F300M['mask'], c = F335M['mask'], alpha = 0.5, cmap = 'rainbow', vmin = 0.001, vmax = 0.1, s =5)
cax = plt.colorbar(im)
# cax.set_label('PAH 3.3 Flux from k method')
cax.set_label('F335M flux ')
plt.xlabel('F360M/F300M', size = 'large')
plt.ylabel('F335M/F300M', size = 'large')

def con_line(x, A, B):
    return A + B*x

x_vals = np.linspace(0, 2, 100)
lai_vals = con_line(x_vals, 0.35, 0.65)
plt.plot(x_vals, lai_vals, label = 'L20')

s23_vals = con_line(x_vals, -0.2, 1.6)
plt.plot(x_vals, s23_vals, label = 'S23')

lam = (F335M['wave'] - F300M['wave'])/(F360M['wave'] - F300M['wave'])
con_vals = con_line(x_vals, 1-lam, lam)
plt.plot(x_vals, con_vals, label = 'Continuum \nA={:3.2f} B={:3.2f}'.format(1-lam, lam))

k = 3.8
m = np.argmax(pah_3['mask'])
incpt = F335M['mask'][m]/F300M['mask'][m] - k * F360M['mask'][m]/F300M['mask'][m]
k_vals = con_line(x_vals, incpt, k)
plt.plot(x_vals, k_vals, label = 'k={:2.1f}, incpt={:3.2f}'.format(k, incpt))

plt.legend(loc = 'best')
plt.xlim(0,2)
plt.ylim(0,2)

savename = 'SextansA_con_ratio_F335M_filter_flux_continuum'
plt.savefig(savepath + savename + '.pdf')


####################### F770W #######################
plt.figure(2)
plt.clf()
# im = plt.scatter(F560W['mask']/F1000W['mask'], F770W['mask']/F1000W['mask'], c = pah_7['mask'], alpha = 0.5, cmap = 'rainbow', vmin = 0.001, vmax = 0.1, s =5)
im = plt.scatter(F560W['mask']/F1000W['mask'], F770W['mask']/F1000W['mask'], c = F770W['mask'], alpha = 0.5, cmap = 'rainbow', vmin = 0.001, vmax = 0.1, s =5)
cax = plt.colorbar(im)
# cax.set_label('PAH 7.7 Flux from k method')
cax.set_label('F770W flux')
plt.xlabel('F560W/F1000W', size = 'large')
plt.ylabel('F770W/F1000W', size = 'large')

def con_line(x, A, B):
    return A + B*x

k = 5.84
x_vals = np.linspace(0, 4, 100)

m = np.argmax(pah_7['mask'])
incpt = F770W['mask'][m]/F1000W['mask'][m] - k * F560W['mask'][m]/F1000W['mask'][m]
k_vals = con_line(x_vals, incpt, k)
plt.plot(x_vals, k_vals, label = 'k={:2.1f}, incpt={:3.2f}'.format(k, incpt), c = 'tomato')

lam = (F770W['wave'] - F560W['wave'])/(F1000W['wave'] - F560W['wave'])
con_vals = con_line(x_vals, 1-lam, lam)
plt.plot(x_vals, con_vals, label = 'Continuum \nA={:3.2f} B={:3.2f}'.format(1-lam, lam), c = 'cornflowerblue')

plt.legend(loc = 'best')
plt.xlim(0,4)
plt.ylim(0,4)

savename = 'SextansA_con_ratio_F770W_filter_flux_continuum'
plt.savefig(savepath + savename + '.pdf')



####################### F1130W #######################
plt.figure(3)
plt.clf()
# im = plt.scatter(F1500W['mask']/F1000W['mask'], F1130W['mask']/F1000W['mask'], c = pah_11['mask'], alpha = 0.5, cmap = 'rainbow', vmin = 0.001, vmax = 0.1, s =5)
im = plt.scatter(F1500W['mask']/F1000W['mask'], F1130W['mask']/F1000W['mask'], c = F1130W['mask'], alpha = 0.5, cmap = 'rainbow', vmin = 0.001, vmax = 0.1, s =5)
cax = plt.colorbar(im)
cax.set_label('PAH 11.3 Flux from k method')
cax.set_label('F1130W')
plt.xlabel('F1500W/F1000W', size = 'large')
plt.ylabel('F1130W/F1000W', size = 'large')

def con_line(x, A, B):
    return A + B*x

k = 5.24
x_vals = np.linspace(0, 10, 100)
m = np.argmax(pah_11['mask'])
incpt = F1130W['mask'][m]/F1000W['mask'][m] - k * F1500W['mask'][m]/F1000W['mask'][m]
k_vals = con_line(x_vals, incpt, k)
plt.plot(x_vals, k_vals, label = 'k={:2.1f}, incpt={:3.2f}'.format(k, incpt), c = 'tomato')

lam = (F1130W['wave'] - F1000W['wave'])/(F1500W['wave'] - F1000W['wave'])
con_vals = con_line(x_vals, 1-lam, lam)
plt.plot(x_vals, con_vals, label = 'Continuum \nA={:3.2f} B={:3.2f}'.format(1-lam, lam), c='cornflowerblue')

plt.legend(loc = 'best')
plt.xlim(0, 10)
plt.ylim(0, 10)

savename = 'SextansA_con_ratio_F1130W_filter_flux_continuum'
plt.savefig(savepath + savename + '.pdf')

# make plot of the outliers
y, x = ma.where(F560W['mask_arr']/F1000W['mask_arr'] > 1.5)
mask_data = np.zeros_like(F770W['mask_arr'], dtype = bool)
mask_data[y,x] = True

copy_data = np.copy(F770W['mask_arr'])
copy_data[~mask_data] = np.nan

plt.figure(4)
plt.clf()
plt.imshow(copy_data,  origin = 'lower', vmin = 0, vmax = 0.1)
plt.xlim(x1, x2)
plt.ylim(y1, y2)

plt.figure(5)
plt.clf()
plt.imshow(F770W['mask_arr'], origin = 'lower', vmin = 0, vmax = 0.1)
plt.xlim(x1, x2)
plt.ylim(y1, y2)


y, x = ma.where((F360M['mask_arr']/F300M['mask_arr'] > 1) & (F335M['mask_arr']/F300M['mask_arr'] <1.4))
mask_data = np.zeros_like(F335M['mask_arr'], dtype = bool)
mask_data[y,x] = True

copy_data = np.copy(F335M['mask_arr'])
copy_data[~mask_data] = np.nan

plt.figure(4)
plt.clf()
plt.imshow(copy_data,  origin = 'lower', vmin = 0, vmax = 0.1)
plt.xlim(x1, x2)
plt.ylim(y1, y2)

plt.figure(5)
plt.clf()
plt.imshow(F335M['mask_arr'], origin = 'lower', vmin = 0, vmax = 0.1)
plt.xlim(x1, x2)
plt.ylim(y1, y2)

