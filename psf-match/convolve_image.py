#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Convolves kernel with image to create smoothed image

@author: Eliz
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian1DKernel, Box1DKernel
import scipy.signal
import numpy as n
import scipy.interpolate
import scipy.ndimage
import glob
from congrid import congrid
from regrid import regrid
from matplotlib.colors import LogNorm

plt.ion()

def pad_zeros(pad, array):
    shape = np.shape(array)
    
    pad_val = shape[0] + pad

    new_array = np.zeros((pad_val, pad_val), dtype = float)
    new_shape = np.shape(new_array)
        
    xstart_ker = int((new_shape[0] - shape[0])/2.)
    xend_ker = xstart_ker + shape[0]
    ystart_ker = int((new_shape[1] - shape[1])/2.)
    yend_ker = ystart_ker + shape[1]
    
    new_array[xstart_ker:xend_ker, ystart_ker:yend_ker] = array
    
    return new_array, xstart_ker, xend_ker, ystart_ker, yend_ker

# filt_list = ['F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W' ]
filt_list = ['F200W']
# filt_list = ['F560W',  'F770W', 'F1000W', 'F1130W' ]


convolve_err = False

for filt in filt_list:

    # miri_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/miri_aug_mosaics/'
    # nircam_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/nircam_may_mosaics/'
    
        
    miri_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'
    nircam_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'
    
    # miri_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/old/align_north_old/'
    # nircam_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/old/align_north_old/'
    
    miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
    nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
    
    if filt in miri_filts:
        path = miri_path
        # image = glob.glob(miri_path + '*{}_*skysub_prop_i2d.fits'.format(filt))[0]
        image = glob.glob(miri_path + '*{}_align_north.fits'.format(filt))[0]
    else:
        path = nircam_path
        # image = glob.glob(nircam_path + '*{}_*modb_bkgrsub_i2d.fits'.format(filt.lower()))[0]
        image = glob.glob(miri_path + '*{}_align_north.fits'.format(filt))[0]
    
    
    filepath_ker = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/kernels/rot_ker/'
    filename_ker = glob.glob(filepath_ker + '*{}*'.format(filt))[0]
    
    hdulist_data = fits.open(image)
    
    if convolve_err:
        data = hdulist_data['ERR'].data
        data_header = hdulist_data['SCI'].header
        orig_data = data
        data = data * data
    
    else:
        # data = hdulist_data['SCI'].data
        # data_header = hdulist_data['SCI'].header
        data = hdulist_data[0].data
        data_header = hdulist_data[0].header
    
    nan_map = np.isnan(data)
    data[nan_map] = 0.0
    # err[np.isnan(err)] = 0.0
    
    pad_data = data
    
    hdulist_ker = fits.open(filename_ker)
    kernel = hdulist_ker[0].data
    kernel_header = hdulist_ker[0].header
    
    data_shape = np.shape(data)
    ker_shape = np.shape(kernel)
    
    print(data_shape, ker_shape)
    
    pad = 200
    pad_data, xstart, xend, ystart, yend = pad_zeros(pad, data)
    
    pad_ker = kernel
    
    vmin = 0
    vmax = np.percentile(data[~np.isnan(data)], 99.5)
    
    plt.figure(1)
    plt.clf()
    plt.imshow(data, origin = 'lower', vmin = vmin, vmax = vmax)
    
    plt.figure(2)
    plt.clf()
    plt.imshow(kernel, origin = 'lower', norm=LogNorm(vmin=1e-12, vmax=1e-6))
    
    plt.figure(3)
    plt.clf()
    plt.imshow(pad_data, origin = 'lower', vmin = vmin, vmax = vmax)
    
    # PIXSCALE is given in arcsec
    pixel_scale_kernel = abs(kernel_header['PIXELSCL'])
    # CDELT is given in degrees
    pixel_scale = abs(data_header['CDELT1'])*3600.0
    
    print('pad_data', np.shape(pad_data))
    
    new_data = regrid(pad_data, pixel_scale, pixel_scale_kernel, method = 'spline', center = False)
    print('new_data', np.shape(new_data))

    conv_data = convolve_fft(new_data, kernel, psf_pad = True, allow_huge=True, nan_treatment = 'fill', fill_value = 0.0)
    print('conv_data', np.shape(conv_data))
    
    plt.figure(4)
    plt.imshow(conv_data, origin = 'lower', vmin = vmin, vmax = vmax)
    
    size_data = len(conv_data)
    size_new = int(round( float(size_data) * pixel_scale_kernel / pixel_scale))
    
    print(size_new, pixel_scale_kernel, pixel_scale)
    
    conv_data = regrid(conv_data, pixel_scale_kernel, pixel_scale)
    print('conv_data', np.shape(conv_data))
    
    # final_data = conv_data
    
    final_data = conv_data[xstart:xend, ystart:yend]
    print('final_data', np.shape(final_data))

    # add back in the nans that were removed before
    final_data[nan_map] = np.nan
    
    if convolve_err:
        final_data = np.sqrt(conv_data)
    else:
        final_data = final_data
    
    plt.figure(5)
    plt.clf()
    plt.imshow(final_data, origin = 'lower', vmin = vmin, vmax = vmax)
    
    
    hdulist_data[0].data = final_data
    
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/convolve_rot/'
    
    if convolve_err:
        filename = f'{filt}_convolve_to_F1500W_ERR'
    else:
        filename = f'{filt}_convolve_to_F1500W_rot'
    # fits.writeto(filepath + filename + '.fits', final_data, data_header, overwrite = True)
    
    # checking flux conservation
    flux1 = np.nansum(data)
    flux2 = np.nansum(final_data)
    
    print(f'{filt} Flux Conservation')
    print(flux1, flux2)
    print(100*((flux1-flux2)/flux1))
    
    print(np.shape(final_data))
    
    
    # try a smaller area
    x1 = 3000
    x2 = 6000
    y1 = 4000
    y2 = 6000
    
    # checking flux conservation
    flux1 = np.nansum(data[y1:y2, x1:x2])
    flux2 = np.nansum(final_data[y1:y2, x1:x2])
    
    print(f'{filt} Flux Conservation')
    print(flux1, flux2)
    print(100*((flux1-flux2)/flux1))
