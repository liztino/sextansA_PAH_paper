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
from fwhm import fwhm 
import copy

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
filt_list = ['F300M', 'F335M', 'F360M']

for filt in filt_list:

    miri_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'
    nircam_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'
    
    miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
    nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
    
    if filt in miri_filts:
        path = miri_path
        # image = glob.glob(miri_path + '*{}_*skysub_prop_i2d.fits'.format(filt))[0]
        image = glob.glob(miri_path + '*{}_err_align_north.fits'.format(filt))[0]
    else:
        path = nircam_path
        # image = glob.glob(nircam_path + '*{}_*modb_bkgrsub_i2d.fits'.format(filt.lower()))[0]
        image = glob.glob(miri_path + '*{}_err_align_north.fits'.format(filt))[0]
    
    
    filepath_ker = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/kernels/rot_ker/'
    filename_ker = glob.glob(filepath_ker + '*{}*'.format(filt))[0]
    
    hdulist_data = fits.open(image)
    
    data = hdulist_data[0].data
    data_header = hdulist_data[0].header
    orig_data = copy.deepcopy(data)
    data = data * data
    
    nan_map = np.isnan(data)
    data[nan_map] = 0.0
        
    hdulist_ker = fits.open(filename_ker)
    kernel = copy.deepcopy(hdulist_ker[0].data)
    kernel_header = hdulist_ker[0].header
    
    # Klein 2021 says to square the kernel
    kernel = kernel * kernel
    
    # normalize
    kernel = kernel/np.sum(kernel)
    
    #
    
    data_shape = np.shape(data)
    ker_shape = np.shape(kernel)
    
    print(data_shape, ker_shape)
    
    pad = 200
    pad_data, xstart, xend, ystart, yend = pad_zeros(pad, data)
    
    pad_ker = kernel
    
    vmin = 0
    vmax = np.percentile(orig_data[~np.isnan(data)], 99.5)
    
    plt.figure(1)
    plt.clf()
    plt.imshow(orig_data, origin = 'lower', vmin = vmin, vmax = vmax)
    
    plt.figure(2)
    plt.clf()
    plt.imshow(kernel, origin = 'lower', vmin = 0, vmax = 0.01)
    
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
    
    conv_data = convolve_fft(new_data, kernel, psf_pad = True, allow_huge=True, nan_treatment = 'fill', fill_value = 0.0, normalize_kernel = False)
    print('conv_data', np.shape(conv_data))
    
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
        
    plt.figure(6)
    plt.clf()
    plt.imshow(final_data, origin = 'lower', vmin = vmin, vmax = vmax)
            
    
    # add the factor for correlated data from Klein+ 2021 (research note)
    psf1_file = f'{filt}_webbpsfv121_rot'
    psf2_file = 'F1500W_webbpsfv121_rot'
    
    psf_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/psfs/'

    # open psf files
    psf1_hdu = fits.open(psf_path + psf1_file + '.fits')
    psf1 = psf1_hdu[0].data
    psf1_head = psf1_hdu[0].header

    psf2_hdu = fits.open(psf_path + psf2_file + '.fits')
    psf2 = psf2_hdu[0].data
    psf2_head = psf2_hdu[0].header
    
    width1 = fwhm(psf1, psf1_head)
    width2 = fwhm(psf2, psf2_head)
    
    # factor for correlated data
    factor = ((2 * np.sqrt(np.pi) * width1['xarcsec'] * width2['xarcsec'])/(np.sqrt(width1['xarcsec']**2 + width2['xarcsec']**2)))
    
    final_data = factor * final_data
    # final_data = final_data
    
    final_data = np.sqrt(final_data)

    
    plt.figure(5)
    plt.clf()
    plt.imshow(final_data, origin = 'lower', vmin = vmin, vmax = vmax)
    
    
    hdulist_data[0].data = final_data
    
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/convolve_rot/'
    
    filename = f'{filt}_convolve_to_F1500W_rot_ERR_kernelsq_norm'
    fits.writeto(filepath + filename + '.fits', final_data, data_header, overwrite = True)
    
    # checking flux conservation
    flux1 = np.nansum(orig_data)
    flux2 = np.nansum(final_data)
    
    print(f'{filt} Flux Conservation')
    print(flux1, flux2)
    print(100*((flux1-flux2)/flux1))
    
    print(np.shape(final_data))