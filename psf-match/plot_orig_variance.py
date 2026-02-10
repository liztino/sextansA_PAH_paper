#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Looks at each variance image for the original data  

@author: etarantino
"""

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve_fft, Gaussian1DKernel, Box1DKernel
import scipy.signal
import numpy as n
import scipy.interpolate
import scipy.ndimage
import glob
from astropy.stats import sigma_clip
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5 
import astropy.units as u

plt.ion()

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'
nircam_dir = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/nircam_nov23/'
miri_dir = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/miri_nov23/'

# filts = ['F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W']
filts = ['F300M']

err = True

# code for the error arrays
for filt in filts:
    # if filt in miri_filts:
    #     path = miri_dir + f'{filt}/stage3/'
    #     infile = glob.glob(path + '*_skysub_i2d.fits')[0]
        
    # elif filt in nircam_filts:
    #     path = nircam_dir + f'{filt}/'
    #     infile = glob.glob(path + '*{:s}_asn_modb_bkgrsub_i2d.fits'.format(filt.lower()))[0]
    
    inpath = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/nircam_nov23/F300M/'
    infile = inpath + 'jw02391007001_03101_00003_nrcblong_bkgrsub_a3001_crf.fits'
    
    hdu = fits.open(infile)
    hdu.info()
    
    data = hdu['SCI'].data
    err =  hdu['ERR'].data
    poisson = hdu['VAR_POISSON'].data
    rnoise = hdu['VAR_RNOISE'].data
    flat = hdu['VAR_FLAT'].data
    
    plt.figure(2)
    plt.clf()
    fig, axs = plt.subplots(2,3, num = 2)
    axs = axs.flatten()
    
    vmax = np.nanpercentile(data, 95)
    vmin = np.nanpercentile(data, 1)
    im = axs[0].imshow(data, origin = 'lower', vmin = vmin, vmax = vmax)
    axs[0].set_title('DATA')
    plt.colorbar(im, ax = axs[0], fraction=0.046, pad=0.04)
    
    vmax = np.nanpercentile(err, 95)
    vmin = np.nanpercentile(err, 1)
    im = axs[1].imshow(err, origin = 'lower', vmin = vmin, vmax = vmax)
    axs[1].set_title('ERR')
    plt.colorbar(im, ax = axs[1], fraction=0.046, pad=0.04)
    
    vmax = np.nanpercentile(poisson, 95)
    vmin = np.nanpercentile(poisson, 1)
    im = axs[2].imshow(poisson, origin = 'lower', vmin = vmin, vmax = vmax)
    axs[2].set_title('POISSON')
    plt.colorbar(im, ax = axs[2], fraction=0.046, pad=0.04)
    
        
    vmax = np.nanpercentile(rnoise, 95)
    vmin = np.nanpercentile(rnoise, 1)
    im = axs[3].imshow(rnoise, origin = 'lower', vmin = vmin, vmax = vmax)
    axs[3].set_title('RNOISE')
    plt.colorbar(im, ax = axs[3], fraction=0.046, pad=0.04)
    
        
    vmax = np.nanpercentile(flat, 95)
    vmin = np.nanpercentile(flat, 1)
    im = axs[4].imshow(flat, origin = 'lower', vmin = vmin, vmax = vmax)
    axs[4].set_title('FLAT')
    plt.colorbar(im, ax = axs[4], fraction=0.046, pad=0.04)
    
    fig.delaxes(axs[5])
    
    
    plt.tight_layout()
    
    savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/plots/'
    savename = f'{filt}_variance_tile3.pdf'
    
    plt.savefig(savepath + savename)

    
    
    
    
    
    
    
    