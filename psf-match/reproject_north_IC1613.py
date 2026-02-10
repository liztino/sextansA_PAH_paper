#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 12:16:58 2024

@author: etarantino
"""

from astropy.wcs import WCS
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_exact
from astropy.coordinates import ICRS
from astropy.io import fits, ascii
import numpy as np
import glob

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

nircam_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613/aug23_reduction/NIRCam/'
miri_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613/aug23_reduction/MIRI/'

savedir =  '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613/aug23_reduction/reproject_north/'

filts = ['F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M', 'F560W', 'F770W', 'F1000W', 'F1130W', 'F1500W']

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

for filt in filts:
    if filt in miri_filts:
        infile = miri_dir + f'miri_{filt}_stage3_asn_skysub_prop_i2d.fits'
        
    elif filt in nircam_filts:
        infile =  nircam_dir + f'jw02391-ic1613_nircam_{filt}_i2d.fits'
    
    hdu = fits.open(infile)['SCI']
    frame = ICRS()

    wcs, shape = find_optimal_celestial_wcs([hdu], frame=frame)
    
    data, _ = reproject_exact(hdu, wcs, shape_out = shape)
    
    header = wcs.to_header()
    
    savename = f'jw2391_{filt}_align_north.fits'
    
    fits.writeto(savedir + savename, data, header, overwrite = True )