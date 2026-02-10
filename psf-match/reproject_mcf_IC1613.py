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

# just do some of the MIRI filters for now 
filts = ['F300M', 'F335M', 'F360M']

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

    
    
# reproject everything to F300M
miri_file = nircam_dir + 'jw02391-ic1613_nircam_F300M_i2d.fits'

out_head = fits.open(miri_file)['SCI'].header
out_data = fits.open(miri_file)['SCI'].header
out_hdu = fits.open(miri_file)['SCI']
s = np.shape(out_data)

frame = ICRS()

north_wcs, north_shape = find_optimal_celestial_wcs([out_hdu], frame=frame)
north_head = north_wcs.to_header()


for filt in filts:
    if filt in miri_filts:
        infile = miri_dir + f'miri_{filt}_stage3_asn_skysub_prop_i2d.fits'
        
    elif filt in nircam_filts:
        infile =  nircam_dir + f'jw02391-ic1613_nircam_{filt}_i2d.fits'
    
    in_hdu = fits.open(infile)['SCI']
    
    reproj, foot = reproject_exact(in_hdu, output_projection = north_wcs, shape_out = north_shape)
    
    
    filename = f'IC1613_jw2391_{filt}_reproject_mcf_north'
    fits.writeto(savedir + filename + '.fits', reproj, north_head, overwrite = True)