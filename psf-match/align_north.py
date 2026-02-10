#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Uses the reproject package to align each image to north before PSF matching. 
Required because PSFs are rotated to match the PA of each image

@author: etarantino
"""

from astropy.wcs import WCS
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_exact
from astropy.coordinates import ICRS
from astropy.io import fits, ascii
import glob
import numpy as np 


miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'
nircam_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/nircam/'
miri_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/miri/'

# filts = ['F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W']
filts = ['F115W', 'F150W', 'F200W']
err = True

# code for the error arrays
for filt in filts:
    if filt in miri_filts:
        infile = glob.glob(miri_dir + f'*{filt}_skysub_i2d.fits')[0]
        
    elif filt in nircam_filts:
        infile = glob.glob(nircam_dir + '*{:s}_*i2d.fits'.format(filt.lower()))[0]
    
    sci_hdu = fits.open(infile)['SCI']
    err_data = fits.open(infile)['ERR'].data
    frame = ICRS()

    wcs, shape = find_optimal_celestial_wcs([sci_hdu], frame=frame)
    
    data, _ = reproject_exact((err_data, sci_hdu.header), wcs, shape_out = shape)
    
    header = wcs.to_header()
    
    savename = f'jw2391_{filt}_err_align_north.fits'
    
    fits.writeto(savedir + savename, data, header, overwrite = True )

### code for the regular images

for filt in filts:
    if filt in miri_filts:
        infile = glob.glob(miri_dir + f'*{filt}_skysub_i2d.fits')[0]
        
    elif filt in nircam_filts:
        infile = glob.glob(nircam_dir + '*{:s}_*_i2d.fits'.format(filt.lower()))[0]
        # infile = glob.glob(nircam_dir + '*{:s}_*jhat_i2d.fits'.format(filt.lower()))[0]
    
    hdu = fits.open(infile)['SCI']
    frame = ICRS()

    wcs, shape = find_optimal_celestial_wcs([hdu], frame=frame)
    
    data, _ = reproject_exact(hdu, wcs, shape_out = shape)
    
    header = wcs.to_header()
    
    savename = f'jw2391_{filt}_align_north.fits'
    
    # checking flux conservation
    flux1 = np.nansum(hdu.data)
    flux2 = np.nansum(data)
    
    print(f'{filt} Flux Conservation')
    print(flux1, flux2)
    print(100*((flux1-flux2)/flux1))
    
    fits.writeto(savedir + savename, data, header, overwrite = True )
    
    