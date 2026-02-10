#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Reprojects the PSF-matched/convolved filter images to the WCS of the F1500W filter

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
import glob


convolve_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/convolve_rot/'
# filt_list = ['F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W' ]
filt_list = ['F300M', 'F335M', 'F360M']

F1500W_file = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/jw2391_F1500W_align_north.fits'

savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'
nircam_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/nircam/'
miri_dir = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/miri/'

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

err = False

if err == True:
        
    for filt in filt_list:
    
        infile = f'{filt}_convolve_to_F1500W_rot_ERR_kernelsq_norm'
        in_hdu = fits.open(convolve_dir + infile + '.fits')
        
        out_head = fits.open(F1500W_file)[0].header
        out_data = fits.open(F1500W_file)[0].header
        s = np.shape(out_data)
        
        # square the uncertainties to get a variance
        data = in_hdu[0].data
        data = data * data
        in_hdu[0].data = data
        
        reproj, foot = reproject.reproject_exact(in_hdu, output_projection = out_head, shape_out = s)
        
        # squareroot to get back to a standard deviation 
        reproj = np.sqrt(reproj)
        
        filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_orig/'
        filename = f'{filt}_reproject_to_F1500W_rot_ERR_kernelsq_norm'
        fits.writeto(filepath + filename + '.fits', reproj, out_head, overwrite = True)
    
else:
    
    for filt in filt_list:
        
        if filt in miri_filts:
            infile = glob.glob(miri_dir + f'*{filt}_skysub_i2d.fits')[0]
            
        elif filt in nircam_filts:
            infile = glob.glob(nircam_dir + '*{:s}_*i2d.fits'.format(filt.lower()))[0]
        
        sci_hdu = fits.open(infile)['SCI']
        err_data = fits.open(infile)['ERR'].data
    
        
        out_head = fits.open(F1500W_file)[0].header
        out_data = fits.open(F1500W_file)[0].header
        s = np.shape(out_data)
        
        reproj, foot = reproject.reproject_exact(sci_hdu, output_projection = out_head, shape_out = s)
        
        filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_orig/'
        filename = f'{filt}_reproject_orig_to_F1500W'
        fits.writeto(filepath + filename + '.fits', reproj, out_head, overwrite = True)
        
        
    
# just make an F1500W file
infile = 'jw2391_F1500W_align_north'
path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'

in_hdu = fits.open(path + infile + '.fits')

data = in_hdu[0].data
header = in_hdu[0].header
# err = in_hdu['ERR'].data

filename = f'F1500W_reproject_to_F1500W'
filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
fits.writeto(filepath + filename + '.fits', data, header, overwrite = True)

# filename = 'F1500W_reproject_to_F1500W_ERR'
# filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject/'
# fits.writeto(filepath + filename + '.fits', err, header, overwrite = True)