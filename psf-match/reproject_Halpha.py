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

filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
# filename = filepath + 'MUSE_SexA_pipeline_Halpha'
filename = filepath + 'Sextans_A_I_Ha_d2009'


F1500W_file = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/jw2391_F1500W_align_north.fits'

in_hdu = fits.open(filename + '.fits')
        
out_head = fits.open(F1500W_file)[0].header
out_data = fits.open(F1500W_file)[0].header
s = np.shape(out_data)

reproj, foot = reproject.reproject_exact(in_hdu, output_projection = out_head, shape_out = s)

filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
filename = 'Sextans_A_I_Ha_d2009_proj_F1500W'
fits.writeto(filepath + filename + '.fits', reproj, out_head, overwrite = True)
        
        
    
# # just make an F1500W file
# infile = 'jw2391_F1500W_align_north'
# path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'

# in_hdu = fits.open(path + infile + '.fits')

# data = in_hdu[0].data
# header = in_hdu[0].header
# # err = in_hdu['ERR'].data

# filename = f'F1500W_reproject_to_F1500W'
# filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
# fits.writeto(filepath + filename + '.fits', data, header, overwrite = True)

# filename = 'F1500W_reproject_to_F1500W_ERR'
# filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject/'
# fits.writeto(filepath + filename + '.fits', err, header, overwrite = True)