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

savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/align_north/'
nircam_dir = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/nircam_nov23/'
miri_dir = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/miri_nov23/'

savepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_mcf/'

filts = ['F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M', 'F560W', 'F770W', 'F1000W', 'F1130W', 'F1500W']

# reproject everything to F200W
F200W_file = nircam_dir + '/F200W/jw02391007001_nircam_f200w_asn_modb_bkgrsub_i2d.fits'

out_head = fits.open(F200W_file)['SCI'].header
out_data = fits.open(F200W_file)['SCI'].header
out_hdu = fits.open(F200W_file)['SCI']
s = np.shape(out_data)

frame = ICRS()

north_wcs, north_shape = find_optimal_celestial_wcs([out_hdu], frame=frame)
north_head = north_wcs.to_header()


for filt in filts:
    if filt in miri_filts:
        path = miri_dir + f'{filt}/stage3/'
        infile = glob.glob(path + '*_skysub_i2d.fits')[0]
        
    elif filt in nircam_filts:
        path = nircam_dir + f'{filt}/'
        infile = glob.glob(path + '*{:s}_asn_modb_bkgrsub_i2d.fits'.format(filt.lower()))[0]
    
    in_hdu = fits.open(infile)['SCI']
    
    reproj, foot = reproject_exact(in_hdu, output_projection = north_wcs, shape_out = north_shape)
    
    
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_mcf/'
    filename = f'SextansA_jw2391_{filt}_reproject_mcf_north'
    fits.writeto(filepath + filename + '.fits', reproj, north_head, overwrite = True)
    

# # sum the PAH filters
# F770W_hdu = fits.open(savepath + 'SextansA_jw2391_F770W_reproject_mcf.fits')
# F1130W_hdu = fits.open(savepath + 'SextansA_jw2391_F1130W_reproject_mcf.fits')


# pah_sum = F770W_hdu[0].data + F1130W_hdu[0].data

# filename = 'SextansA_jw2391_F770WF1130WSUM_reproject_mcf'
# fits.writeto(savepath + filename + '.fits', pah_sum, F770W_hdu[0].header, overwrite = True)