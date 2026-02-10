#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 21:40:44 2023

@author: etarantino
"""
import aplpy
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits

plt.ion()

con_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images/'
reproj_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject/'

pah_3_file = 'SextansA_F335M_sub_match_pah.fits'
pah_7_file = 'SextansA_F770W_sub_match_pah.fits'
pah_11_file = 'SextansA_F1130W_sub_match_pah.fits'

fig = plt.figure(2, figsize = (8,9))
plt.clf()

filt = 'F770W'

f = aplpy.FITSFigure(con_path + pah_7_file, figure=fig, north = True)

vmin = -0.01
vmax = 0.1


# cen_x = 152.7772917
# cen_y = -4.7013056

# f.recenter(cen_x, cen_y, width=0.008,height=0.008)


fig.canvas.draw()

f.add_scalebar(1 * u.arcsec, corner = 'bottom left', color = 'white')
f.scalebar.set_label("1'' = 7 pc")

reg_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/dolphot/reg/sw_lw_v2/'
reg_name = 'YSO_line_cut'
f.show_regions(reg_path + reg_name + '.reg')

# load the HI in contours
ancillary = '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
HI_file = 'DDO75_R_X0_P_R.FITS'

# hdu = fits.open(ancillary + HI_file)
# data = hdu[0].data
# data = data[0,0,:,:]
# header = hdu[0].header

# header['NAXIS'] = 2
# del header['NAXIS3']
# del header['CTYPE3']
# del header['CRVAL3']
# del header['CDELT3']
# del header['CRPIX3']
# del header['CROTA3']

# del header['NAXIS4']
# del header['CTYPE4']
# del header['CRVAL4']
# del header['CDELT4']
# del header['CRPIX4']
# del header['CROTA4']

# fits.writeto(ancillary + 'DDO75_R_X0_P_R_FLATTEN.FITS', data, header = header)

HI_file = 'DDO75_R_X0_P_R_FLATTEN.FITS'

levels = [100, 150, 200, 230]
f.show_contour(ancillary + HI_file, levels = levels, colors = 'b')
# f.add_beam()



f.show_grayscale()

plt.ylim(585, 842)
plt.xlim(742, 978)

plt.savefig('HI_YSO_map.pdf')