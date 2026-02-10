#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Create lupton RGB images for both Sextans A and IC 1613

@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
import numpy.ma as ma
import astropy.visualization.wcsaxes as vis
import astropy.units as u 
import matplotlib
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
import seaborn as sns
from regions import Regions

from astropy.visualization import make_lupton_rgb


plt.ion()

from get_pivot_wave import get_pivot_wave


savepath = '/Users/etarantino/Documents/proposals/ALMA_C11/SexA_IC1613/figures/'

### First the Sextans A data ####


def load_filter(filt):
    
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_mcf/'
    
    # load the middle filter we will be continuum subtracting from
    filename = filepath + f'SextansA_jw2391_{filt}_reproject_mcf_north'
    hdu = fits.open(filename + '.fits')
    header = hdu[0].header
    data = hdu[0].data
    pivot = get_pivot_wave(filt)
    
    filt_dict = {'name': filt, 'data': data, 'header': header, 'wave': pivot}
    
    return filt_dict

F560W = load_filter('F560W')
F770W = load_filter('F770W')
F1000W = load_filter('F1000W')

min_vals = [0,0,0]
image = make_lupton_rgb(F1000W['data'], F770W['data'], F560W['data'], stretch=0.1, minimum = min_vals)

w = WCS(F560W['header'])
fig1 = plt.figure(1)
plt.clf()
ax = fig1.add_subplot(projection = w)

im = ax.imshow(image, origin = 'lower')
plt.xlabel('RA')
plt.ylabel('DEC')

x1 = 2881
x2 = 6593
y1 = 3366
y2 = 6643

plt.xlim(x1, x2)
plt.ylim(y1, y2)

# HI data
HI_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
HI_file = 'DDO75_R_X0_P_R.FITS'

HI_hdu = fits.open(HI_path + HI_file)[0]
HI_data = HI_hdu.data
HI_header = HI_hdu.header
HI_w = WCS(HI_header, naxis = [1,2])

# contours for HI
levels = [50, 100, 150, 200, 250]
ax.contour(HI_data[0,0,:,:], levels = levels, colors = 'gainsboro', transform = ax.get_transform(HI_w), linewidths =0.6, alpha = 0.7)


HI_bmaj = 2.0980E-03* u.degree
HI_bmin = 1.7921E-03 * u.degree
HI_bpa =  13.70  * u.degree

arcsecond = 1
d_ic1613 = 1.3
res_pc = (arcsecond/206265) * (d_ic1613 * 1e6)
print(res_pc)

vis.add_beam(ax = ax, header = None, major = HI_bmaj,  minor = HI_bmin, angle = HI_bpa, color = 'gainsboro', frame = True, label = 'HI beam')
vis.add_scalebar(ax = ax, length = 10 * u.arcsec, corner = 'bottom right', label = '10" = 70 pc', color = 'white')

CO_bmaj = 0.3 * u.arcsecond
CO_bmin = 0.3 * u.arcsecond
CO_bpa =  0 * u.degree
# vis.add_beam(ax = ax, header = None, major =CO_bmaj,  minor = CO_bmin, angle = CO_bpa, color = 'cyan', frame = True)

# reg_path = '/Users/etarantino/Documents/proposals/ALMA_C11/SexA_IC1613/reg/'
# reg_name = 'SextansA_ALMA_C11.reg'
# reg = Regions.read(reg_path + reg_name, format='ds9')[0]

# pix_reg = reg.to_pixel(w)

# pix_reg.plot()

# add the PAH region

# reg_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/reg/'
# reg_name = 'pah_box.reg'
# reg = Regions.read(reg_path + reg_name, format='ds9')[0]

# pix_reg = reg.to_pixel(w)

# pix_reg.plot()

# save 
savename = 'SextansA_luptonRGB_Present.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)


# ### Now IC 1613 data ####


def load_filter(filt):
    
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613/aug23_reduction/reproject_north/'
    
    # load the middle filter we will be continuum subtracting from
    filename = filepath + f'IC1613_jw2391_{filt}_reproject_mcf_north'
    hdu = fits.open(filename + '.fits')
    header = hdu[0].header
    data = hdu[0].data
    pivot = get_pivot_wave(filt)
    
    filt_dict = {'name': filt, 'data': data, 'header': header, 'wave': pivot}
    
    return filt_dict

F560W = load_filter('F560W')
F770W = load_filter('F770W')
F1000W = load_filter('F1000W')

min_vals = [0,0,0]
image = make_lupton_rgb(F1000W['data'], F770W['data'], F560W['data'], stretch=0.25, minimum = min_vals)

w = WCS(F560W['header'])
fig1 = plt.figure(2)
plt.clf()
ax = fig1.add_subplot(projection = w)

im = ax.imshow(image, origin = 'lower')
plt.xlabel('RA')
plt.ylabel('DEC')

x1 = 156
x2 = 1352
y1 = 148
y2 = 1189

plt.xlim(x1, x2)
plt.ylim(y1, y2)

# HI data
HI_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613_ANCILLARY/'
HI_file = 'IC1613_R_X0_P_R.FITS'

HI_hdu = fits.open(HI_path + HI_file)[0]
HI_data = HI_hdu.data
HI_header = HI_hdu.header
HI_w = WCS(HI_header, naxis = [1,2])

# contours for HI
levels = [50, 110, 140]
ax.contour(HI_data[0,0, :, :], levels = levels, colors = 'gainsboro', transform = ax.get_transform(HI_w), linewidths =0.6, alpha = 0.7)

HI_bmaj = 2.1363E-03 * u.degree
HI_bmin = 1.7989E-03 * u.degree
HI_bpa =  6.69 * u.degree

arcsecond = 1
d_ic1613 = 0.720
res_pc = (arcsecond/206265) * (d_ic1613 * 1e6)

vis.add_beam(ax = ax, header = None, major = HI_bmaj,  minor = HI_bmin, angle = HI_bpa, color = 'gainsboro', frame = True, label = 'HI beam')
vis.add_scalebar(ax = ax, length = 10 * u.arcsec, corner = 'bottom right', label = '10" = 35 pc', color = 'white')

CO_bmaj = 0.3 * u.arcsecond
CO_bmin = 0.3 * u.arcsecond
CO_bpa =  0 * u.degree
# vis.add_beam(ax = ax, header = None, major =CO_bmaj,  minor = CO_bmin, angle = CO_bpa, color = 'cyan', frame = True)

# reg_path = '/Users/etarantino/Documents/proposals/ALMA_C11/SexA_IC1613/reg/'
# reg_name = 'IC1613_ALMA_C11.reg'
# reg = Regions.read(reg_path + reg_name, format='ds9')[0]

# pix_reg = reg.to_pixel(w)

# pix_reg.plot()

savename = 'IC1613_luptonRGB_present.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)




