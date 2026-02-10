#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates dendrograms of the continuum subtracted PAH emission

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from astrodendro import Dendrogram, pp_catalog
from astrodendro.analysis import PPStatistic
from astropy import units as u

# try the 7.7 first as an intermediate test
filt_mid = 'F770W'
datapath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject/'
dataname = f'F770W_reproject_to_F1500W.fits'

pah_hdr = fits.open(datapath + dataname)[0]
pah_data = pah_hdr.data
pah_head = pah_hdr.header

pah_sigma = 0.004589

x1 = 800
x2 = 940
y1 = 650
y2 = 782

pah_data = pah_data[y1:y2, x1:x2]

d = Dendrogram.compute(pah_data, min_value = pah_sigma, min_delta = pah_sigma, min_npix = 4)

v = d.viewer()
v.show()

# metadata for statistics
metadata = {}
metadata['data_unit'] = 1e6 * u.Jy / u.sr
metadata['spatial_scale'] =  3.0555555555555E-05 * u.deg
metadata['beam_major'] =  0.488 * u.arcsec # FWHM
metadata['beam_minor'] =  0.488 * u.arcsec # FWHM

# calculate statistics 
cat = pp_catalog(d, metadata)

# # save as a region file so we can remove the clumps that are not real PAH clumps 
# clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'

# # create region file
# reg_name = f'{filt_mid}_dendro_cat_ds9'
# reg_file = open(clump_path + reg_name + '.reg', 'w')

# reg_file.write('# Region file format: DS9 version 4.1\n')
# reg_file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
# reg_file.write('image\n')


# for i in range(len(cat)):
#     reg_file.write('point({},{}) # point=x text={{{}}} color=magenta \n'.format(cat['x_cen'][i], cat['y_cen'][i], cat['_idx'][i]))
    
# reg_file.close()

# get the clumps that are only in the previous catalog with detections in all PAH features
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_file = 'clump_catalog_with_dendro.txt'

clumps = ascii.read(clump_path + clump_file)
filt_idx = f'{filt_mid}_dendro'

new_cat = cat[clumps[filt_idx]]

savefile = f'{filt_mid}_dendro_clump_cat.txt'
# ascii.write(new_cat, clump_path + savefile)

plt.figure(2)
plt.clf()
plt.scatter(new_cat['radius'].to(u.arcsec), new_cat['flux'], alpha = 0.7, edgecolor = 'k', s = 50)
plt.xlabel('Radius (arcseconds)')
plt.ylabel('Clump Flux (Jy)')


plt.figure(3)
plt.clf()
plt.scatter(clumps['Clump'], new_cat['radius'].to(u.arcsec), alpha = 0.7, edgecolor = 'k', s = 50)
plt.xlabel('Clump Number')
plt.ylabel('Clump Radius (arcseconds)')

plt.figure(4)
plt.clf()
plt.scatter(clumps['Clump'], new_cat['flux'], alpha = 0.7, edgecolor = 'k', s = 50)
plt.xlabel('Clump Number')
plt.ylabel('Clump Flux (Jy)')

plt.figure(5)
plt.clf()
plt.scatter(clumps['Clump'], new_cat['major_sigma'].to(u.arcsec), alpha = 0.7, edgecolor = 'k', s = 50)
plt.xlabel('Clump Number')
plt.ylabel('Major Axis (arcseconds)')


