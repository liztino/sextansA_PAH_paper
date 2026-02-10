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

ds9_save = False

# try the 7.7 first as an intermediate test
filt_mid = 'F770W'
k = 4.33
datapath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
dataname = 'SextansA_{:s}_k_method_pah_k_{:3.2f}.fits'.format(filt_mid, k)

pah_hdr = fits.open(datapath + dataname)[0]
pah_data = pah_hdr.data
pah_head = pah_hdr.header

# # quick ratio map, don't mind me
# conname = f'SextansA_{filt_mid}_k_method_con.fits'
# con_hdr = fits.open(datapath + conname)[0]
# con_data = con_hdr.data

# ratio = pah_data/con_data

# savename = f'SextansA_{filt_mid}_k_method_ratio.fits'
# fits.writeto(datapath + savename, ratio, header = pah_head, overwrite = True)

if filt_mid == 'F335M':
	pah_sigma =  0.0036
elif filt_mid == 'F770W':
	pah_sigma =  0.0066
elif filt_mid == 'F1130W':
	pah_sigma = 0.0138
    
# # for the 3.3 micron feature, cut out some of the area
# x1 = 600
# x2 = 1100
# y1 = 300
# y2 = 930
# pah_data = pah_data[y1:y2, x1:x2]

d = Dendrogram.compute(pah_data, min_value = pah_sigma, min_delta = 3*pah_sigma, min_npix = 4)

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


# if ds9_save:
#  	# save as a region file so we can remove the clumps that are not real PAH clumps 
#  	clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/dendro_clump_cats/'

#  	# create region file
#  	reg_name = '{:s}_k_{:3.2f}_dendro_cat_ds9_v2'.format(filt_mid, k)
#  	reg_file = open(clump_path + reg_name + '.reg', 'w')

#  	reg_file.write('# Region file format: DS9 version 4.1\n')
#  	reg_file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
#  	reg_file.write('image\n')


#  	for i in range(len(cat)):
#  	    reg_file.write('point({},{}) # point=x text={{{}}} color=green \n'.format(cat['x_cen'][i] , cat['y_cen'][i] , cat['_idx'][i]))
 	    
#  	reg_file.close()

# else:

#  	# get the clumps that are only in the previous catalog with detections in all PAH features
#  	clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
#  	clump_file = 'master_clump_cat.csv'

#  	clumps = ascii.read(clump_path + clump_file)
#  	filt_idx = '{:s}_k_{:3.2f}'.format(filt_mid, k)
#  	new_cat = cat[clumps[filt_idx]]
    
#  	savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/clump_prop_cats/'
#  	savefile = '{:s}_k_{:3.2f}_dendro_clump_props_v2.txt'.format(filt_mid, k)
#  	ascii.write(new_cat, savepath + savefile)

#  	plt.figure(2)
#  	plt.clf()
#  	plt.scatter(new_cat['radius'].to(u.arcsec), new_cat['flux'], alpha = 0.7, edgecolor = 'k', s = 50)
#  	plt.xlabel('Radius (arcseconds)')
#  	plt.ylabel('Clump Flux (Jy)')


#  	plt.figure(3)
#  	plt.clf()
#  	plt.scatter(clumps['Clump'], new_cat['radius'].to(u.arcsec), alpha = 0.7, edgecolor = 'k', s = 50)
#  	plt.xlabel('Clump Number')
#  	plt.ylabel('Clump Radius (arcseconds)')

#  	plt.figure(4)
#  	plt.clf()
#  	plt.scatter(clumps['Clump'], new_cat['flux'], alpha = 0.7, edgecolor = 'k', s = 50)
#  	plt.xlabel('Clump Number')
#  	plt.ylabel('Clump Flux (Jy)')

#  	plt.figure(5)
#  	plt.clf()
#  	plt.scatter(clumps['Clump'], new_cat['major_sigma'].to(u.arcsec), alpha = 0.7, edgecolor = 'k', s = 50)
#  	plt.xlabel('Clump Number')
#  	plt.ylabel('Major Axis (arcseconds)')

# finding the mask for each clump 
# 	# get the clumps that are only in the previous catalog with detections in all PAH features
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_file = 'template_clump_F770W.csv'

clumps = ascii.read(clump_path + clump_file, format = 'csv')
filt_idx = '{:s}_k_{:3.2f}'.format(filt_mid, k)
clump_list = clumps[filt_idx]

# loop through the clumps to create a mask for each
for i in range(len(clump_list)):
    mask = np.where(d.index_map == clump_list[i])
    
    # add the features for the trunks/branches to ensure the leaves are included
    if clump_list[i] == 28:
        union1 = np.logical_or((d.index_map == 28), (d.index_map == 29))
        union2 = np.logical_or(union1, (d.index_map == 31))
        
        mask = np.where(union2 == True)
        plt.figure(2)
        plt.clf()
        plt.imshow(union2, origin = 'lower')
        
    if clump_list[i] == 35:
        union1 = np.logical_or((d.index_map == 35), (d.index_map == 36))
        union2 = np.logical_or(union1, (d.index_map == 37))
        
        mask = np.where(union2 == True)
        plt.figure(3)
        plt.clf()
        plt.imshow(union2, origin = 'lower')
        
    save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/clump_masks/'
    save_name = '{:s}_clump{:d}_mask.txt'.format(filt_idx, clumps['Clump'][i])
    ascii.write(mask, save_path + save_name, overwrite = True)
    
# grab and save the catalog file so we can use it in the clump properties table

# get the clumps that are only in the previous catalog with detections in all PAH features
new_cat = cat[clumps[filt_idx]]
new_cat.add_column(clumps['Clump'], name = 'clump_num', index = 0)



# add a column describing if object is a leaf or a trunk
dendro_struct = ['leaf'] * len(new_cat)
dendro_struct[-1] = 'branch'
dendro_struct[-2] = 'branch'

new_cat.add_column(dendro_struct, name = 'dendro_struct')

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
savefile = 'template_clump_F770W_props.txt'
ascii.write(new_cat, savepath + savefile)

    
