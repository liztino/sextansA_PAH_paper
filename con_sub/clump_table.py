#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates table of clump properties-- this is the ascii table that we will use to build the Latex table

TODO: maybe add a third file that will calculate the value for each filter?

@author: etarantino
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import fits, ascii
from astropy.table import Table, hstack
from astropy.wcs import WCS

clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'

# get the dendrogram file for the 7.7 clumps that made the mask 
template_file = 'template_clump_F770W_props.txt'
template = ascii.read(clump_path + template_file)

# get the table that used the masks and calculate all the PAH stuff
pah_flux_file = 'dendro_mask_clump_flux_final.txt'
pah_flux = ascii.read(clump_path + pah_flux_file)

# get the table with the flux values for each clump
pah_flux_file = 'dendro_mask_clump_flux_final.txt'
pah_flux = ascii.read(clump_path + pah_flux_file)

# get the table with the filter fluxes
filt_flux_file = 'filter_flux_dendro_clump_mask_corrected.txt'
filt_flux = ascii.read(clump_path + filt_flux_file)

# check to make sure the clumps all have the same columns
clump_order = template['clump_num'] == pah_flux['clump_num']
print(clump_order.all())

if clump_order.all():
    pah_flux.remove_column('clump_num')
else:
    print('Fix the clump order!!!!')

clump_order = template['clump_num'] == filt_flux['clump_num']
print(clump_order.all())

if clump_order.all():
    filt_flux.remove_column('clump_num')
else:
    print('Fix the clump order!!!!')

# combine tables
final_table = hstack([template, pah_flux, filt_flux])

# rename the "flux" keywords to represent the template flux
final_table.rename_column('flux', 'F770W_template_flux')

# add RA and DEC values to the table
# load F770W since it was the template
filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
filename = 'F770W_reproject_to_F1500W_rot.fits'
hdu = fits.open(filepath + filename)
head = hdu[0].header

# create WCS instance
w = WCS(head)

# convert center values to RA and DEC
RA, DEC = w.all_pix2world(final_table['x_cen'], final_table['y_cen'], 0)

# add columns 11, 12
final_table.add_columns([RA, DEC], names = ['ra_cen', 'dec_cen'], indexes = [11,12])

# save the final table
table_name = 'clump_props_pah_filt_flux_corrected.txt'
ascii.write(final_table, clump_path + table_name, overwrite = True)


# # make a quick ds9 region because I love it
# reg_name = 'clump_table_ra_dec_clump_num_ds9'
# reg_file = open(clump_path + reg_name + '.reg', 'w')

# reg_file.write('# Region file format: DS9 version 4.1\n')
# reg_file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
# reg_file.write('fk5\n')

# for i in range(len(final_table)):
#     reg_file.write('point({},{}) # point=circle text={{{:d}}} \n'.format(final_table['ra_cen'][i], final_table['dec_cen'][i], final_table['clump_num'][i]))
    
# reg_file.close()