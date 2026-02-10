#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates ds9 region of clump catalog

@author: etarantino
"""

from astropy.wcs import WCS
from astropy.io import ascii, fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u
import glob
import os
import matplotlib as mpl
from astropy.table import QTable, Table, Column

clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
filename = 'clump_catalog.txt'

clumps = ascii.read(clump_path + filename)

# make a ds9 region file

# create region file
reg_name = 'clump_catalog_ds9'
reg_file = open(clump_path + reg_name + '.reg', 'w')

reg_file.write('# Region file format: DS9 version 4.1\n')
reg_file.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
reg_file.write('image\n')


for i in range(len(clumps)):
    reg_file.write('point({},{}) # point=x text={{{}}} color=magenta \n'.format(clumps['x'][i], clumps['y'][i], clumps['Clump'][i]))
    
reg_file.close()