#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 18:07:05 2025

@author: etarantino
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.io import ascii
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord

# load in main clump table
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
table_name = 'clump_props_and_flux.txt'
clump_table = ascii.read(clump_path + table_name)

# flux units are in Jy, should convert to micro Janskys!

# # open the textfile for the new table
# tabpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/manuscript/tables/'
# tabname = 'clump_prop_table_main.tex'

# f = open(tabpath + tabname, 'w+')

# # preamble set up
# f.write('\\begin{table}[h] \n')
# f.write('\\caption{Clump Properties}\\label{tab:clump_prop} \n')
# f.write('\\begin{tabular}{ccccccccc} \n')
# f.write('\\toprule \n')
# f.write('Clump & RA & Dec. & Radius & Major Axis & Minor Axis & 3.3\,$\mu$m Flux & 7.7\,$\mu$m Flux & 11.3\,$\mu$m Flux \\\\ \n' )
# f.write('& (10h 11m s) & (-04$^{\circ} \ 42\' \ \'\'$) & (pc) & (pc) & (pc) & ($\mu$Jy) & ($\mu$Jy) & ($\mu$Jy) \\\\ \n' )
# f.write('\\midrule \n')

# - d_sexa = 1.4 Mpc
# - 1 arcsec = 6.787 pc = 6.8 pc

# conversion from degrees to pc
to_pc = 6.787 * 3600

# now loop through and add a row
for i in range(len(clump_table)):
    
    # convert to micro janskys
    flux3 = clump_table['clump_flux_3_k1'][i] * 1e6
    flux3_err = clump_table['clump_flux_3_k1_err_k'][i] * 1e6
    
    # convert to micro janskys
    flux7 = clump_table['clump_flux_7_k1'][i] * 1e6
    flux7_err = clump_table['clump_flux_7_k1_err_k'][i] * 1e6
    
    # convert to micro janskys
    flux11 = clump_table['clump_flux_11_k1'][i] * 1e6
    flux11_err = clump_table['clump_flux_11_k1_err_k'][i] * 1e6
    
    # lenght measures are in degrees, need to convert to pc
    r = clump_table['radius'][i] * to_pc
    major = clump_table['major_sigma'][i] * to_pc
    minor = clump_table['minor_sigma'][i] * to_pc
    
    # convert to hms and dms
    c = SkyCoord(clump_table['ra_cen'][i], clump_table['dec_cen'][i], unit = 'deg')
    ra_str, dec_str = c.to_string('hmsdms', precision = 3).split(' ')

    ra_str = ra_str.split('m')[-1][0:-1]
    dec_str = dec_str.split('m')[-1][0:-1]


#     f.write('{:d} & {:s} & {:s} & {:1.2f} & {:1.2f} & {:1.2f} & {:3.2f} $\pm$ {:1.2f} & {:3.2f} $\pm$ {:1.2f} & {:3.2f} $\pm$ {:1.2f} \\\\\n'.format(clump_table['clump_num'][i], ra_str, dec_str, r, major, minor, flux3, flux3_err, flux7, flux7_err, flux11, flux11_err))


# f.write('\\botrule \n')
# f.write('\\end{tabular} \n')
# f.write('\\end{table} \n')

# f.close()