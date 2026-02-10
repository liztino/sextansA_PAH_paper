#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

PAH species ratios

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

plt.ion()

clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'

clump_file = 'clump_props_pah_filt_flux_corrected.txt'

clumps = ascii.read(clump_path + clump_file)

sigma_pah = clumps['clump_flux_3_k2'] +clumps['clump_flux_7_k1'] + clumps['clump_flux_11_k1']

plt.figure(1)
plt.clf()
plt.scatter(clumps['clump_num'], clumps['clump_flux_3_k2']/sigma_pah, label = '3.3')
plt.scatter(clumps['clump_num'], clumps['clump_flux_7_k1']/sigma_pah, label = '7.7')
plt.scatter(clumps['clump_num'], clumps['clump_flux_11_k1']/sigma_pah, label = '11.3')
plt.legend(loc='best')