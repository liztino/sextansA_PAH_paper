#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Aquires ISO spectra of SMC to see if any regions would be good examples for the 3.3 feature 

@author: etarantino
"""

from astropy.io import fits, ascii
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt

datapath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/ISO_SMC_data/'
dataname = 'NGC104_74803701_sws.tbl'

data = ascii.read(datapath + dataname)

plt.figure(1)
plt.clf()
plt.plot(data['lambda'], data['Flux'])
plt.xlabel('Wavelength ($\mu$m)', size = 'x-large')
plt.ylabel('Flux (Jy)', size = 'x-large')
plt.xlim(0, 20)
plt.ylim(-5, 10)
