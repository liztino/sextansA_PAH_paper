#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Attempts to make a corner plot with results from models

@author: etarantino
"""

import corner
import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
import os
import astropy.units as u

# ndim, nsamples = 4, 50000
# np.random.seed(1234)
# data1 = np.random.randn(ndim * 4 * nsamples // 5).reshape(
#     [4 * nsamples // 5, ndim]
# )
# mean = 4 * np.random.rand(ndim)
# data2 = mean[None, :] + np.random.randn(ndim * nsamples // 5).reshape(
#     [nsamples // 5, ndim]
# )
# samples = np.vstack([data1, data2])

# figure = corner.corner(samples)

plt.ion()
savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/D21_new_phoenix/'
filename = 'D21_combo_SexA_pah_small_box_models_master_table_v2.txt'
plotdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/D21_new_phoenix/'

modeldir = savedir + 'spec/'
throughdir = savedir + 'throughput/'

if not os.path.exists(throughdir):
    os.makedirs(throughdir)    

if not os.path.exists(modeldir):
    os.makedirs(modeldir)  

data = ascii.read(savedir + filename)

ind = np.where(data['chisq_Jy'] < 1000)
data = data[ind]

# samples = np.array([data['U'], data['qpah'], data['const'], data['chisq_Jy']]).T

samples = np.array([data['U'], data['qpah'], data['const'], data['chisq_Jy']]).T

figure = corner.corner(samples, bins = 200)

