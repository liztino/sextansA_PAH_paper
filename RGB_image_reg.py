#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:14:29 2023

@author: etarantino
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from astropy.visualization import make_lupton_rgb

nircam_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/nircam_mar/mosaic/'
miri_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/nircam_mar/mosaic/'

filts = 

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W'}
nircam_filts = {'f115w', 'f150w', 'f200w', 'f300m', 'f335m', 'f360m'}




fig = plt.figure(1)
plt.clf()
ax = plt.subplot(projection=wcs, figure = fig)
