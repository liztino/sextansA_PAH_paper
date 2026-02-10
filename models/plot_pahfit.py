#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plots the PAHFIT results Dries has

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import scipy.interpolate as inter

os.environ['PYSYN_CDBS']
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import synphot

from pahfit.model import Model
import pahfit
from specutils import Spectrum1D
from astropy.units import Quantity
from astropy.nddata import NDData
from astropy.nddata import StdDevUncertainty

plt.ion()

# PDR_spec_arr = ['HII', 'Atomic', 'DF1', 'DF2', 'DF3']
PDR_spec_arr = ['DF2']

pahfit_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/PDRs4ALL_data/orion_template_fits_and_clean_specs/'


for spec in PDR_spec_arr:
    
    spec_name = f'spec_orion_{spec}.ecsv'
    pahfit_name = f'm_orion_{spec}.ecsv'
    
    spec = ascii.read(pahfit_path + spec_name)
    
    mod = Model.from_saved(pahfit_path + pahfit_name)
    
    unc = StdDevUncertainty(spec['col2'].value)
    flux = NDData(data = spec['flux'].value, uncertainty = unc, unit = spec['flux'].unit)
    flux_q = Quantity(spec['flux'])
    wave = Quantity(spec['wavelength'])
    fit_spec = Spectrum1D(flux = flux_q, spectral_axis = wave, redshift = 0, uncertainty = flux.uncertainty)    

    # plt.figure(1, figsize = (10, 5))
    # plt.clf()
    # ax = plt.gca()
    # plt.plot(spec['wavelength'], spec['flux'], drawstyle = 'steps-mid', lw = 2, c = 'k')
    # plt.xlim(2.8, 20)
    
    # pahfit.plot(axs = ax, x = spec['wavelength'], y = spec['flux'], yerr = spec['col2'], model = mod )
    
    fit_spec.meta['instrument'] = 'jwst.*'
    
    mod.plot(fit_spec)

    



