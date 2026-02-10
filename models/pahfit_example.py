#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 12:28:56 2024

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u


from pahfit.model import Model
from specutils import Spectrum1D
from astropy.units import Quantity
from astropy.nddata import NDData
from astropy.nddata import StdDevUncertainty

plt.ion()

pdrs_filepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/PDRs4ALL_data/'
pahfit_file = 'T1'
orig_spec = 'templates_1147pmap_dries'

fit = ascii.read(pdrs_filepath + pahfit_file + '.ecsv')
spec =  ascii.read(pdrs_filepath + orig_spec + '.ecsv')

# an attempt to use the built in PAHFIT tools to plot, then giving up 
mod = Model.from_saved(pdrs_filepath + pahfit_file + '.ecsv')

unc = StdDevUncertainty(spec['DF1_unc'].value)
flux = NDData(data = spec['DF1_flux'].value, uncertainty = unc, unit = spec['DF1_flux'].unit)
flux_q = Quantity(spec['DF1_flux'])
wave = Quantity(spec['wavelength'])
fit_spec = Spectrum1D(flux = flux_q, spectral_axis = wave, redshift = 0, uncertainty = flux.uncertainty)

fit_spec.meta['instrument'] = 'jwst.nirspec.g140.high'

mod = Model.from_saved(pdrs_filepath + pahfit_file + '.ecsv')

# notes from Dries
tot = mod.tabulate(wavelengths=wave, instrumentname='jwst.nirspec.g140.high', redshift=0)
pah = mod.tabulate(wavelengths=wave, instrumentname='jwst.nirspec.g140.high', redshift=0, feature_mask = mod.features['kind']=='dust_feature') 
con = mod.tabulate(wavelengths=wave, instrumentname='jwst.nirspec.g140.high', redshift=0, feature_mask = mod.features['kind']=='dust_continuum')
star = mod.tabulate(wavelengths=wave, instrumentname='jwst.nirspec.g140.high', redshift=0, feature_mask = mod.features['kind']=='starlight')


plt.figure(1)
plt.clf()
# plt.plot(tot.wavelength, tot.flux)
# plt.plot(pah.wavelength, pah.flux)
# plt.plot(con.wavelength, con.flux)
plt.plot(star.wavelength, star.flux)

# create an instrument pack for the Draine 2021 models

U = 1.0
ion = 'hi'
size = 'lrg'

# Max integrated HI flux near reg in Jy/beam * m/s 
HI_flux = 250

# beam in arcseconds
bmaj = 7.6
bmin = 6.5

# calculate the column density 
N_HI = (1.105e21 * HI_flux)/(bmaj * bmin)

# convert to MJy
erg_MJy = 1e-23

c = 3e14        # in microns/sec


# load Draine model
drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
model = 'BC03_Z0.0004_10Myr'
# drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
drainename = 'pahspec.out_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}'.format(U, ion, size)
header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']

# D21 models are given in 4pi nu jnu -> erg s^-1 H^-1
draine_data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)
print('wave start', draine_data['wave'][0])

ind = np.where(draine_data['wave'] > 40)[0][0]

print(ind)
print('wave end', draine_data['wave'][ind])

# convert to useable units 
micron = draine_data['wave'][:ind]
Hz = c/(micron)
total = ((draine_data['total'][:ind]/(Hz * 4 * np.pi)) * N_HI)/erg_MJy

plt.figure(2)
plt.clf()
plt.plot(micron, total )

# plot the spectral resolution 
R = (micron[:-1])/np.diff(micron)

print('R', np.mean(R))

plt.figure(3)
plt.clf()
plt.plot(micron[:-1], R)
plt.xlabel('Wavelength')
plt.ylabel('R')

# # attempt to fit with new instrument pack 
unc = StdDevUncertainty(0.01 * total)
# flux = NDData(data = total, uncertainty = unc, unit = u.Jy / u.sr )
flux_q = Quantity(total * u.Jy / u.sr )
wave = Quantity(micron * u.micron)
fit_spec = Spectrum1D(flux = flux_q, spectral_axis = wave, redshift = 0, uncertainty = unc)

fit_spec.meta['instrument'] = 'draine21*'
model = Model.from_yaml('classic_nir.yaml')
model.guess(fit_spec)
model.fit(fit_spec, redshift = 0, verbose = True)
model.plot(fit_spec)

tot = model.tabulate(wavelengths=wave, instrumentname='draine21*', redshift=0)
pah = model.tabulate(wavelengths=wave, instrumentname='draine21*', redshift=0, feature_mask = model.features['kind']=='dust_feature') 
con = model.tabulate(wavelengths=wave, instrumentname='draine21*', redshift=0, feature_mask = model.features['kind']=='dust_continuum')
star = model.tabulate(wavelengths=wave, instrumentname='draine21*', redshift=0, feature_mask = model.features['kind']=='starlight')

plt.figure(5)
plt.clf()
plt.plot(tot.wavelength, tot.flux, c = 'k', label = 'Total')
plt.plot(pah.wavelength, pah.flux, c = 'm', label = 'PAH')
plt.plot(con.wavelength, con.flux, c = 'g', label = 'Continuum')
plt.plot(star.wavelength, star.flux, c = 'b', label = 'Stellar')

plt.legend(loc='best')
