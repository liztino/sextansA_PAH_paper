#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script to make figure for paper that shows the contamination on the filters with D21 models and PDRs4All

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter
from astropy.table import Table
import os, glob
import seaborn as sns
import astropy.units as u

os.environ['PYSYN_CDBS']
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import synphot

# constants
c = 3e14        # in microns/sec
c_A = 3e18      # in Angstrom/sec
h = 6.62e-27    # in erg s 

# column density in Sextans A to normalize the Draine models
# Max integrated HI flux near reg in Jy/beam * m/s 
HI_flux = 250

# beam in arcseconds
bmaj = 7.6
bmin = 6.5

# calculate the column density 
N_HI = (1.105e21 * HI_flux)/(bmaj * bmin)


filt_list = ['F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W']

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W', 'F2100W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
filt_dict = {}
filt_wave = np.zeros(len(filt_list))

# load the filter transmission curves and info 
for i, filt_name in enumerate(filt_list):
    # load filter info
    if filt_name in nircam_filts:
        filts_dir = '/Users/etarantino/Documents/JWST/filts/nircam/'
    else:
        filts_dir = '/Users/etarantino/Documents/JWST/filts/miri/'
        
    filt_file = f'{filts_dir}/{filt_name}_mean_system_throughput.txt'
    filt = ascii.read(filt_file)

    # create filter object in synphot
    filt_feature = synphot.SpectralElement(Empirical1D, points = filt['Microns']*1e4, lookup_table = filt['Throughput'], keep_neg = True)

    # save to the filter dictionary
    filt_dict[filt_name] = filt_feature

    filt_wave[i] = filt_feature.pivot().value
    
# load in a D21 model
drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
model = 'BC03_Z0.0004_10Myr'

U = 1
ion = 'st'
size = 'std'
drainename = 'pahspec.out_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}'.format(U, ion, size)
header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']

# D21 models are given in 4pi nu jnu -> erg s^-1 H^-1
draine_data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)

# convert from nu * P_nu to P_nu
Hz = c/(draine_data['wave'])

# extract and convert from erg/s/H to erg/s/cm^2/Hz/sr
total = (draine_data['total'].value/(Hz * 4 * np.pi)) * N_HI
(total * u.erg / u.s / u.cm**2 / u.Hz).to(u.Jy)

# load a PDRs4All spectrum 
pdrs_filepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/PDRs4ALL_data/'
orig_spec = 'templates_1147pmap_dries'
pdrs_data = ascii.read(pdrs_filepath + orig_spec + '.ecsv')

PDR_spec_arr = ['HII_flux', 'Atomic_flux', 'DF1_flux', 'DF2_flux', 'DF3_flux']
wave = pdrs_data['wavelength']
spec = pdrs_data['DF1_flux']

###########################
####### Now plot ##########
###########################

plt.figure(1, figsize = (13, 20))
plt.clf()
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True, num = 1)

plt.subplots_adjust(hspace = 0)
xlim = (2, 20)

ax1.plot(draine_data['wave'], total)
ax1.set_xlim(xlim)
ax1.set_ylim()



