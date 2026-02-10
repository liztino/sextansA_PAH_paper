#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Tests the subtraction from Julia on known models/spectra

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import os

os.environ['PYSYN_CDBS']
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import synphot
from specutils import Spectrum1D
from pahfit.model import Model
from astropy.nddata import StdDevUncertainty
import seaborn as sns
import matplotlib.lines as mlines

# custom functions
from get_pivot_wave import get_pivot_wave
import k_eq
import k_eq_with_err

plt.ion()


# constants
c = 3e14        # in microns/sec
c_A = 3e18      # in Angstrom/sec
h = 6.62e-27    # in erg s 

N_HI = 1e21


# filt_pah = 'F770W'
# filt_list = ['F560W', 'F770W', 'F1000W']
# k_val = 5.84

# filt_pah = 'F335M'
# filt_list = ['F300M', 'F335M', 'F360M']
# k_val = 23

filt_pah = 'F1130W'
filt_list = ['F770W', 'F1000W', 'F1130W']
# k_val = 20
    
# gets the filters 
miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
filt_dict = {}
filt_wave = np.zeros(len(filt_list))

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

# start arrays that we will populate
pah_k_method = []
pah_synphot = []
contam_synphot = []

f1_array = []
f2_array = []
f3_array = []

U_list = np.arange(0, 7, 0.5)
ion_list = ['lo', 'st', 'hi']
size_list = ['sma', 'std', 'lrg']

D21_sp = 'BC03_Z0.0004_10Myr'
time = '1e7'
z = '0.0004'

# U_list = [7,3]
# ion_list = ['hi','lo']
# size_list = ['lrg','sma']

U_master = []
ion_master = []
size_master = []

Flam_master =  []
Jy_master = []

U = 3.0
size  = 'std'
ion = 'st'

# load Draine model
drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
# drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
drainename = 'pahspec.out_bc03_z{:s}_{:s}_{:01.2f}_{:s}_{:s}'.format(z, time, U, ion, size)
# drainename = 'pahspec.out_mmpisrf_{:01.2f}_{:s}_{:s}.gz'.format(U, ion, size)
header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']

# D21 models are given in 4pi nu jnu -> erg s^-1 H^-1
draine_data = ascii.read(drainepath + D21_sp + '/' + drainename, names = header, data_start = 7)

# convert from nu * P_nu to P_nu
Hz = c/(draine_data['wave'])

# extract and convert from erg/s/H to erg/s/cm^2/Hz/sr
# only take the total part, we will use the seperate PAH parts later in the code
total = (draine_data['total'].value/(Hz * 4 * np.pi)) * N_HI
# pah_pos = (draine_data['PAH^+'].value/(Hz * 4 * np.pi)) * N_HI
# pah_neut = (draine_data['PAH^0'].value/(Hz * 4 * np.pi)) * N_HI
# total = pah_pos + pah_neut

# save as an ETC file
draine_mjy = (total * u.erg/u.s/u.cm/u.cm/u.Hz).to(u.Jy)*1e3
ETC = [draine_data['wave'] , draine_mjy]
ETC_name = 'ETC_pahspec.out_bc03_z{:s}_{:s}_{:01.2f}_{:s}_{:s}.txt'.format(z, time, U, ion, size)
names = ['microns', 'mJy']
ascii.write(ETC, drainepath + ETC_name, names = names, overwrite = True)


# only use the MIR
draine_wave = draine_data['wave']            
cutoff = np.where( (draine_wave < 20))[0]
draine_wave = draine_wave[cutoff]
total = total[cutoff]

# convert to Flambda 
angstrom = draine_data['wave'][cutoff].value*1e4
Flambda = total * c_A/(angstrom)**2

# turn the total spectrum into a source spectrum object in synphot 
sp = synphot.SourceSpectrum(Empirical1D, points = angstrom * u.AA, lookup_table=Flambda * units.FLAM, keep_neg = True)

# loop through all filters
filt_through = np.zeros(len(filt_list))
filt_through_Jy = np.zeros(len(filt_list))
count_rate = np.zeros(len(filt_list))

for f, filt_name in enumerate(filt_list):
    # create an observation and calculate the effective stimulation
    obs = synphot.Observation(sp, filt_dict[filt_name])
    filt_through[f] = obs.effstim(flux_unit='flam').value

    
    filt_through_Jy[f] = (filt_through[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
    filt_through_Jy[f] = filt_through_Jy[f]/1e6
    
    # count_rate[f] = obs.effstim('count', area = (25 * u.m * u.m))
    

f1 = filt_through_Jy[0]
f2 = filt_through_Jy[1]
f3 = filt_through_Jy[2]

lam1 = get_pivot_wave(filt_list[0])
lam2 = get_pivot_wave(filt_list[1])
lam3 = get_pivot_wave(filt_list[2])


# # use equation 5 in Gordon+ 2022

# filt_name = 'F770W'
# if filt_name in nircam_filts:
#     filts_dir = '/Users/etarantino/Documents/JWST/filts/nircam/'
# else:
#     filts_dir = '/Users/etarantino/Documents/JWST/filts/miri/'
    
# filt_file = f'{filts_dir}/{filt_name}_mean_system_throughput.txt'
# filt = ascii.read(filt_file)
# filt_wave = (filt['Microns'] * u.micron).to(u.AA).value
# filt_trans = filt['Throughput']

# f_weighted = np.nansum(Flambda*filt_trans*filt_wave)

