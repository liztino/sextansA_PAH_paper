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
# k_val = 4.33

filt_pah = 'F335M'
filt_list = ['F300M', 'F335M', 'F360M']
k_val = 2.07

# filt_pah = 'F1130W'
# filt_list = ['F1000W', 'F1130W', 'F1500W']
# k_val = 7.21
    
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


PDR_spec_arr = ['HII', 'Atomic', 'DF1', 'DF2', 'DF3']
# PDR_spec_arr = ['DF2']

pahfit_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/PDRs4ALL_data/orion_template_fits_and_clean_specs/'


for spec_val in PDR_spec_arr:
    
    spec_name = f'spec_orion_{spec_val}.ecsv'
    pahfit_name = f'm_orion_{spec_val}.ecsv'
    
    data = ascii.read(pahfit_path + spec_name)
    wave = data['wavelength']
    spec = data['flux']
    
    # remove nans
    nans = np.isnan(spec)
    spec[nans] = 0.0
    
    # convert to Flambda 
    angstrom = wave.to(u.AA)
    spec_val = spec.value * 1e6 * u.Jy
    wave_val = wave.value * u.micron
    Flambda = spec_val.to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(wave_val))
    
    sp = synphot.SourceSpectrum(Empirical1D, points = angstrom, lookup_table=Flambda.value * units.FLAM, keep_neg = True)

    # loop through all filters
    filt_through = np.zeros(len(filt_list))
    filt_through_Jy = np.zeros(len(filt_list))

    for f, filt_name in enumerate(filt_list):
        # create an observation and calculate the effective stimulation
        obs = synphot.Observation(sp, filt_dict[filt_name], force = 'extrap')
        filt_through[f] = obs.effstim(flux_unit='flam').value
        
        filt_through_Jy[f] = (filt_through[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
        filt_through_Jy[f] = filt_through_Jy[f]/1e6
        
    
    f1 = filt_through_Jy[0]
    f2 = filt_through_Jy[1]
    f3 = filt_through_Jy[2]
    
    lam1 = get_pivot_wave(filt_list[0])
    lam2 = get_pivot_wave(filt_list[1])
    lam3 = get_pivot_wave(filt_list[2])
    
    if filt_pah == 'F770W':
        consub = k_eq_with_err.get_pah_low(f1, f2, f3, lam1, lam2, lam3, k_val, 0, 0, 0, 0)
        print(f1, f2, f3)
        print(lam1, lam2, lam3)
        print(consub)
        # print('pah', consub['pah'])
        # print('con', consub['con'])
    else: 
        consub = k_eq_with_err.get_pah_up(f1, f2, f3, lam1, lam2, lam3, k_val, 0, 0, 0, 0)
        print('pah', consub['pah'])
        print('con', consub['con'])
        
    # plt.figure()
    # plt.plot(angstrom, Flambda)
    # plt.title('pah {:5.4f} con {:5.4f}'.format(consub['pah'], consub['con']))


    pah_k_method.append(consub['pah'])
    
    f1_array.append(f1)
    f2_array.append(f2)
    f3_array.append(f3)
            

            
            

# save the files
pah_k_method = np.array(pah_k_method)
pah_synphot = np.array(pah_synphot)
f1_array = np.array(f1_array)
f2_array = np.array(f2_array)
f3_array = np.array(f3_array)


final_table = np.array([PDR_spec_arr, pah_k_method,  f1_array, f2_array, f3_array]).T
# final_table = np.hstack((param_table, pah_k_method, pah_synphot, f1_array, f2_array, f3_array))

name1 = ['reg']
name2 = ['pah_k_method']

# filt_Jy = []
# for filt in filt_list:
#     filt_Jy.append(filt + '_Jy')
    
names = name1 + name2 + filt_list 

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/D21_models/'
final_table_name = 'PDRs4All_' + filt_pah + '_D21_models_with_k_consub_new_pdrs4all_k.txt'
ascii.write(final_table, savepath + final_table_name, names = names, overwrite = True)            
                



# pah_k_method = np.array(pah_k_method)
# pah_synphot = np.array(pah_synphot)
# contam_synphot = np.array(contam_synphot)

# colors = sns.color_palette("hls", 3)

# # translate size text into actual size 
# d = {'sma': 20, 'lrg': 80, 'std': 40}
# size_arr = [d[x] for x in size_master] 

# # translate ionization into colors
# d = {'hi': colors[0], 'lo': colors[2], 'st': colors[1]}
# ion_arr = [d[x] for x in ion_master] 

# plt.figure(1)
# plt.clf()
# plt.scatter(pah_synphot, pah_k_method, s = size_arr, alpha = 0.7, c = ion_arr)
# plt.plot(pah_synphot, pah_synphot)
# plt.xlabel('PAH Flux from Spectra')
# plt.ylabel('PAH Flux from Continuum Sub Math')
# plt.semilogy()
# plt.semilogx()

# hi = mlines.Line2D([], [], color=colors[0], marker='o', ls='', label='Hi', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
# st = mlines.Line2D([], [], color=colors[1], marker='o', ls='', label='St', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
# lo = mlines.Line2D([], [], color=colors[2], marker='o', ls='', label='Lo', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
# # pdrs_label = mlines.Line2D([], [], color='k',  ls='-', label='Empirical: PDRS4ALL', alpha = 1.0)
# # mean_label = mlines.Line2D([], [], color='blue',  ls='-', label='Mean U<4', alpha = 1.0)

# # etc etc
# plt.legend(handles=[ hi, st, lo], loc = 'best')

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/plots/k_method/'
# savename = filt_contam + '_test_k_method_synphot_compare.pdf'
# plt.savefig(savepath + savename)

# perc = (pah_synphot -  pah_k_method)/pah_synphot

# plt.figure(2)
# plt.clf()
# plt.scatter(U_master, perc * 100, s = size_arr, alpha = 0.7, c = ion_arr)
# plt.axhline(0, c = 'k', lw = 0.5)
# plt.xlabel('Ionization Parameter')
# plt.ylabel('Percent Difference')

# plt.legend(handles=[ hi, st, lo], loc = 'best')

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/plots/k_method/'
# savename = filt_contam + '_test_k_method_synphot_percdiff.pdf'
# plt.savefig(savepath + savename)


# k_calc = pah_synphot/contam_synphot

# plt.figure(3)
# plt.clf()
# plt.scatter(U_master, k_calc, s = size_arr, alpha = 0.7, c = ion_arr)
# plt.axhline(k_val, c = 'b', lw = 2)
# plt.xlabel('Ionization Parameter')
# plt.ylabel('k from synphot models')

# plt.legend(handles=[ hi, st, lo], loc = 'best')

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/plots/k_method/'
# savename = filt_contam + '_test_k_method_synphot_k_vals.pdf'
# plt.savefig(savepath + savename)

