#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Uses the PDRs4ALL data that Dries sent me to test the PAH contamination in photometric filters.
Will compare to the D21 model results

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import scipy.interpolate as inter
import astropy.stats as stats

os.environ['PYSYN_CDBS']
import pysynphot as S
import stsynphot as stsyn
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import synphot

from pahfit.model import Model
from specutils import Spectrum1D
from astropy.units import Quantity
from astropy.nddata import NDData
from astropy.nddata import StdDevUncertainty

plt.ion()

def split_con(wave, con1, con2, ind):
    ind[np.where((wave < con2) & ((wave > con1)))[0]] = True
    return ind


filt_contam = 'F335M'

if filt_contam == 'F770W':
    wave_low = 4.5    
    wave_high = 10
    
    pah_wave1 = 5
    pah_wave2 = 9.4
    
    filt_list = ['F560W', 'F770W']
    
    p_deg = 1
    
elif filt_contam == 'F335M':
    wave_low = 2
    wave_high = 5
    
    pah_wave1 = 3
    pah_wave2 = 4
    
    filt_list = ['F360M', 'F335M', 'F300M']
    
    p_deg = 1
    
elif filt_contam == 'F1500W':
    wave_low = 10
    wave_high = 20
    
    pah_wave1 = 10.4
    pah_wave2 = 14.8
    
    pah_wave3 = 15.8
    pah_wave4 = 18.2
    
    filt_list = ['F1130W', 'F1500W']
    
    p_deg = 6
    


# filt_list = ['F335M', 'F360M', 'F560W', 'F770W']
miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
filt_dict = {}
filt_wave = np.zeros(len(filt_list))

PDR_spec_arr = ['HII_flux', 'Atomic_flux', 'DF1_flux', 'DF2_flux', 'DF3_flux']

pdrs_filepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/PDRs4ALL_data/'
pahfit_file = 'T1'
orig_spec = 'templates_1147pmap_dries'

fit = ascii.read(pdrs_filepath + pahfit_file + '.ecsv')
data =  ascii.read(pdrs_filepath + orig_spec + '.ecsv')

# loop through filters to setup synthetic photometry later
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

# # an attempt to use the built in PAHFIT tools to plot, then giving up 
# mod = Model.from_saved(pdrs_filepath + pahfit_file + '.ecsv')

# unc = StdDevUncertainty(spec['DF1_unc'].value)
# flux = NDData(data = spec['DF1_flux'].value, uncertainty = unc, unit = spec['DF1_flux'].unit)
# flux_q = Quantity(spec['DF1_flux'])
# wave = Quantity(spec['wavelength'])
# fit_spec = Spectrum1D(flux = flux_q, spectral_axis = wave, redshift = 0, uncertainty = flux.uncertainty)

# fit_spec.meta['instrument'] = 'jwst.nirspec.g140'

# notes from 

# plt.figure(2)
# plt.clf()
# mod.plot(fit_spec)

def remove_PAH(pah1, pah2, ind, wave):
    ind[np.where((wave > pah1) & ((wave < pah2)))[0]] = False
    return ind

Flam_master =  []
Jy_master = []

orig_filt_master = []
orig_filt_master_Jy = []

for PDR_spec in PDR_spec_arr:

    # # plot the raw spectrum just to see it 
    # plt.figure(1)
    # plt.clf()
    # plt.plot(spec['wavelength'], spec[PDR_spec])
    # plt.xlim(3,4)
    # plt.ylim(0, 2100)
    
    # fit and subtract the continuum similar to how was done with the D21 models
    
    data =  ascii.read(pdrs_filepath + orig_spec + '.ecsv')

    # only use the MIR
    cutoff = np.where( (data['wavelength'] > wave_low) &  (data['wavelength'] < wave_high))[0]
    wave = data['wavelength'][cutoff]
    spec = data[PDR_spec][cutoff]
    
    
    plt.figure(1)
    plt.clf()
    plt.plot(wave, spec)
    # plt.ylim(0, 2100)
    
    # convert to Flambda 
    angstrom = wave.to(u.AA)
    spec_val = spec.value * 1e6 * u.Jy
    wave_val = wave.value * u.micron
    Flambda = spec_val.to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(wave_val))
    
    # fit continuum and subtract to isolate PAH features
    ind = np.ones(len(Flambda), dtype = bool)
    ind = remove_PAH(pah_wave1, pah_wave2, ind, wave)
        
    if filt_contam == 'F1500W':
        ind = remove_PAH(pah_wave3, pah_wave4, ind, wave)
        
    clipped = stats.sigma_clip(data = Flambda[ind].value, sigma = 3)
    clipped_mask = clipped.mask
    
    clipped_data = Flambda[ind].value[~clipped_mask]
    clipped_angstrom = angstrom[ind].value[~clipped_mask]
    
    # fit continuum
    pfit = np.polyfit(clipped_angstrom, clipped_data, deg = p_deg)
    p = np.poly1d(pfit)
    xx = np.linspace(min(clipped_angstrom), max(clipped_angstrom))
    
    # subtract continuum
    pah = Flambda.value - p(angstrom.value)
    
    
    ###########################################
    ######## Plotting to check con sub ########
    ###########################################
    
    plt.figure(2)
    plt.clf()
    plt.plot(angstrom, Flambda, c = 'k', lw = 2, label = 'DF1 Spec PDRs4ALL')
    plt.plot(angstrom[ind], Flambda[ind], '.', c = 'orange', label = 'Continuum indices')
    plt.plot(xx, p(xx), ls = '--', c  = 'purple', label = 'Continuum Fit')
    plt.plot(clipped_angstrom, clipped_data, '.', c = 'blue', label = 'Sigma Clipped')
    plt.xlabel('Angstrom')
    plt.ylabel('F_lambda')
    plt.legend(loc = 'best')
    plt.ylim(-1e-5, 5e-5)
    savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/{}_con_fit_tests/'.format(filt_contam)
    savename = '{}_PDRs4ALL_{}_spec_confit.pdf'.format(filt_contam, PDR_spec)
    # plt.savefig(savepath + savename)
    
    plt.figure(3)
    plt.clf()
    plt.plot(angstrom, pah, c = 'b', lw = 2)
    plt.axhline(0, c = 'k', lw = 0.5)
    plt.ylim(-1e-5, 5e-5)
    
    # plt.xlim(0, 15)
    plt.xlabel('Angstrom')
    plt.ylabel('F_lambda')
    plt.title('Continuum subtracted with deg = {} poly'.format(p_deg))
    savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/{}_con_fit_tests/'.format(filt_contam)
    savename = '{}_PDRs4ALL_{}_spec_consub.pdf'.format(filt_contam, PDR_spec)
    # plt.savefig(savepath + savename)
    
    
                
    ###########################################
    ######## Synthetic photometry #############
    ###########################################
       
    # turn the total spectrum into a source spectrum object in synphot 
    sp = synphot.SourceSpectrum(Empirical1D, points = angstrom, lookup_table=pah * units.FLAM, keep_neg = True)
    
    # loop through all filters
    filt_through = np.zeros(len(filt_list))
    filt_through_Jy = np.zeros(len(filt_list))
    

    for f, filt_name in enumerate(filt_list):
        # create an observation and calculate the effective stimulation
        obs = synphot.Observation(sp, filt_dict[filt_name])
        filt_through[f] = obs.effstim(flux_unit='flam').value
        
        # convert back to Jy
        filt_through_Jy[f] = (filt_through[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
        filt_through_Jy[f] = filt_through_Jy[f]/1e6
        
        # save the filter information
    Flam_master.append(filt_through)
    Jy_master.append(filt_through_Jy)
    
    
    filt_Jy = []
    for filt in filt_list:
        filt_Jy.append(filt + '_Jy')
        
        
        
    # apply to the non-continuum subtracted data so we can test S23 method
    sp = synphot.SourceSpectrum(Empirical1D, points = angstrom, lookup_table=Flambda.value * units.FLAM, keep_neg = True)
    
    # loop through all filters
    orig_filt = np.zeros(len(filt_list))
    orig_filt_Jy = np.zeros(len(filt_list))
    

    for f, filt_name in enumerate(filt_list):
        # create an observation and calculate the effective stimulation
        obs = synphot.Observation(sp, filt_dict[filt_name])
        orig_filt[f] = obs.effstim(flux_unit='flam').value
        
        # convert back to Jy
        orig_filt_Jy[f] = (orig_filt[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
        orig_filt_Jy[f] = orig_filt_Jy[f]/1e6
        
        # save the filter information
    orig_filt_master.append(orig_filt)
    orig_filt_master_Jy.append(orig_filt_Jy)
    
    
final_table = np.hstack((Flam_master, Jy_master))
names = filt_list + filt_Jy

orig_filt_master_Jy = np.array(orig_filt_master_Jy)
Jy_master = np.array(Jy_master)


# apply S23 math 
# array is in F360M, F335M, F300M order 
B_pah = 1.6
A_lai = 0.35
B_lai = 0.65

        
xm = orig_filt_master_Jy[:,0]/orig_filt_master_Jy[:,2]
ym = orig_filt_master_Jy[:,1]/orig_filt_master_Jy[:,2]

yc = B_lai * ( (B_pah*xm - ym + A_lai)/(B_pah - B_lai)) + A_lai

F335M_con = yc * orig_filt_master_Jy[:,2]
F335M_pah = orig_filt_master_Jy[:,1] - F335M_con
            
            
plt.figure(4)
plt.clf()
plt.scatter(F335M_pah, Jy_master[:,1])            

# draw y = x line
xvals = np.linspace(0, 800, 100)
plt.plot(xvals, xvals, c = 'k', lw = 0.7)
plt.xlim(0,800)
plt.ylim(0,800)

plt.xlabel('3.3 PAH from S23')
plt.ylabel('3.3 PAH from Continuum Subtraction')

# percent difference
perc = (F335M_pah - Jy_master[:,1])/(Jy_master[:,1])
print(perc)
            
        
