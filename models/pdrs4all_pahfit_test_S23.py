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
    wave_low = 4    
    wave_high = 10.5
    
    pah_wave1 = 5
    pah_wave2 = 9.4
    
    filt_list = ['F560W', 'F770W']
    
elif filt_contam == 'F335M':
    wave_low = 2
    wave_high = 5
    
    pah_wave1 = 3
    pah_wave2 = 4
    
    filt_list = ['F360M', 'F335M', 'F300M']
    
elif filt_contam == 'F1500W':
    wave_low = 9.4
    wave_high = 20
    
    pah_wave1 = 10.5
    pah_wave2 = 14.5
    
    pah_wave3 = 16
    pah_wave4 = 18
    
    filt_list = ['F1130W', 'F1500W']

# filt_list = ['F335M', 'F360M', 'F560W', 'F770W']
miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
filt_dict = {}
filt_wave = np.zeros(len(filt_list))

PDR_spec_arr = ['HII_flux', 'Atomic_Flux', 'DF1_Flux', 'DF2_flux', 'DF3_flux']
pahfit_file_list = ['T1', 'T2', 'T3', 'T4', 'T5']

pdrs_filepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/PDRs4ALL_data/'
pahfit_file = 'T1'
orig_spec = 'templates_1147pmap_dries'

fit = ascii.read(pdrs_filepath + pahfit_file + '.ecsv')
spec =  ascii.read(pdrs_filepath + orig_spec + '.ecsv')

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


for pahfit_file in pahfit_file_list:
    
    mod = Model.from_saved(pdrs_filepath + pahfit_file + '.ecsv')
    
    print(pahfit_file, mod.features.loc['PAH_3.3']['power'][0])

    tot = mod.tabulate(wavelengths=spec['wavelength'], instrumentname='jwst.nirspec.*', redshift=0)
    pah = mod.tabulate(wavelengths=spec['wavelength'], instrumentname='jwst.nirspec.*', redshift=0, feature_mask = mod.features['kind']=='dust_feature') 
    con = mod.tabulate(wavelengths=spec['wavelength'], instrumentname='jwst.nirspec.*', redshift=0, feature_mask = mod.features['kind']=='dust_continuum')
    star = mod.tabulate(wavelengths=spec['wavelength'], instrumentname='jwst.nirspec.*', redshift=0, feature_mask = mod.features['kind']=='starlight')


    plt.figure(1)
    plt.clf()
    plt.plot(tot.wavelength.to(u.micron), tot.flux, c = 'k', label = 'Total')
    plt.plot(pah.wavelength.to(u.micron), pah.flux, c = 'm', label = 'PAH')
    plt.plot(con.wavelength.to(u.micron), con.flux, c = 'g', label = 'Continuum')
    plt.plot(star.wavelength.to(u.micron), star.flux, c = 'b', label = 'Stellar')
    
    plt.xlim(1.5, 6)
    plt.ylim(0, max(pah.flux.value))

    plt.legend(loc='best')
    
    savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/{}_con_fit_tests/'.format(filt_contam)
    savename = 'PAHFIT_{}_PDRs4ALL_{}_spec_consub.pdf'.format(filt_contam, pahfit_file)
    plt.savefig(savepath + savename)

                
    # just use the fits that Dries made to find the PAH part
    # no need to do continuum subtraction for the F335M filter
    pah_flux_lam = (pah.flux).to(u.erg / u.s / u.AA / u.cm**2 / u.sr, equivalencies=u.spectral_density(pah.wavelength))

    pah_lam = (pah.wavelength).to(u.AA, equivalencies=u.spectral_density(pah.wavelength))
    
    tot_flam = (tot.flux).to(u.erg / u.s / u.AA / u.cm**2 / u.sr, equivalencies=u.spectral_density(pah.wavelength))
    
    
    ###########################################
    ######## Synthetic photometry #############
    ###########################################
       
    # turn the total spectrum into a source spectrum object in synphot 
    sp = synphot.SourceSpectrum(Empirical1D, points = pah_lam, lookup_table=pah_flux_lam.value * units.FLAM, keep_neg = True)
    
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
    sp = synphot.SourceSpectrum(Empirical1D, points = pah_lam, lookup_table=tot_flam.value * units.FLAM, keep_neg = True)
    
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
            
            
plt.figure(5)
plt.clf()
plt.scatter(F335M_pah, Jy_master[:,1])            

# draw y = x line
xvals = np.linspace(0, 1000, 100)
plt.plot(xvals, xvals, c = 'k', lw = 0.7)
plt.xlim(0,1000)
plt.ylim(0,1000)

plt.xlabel('3.3 PAH from S23 [MJy/sr]')
plt.ylabel('3.3 PAH from PAHFIT [MJy/sr]')

# percent difference
perc = (F335M_pah - Jy_master[:,1])/(Jy_master[:,1])
print(perc)
            
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/plots/'   
savename = 'S23_method_on_PDRs4all_pahfit.pdf'
plt.savefig(savepath + savename)            




            
        
