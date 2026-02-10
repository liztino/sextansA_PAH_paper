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
from scipy.interpolate import interp1d, PchipInterpolator

# os.environ['PYSYN_CDBS']
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import synphot

from pahfit.model import Model
from specutils import Spectrum1D
from astropy.units import Quantity
from astropy.nddata import NDData
from astropy.nddata import StdDevUncertainty

import k_eq



plt.ion()

def split_con(wave, con1, con2, ind):
    ind[np.where((wave < con2) & ((wave > con1)))[0]] = True
    return ind


filt_contam = 'F770W'

if filt_contam == 'F770W':
    wave_low = 2.8
    wave_high = 9
    
    pah_wave1 = 5
    pah_wave2 = 9.4
    
filt_list = ['F444W', 'F560W', 'F770W', 'F1000W']

# filt_list = ['F335M', 'F360M', 'F560W', 'F770W']
miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M', 'F444W'}
filt_dict = {}
filt_wave = np.zeros(len(filt_list))

PDR_spec_arr = ['HII', 'Atomic', 'DF1', 'DF2']
# PDR_spec_arr = ['DF2']

pahfit_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/PDRs4ALL_data/orion_template_fits_and_clean_specs/'


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

for spec_val in PDR_spec_arr:
    
    spec_name = f'spec_orion_{spec_val}.ecsv'
    pahfit_name = f'm_orion_{spec_val}.ecsv'
    
    spec = ascii.read(pahfit_path + spec_name)
    
    mod = Model.from_saved(pahfit_path + pahfit_name)
    
    unc = StdDevUncertainty(spec['col2'].value)
    flux = NDData(data = spec['flux'].value, uncertainty = unc, unit = spec['flux'].unit)
    flux_q = Quantity(spec['flux'])
    wave = Quantity(spec['wavelength'])
    fit_spec = Spectrum1D(flux = flux_q, spectral_axis = wave, redshift = 0, uncertainty = flux.uncertainty)    
    
    fit_spec.meta['instrument'] = 'jwst.*'
    
    mod.plot(fit_spec)
    
    fig = plt.gcf()
    ax_list = fig.axes
    
    ylim_xval = np.where((spec['wavelength'] > 7) & (spec['wavelength'] < 8))[0]
    norm_val = np.nanmax(spec['flux'][ylim_xval])
    ax_list[0].set_ylim(0, norm_val + 0.2* norm_val)
    
    # ax_list[0].set_xlim(wave_low, wave_high)
    
    # plot transmission curves
    # ax_list[0].plot(spec['wavelength'], norm_val*filt_dict['F360M'](spec['wavelength'].value * u.micron), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    # ax_list[0].plot(spec['wavelength'], norm_val*filt_dict['F300M'](spec['wavelength'].value * u.micron), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    # ax_list[0].plot(spec['wavelength'], norm_val*filt_dict['F335M'](spec['wavelength'].value * u.micron), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    ax_list[0].plot(spec['wavelength'], 2*norm_val*filt_dict['F560W'](spec['wavelength'].value * u.micron), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    ax_list[0].plot(spec['wavelength'], norm_val*filt_dict['F444W'](spec['wavelength'].value * u.micron), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    ax_list[0].plot(spec['wavelength'], 2*norm_val*filt_dict['F1000W'](spec['wavelength'].value * u.micron), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    ax_list[0].plot(spec['wavelength'], 2*norm_val*filt_dict['F770W'](spec['wavelength'].value * u.micron), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1) 
    # ax_list[0].plot(spec['wavelength'], 2*norm_val*filt_dict['F1500W'](spec['wavelength'].value * u.micron), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1) 
    # ax_list[0].plot(spec['wavelength'], 2*norm_val*filt_dict['F1130W'](spec['wavelength'].value * u.micron), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1) 

    tot = mod.tabulate(wavelengths=spec['wavelength'], instrumentname='jwst*', redshift=0)
    pah = mod.tabulate(wavelengths=spec['wavelength'], instrumentname='jwst*', redshift=0, feature_mask = mod.features['kind']=='dust_feature') 
    con = mod.tabulate(wavelengths=spec['wavelength'], instrumentname='jwst*', redshift=0, feature_mask = mod.features['kind']=='dust_continuum')
    star = mod.tabulate(wavelengths=spec['wavelength'], instrumentname='jwst*', redshift=0, feature_mask = mod.features['kind']=='starlight')
                    
    # just use the fits that Dries made to find the PAH part
    # no need to do continuum subtraction for the F335M filter
    pah_flux_lam = (pah.flux).to(u.erg / u.s / u.AA / u.cm**2 / u.sr, equivalencies=u.spectral_density(pah.wavelength))

    pah_lam = (pah.wavelength).to(u.AA, equivalencies=u.spectral_density(pah.wavelength))
    
    tot_flux_lam = (flux_q).to(u.erg / u.s / u.AA / u.cm**2 / u.sr, equivalencies=u.spectral_density(pah.wavelength))

    tot_lam = (wave).to(u.AA, equivalencies=u.spectral_density(pah.wavelength))
    
    
    ###########################################
    ######## Synthetic photometry #############
    ###########################################
       
    # turn the total spectrum into a source spectrum object in synphot 
    # pah_sp = synphot.SourceSpectrum(Empirical1D, points = pah_lam, lookup_table=pah_flux_lam.value * units.FLAM, keep_neg = True)
    
    # clean out the nans before sending over the spectrum
    nan_ind = np.where(np.isnan(tot_flux_lam))[0]    
    if nan_ind.any():
    
        ind1 = nan_ind[0] - int(len(nan_ind))*2
        ind2 = nan_ind[-1] + int(len(nan_ind))*2
        
        mask = np.zeros(len(wave), dtype = bool)
        mask[ind1:(nan_ind[0]-5)] = True
        mask[(nan_ind[-1]+5):ind2] = True
        
        # spline_inter = PchipInterpolator(wave[mask], tot_flux_lam[mask])
        spline_inter = interp1d(wave[mask], tot_flux_lam[mask], kind = 'nearest')
        
        tot_flux_lam[nan_ind] = spline_inter(wave[nan_ind]) * tot_flux_lam.unit
    
        # vals = np.linspace(wave[ind1], wave[ind2-1], 1000)
        
        # plt.figure()
        # plt.plot(vals, spline_inter(vals), label = 'Interp', alpha = 0.5, ls = '--')
        # plt.plot(wave, tot_flux_lam, label = 'data', alpha = 0.5)
        # plt.legend(loc = 'best')

    
    tot_sp = synphot.SourceSpectrum(Empirical1D, points = tot_lam, lookup_table=tot_flux_lam.value * units.FLAM, keep_neg = True)
    # loop through all filters
    filt_through = np.zeros(len(filt_list))
    filt_through_Jy = np.zeros(len(filt_list))
    syn_flux = {}    

    for f, filt_name in enumerate(filt_list):
        # create an observation and calculate the effective stimulation
        obs = synphot.Observation(tot_sp, filt_dict[filt_name])
        filt_through[f] = obs.effstim(flux_unit='flam').value
        
        # convert back to Jy
        filt_through_Jy[f] = (filt_through[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
        filt_through_Jy[f] = filt_through_Jy[f]/1e6
        
        syn_flux[filt_list[f]] = filt_through_Jy[f]
        
    # # save the filter information
    # Flam_master.append(filt_through)
    # Jy_master.append(filt_through_Jy)
    
    # find the slope just using the pivot wavelengths
    # wavelength range
    
    delta_wave = (filt_dict['F1000W'].pivot().to(u.micron).value  - filt_dict['F444W'].pivot().to(u.micron).value)
    slope = (filt_through_Jy[3] -  filt_through_Jy[0])/delta_wave
    
    con_wave = np.linspace(3.5, 12, 1000)
    con = slope * (con_wave - filt_dict['F444W'].pivot().to(u.micron).value) + filt_through_Jy[0]
    
    ax_list[0].plot(con_wave, con, c = 'm', lw = 3, zorder = 10000)
    
    # now plot what we do for F770W and the k method
    con_wave = np.linspace(4.5, 12, 1000)
    # beta = (con_wave  - filt_dict['F560W'].pivot().to(u.micron).value)/(filt_dict['F1000W'].pivot().to(u.micron).value  - filt_dict['F560W'].pivot().to(u.micron).value)
    beta = (filt_dict['F770W'].pivot().to(u.micron).value  - filt_dict['F560W'].pivot().to(u.micron).value)/(filt_dict['F1000W'].pivot().to(u.micron).value  - filt_dict['F560W'].pivot().to(u.micron).value)
    
    k = 4.33
    con = ((k / (k+beta - 1)) * (syn_flux['F770W'] - (1 - beta) * syn_flux['F560W'] - beta * syn_flux['F1000W'])) - syn_flux['F770W']
    
    # get result from k method
    f1 = syn_flux['F560W']
    f2 = syn_flux['F770W']
    f3 = syn_flux['F1000W']
    lam1 = filt_dict['F560W'].pivot().to(u.micron).value
    lam2 = filt_dict['F770W'].pivot().to(u.micron).value
    lam3 = filt_dict['F1000W'].pivot().to(u.micron).value
    
    fp1 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(1 - k - ((lam2 - lam1)/(lam3 - lam1)))
    fp2 = fp1*k
    fc1 = f1 - fp1
    fc2 = (((f3 - fc1)/(lam3 - lam1)) * (lam2 - lam1)) + fc1
    slope = (((f3 - fc1)/(lam3 - lam1)))
    
    con = (((f3 - fc1)/(lam3 - lam1)) * (con_wave - lam1)) + fc1
    
    ax_list[0].plot(con_wave, con, c = 'lime', lw = 3, zorder = 10000)
    
    # now try powerlaw from Donneley+ 2025
    con_wave = np.linspace(4.5, 12, 1000)
    alpha = (np.log10(con_wave) - np.log10(filt_dict['F444W'].pivot().to(u.micron).value))/(np.log10(filt_dict['F1000W'].pivot().to(u.micron).value)  - np.log10(filt_dict['F444W'].pivot().to(u.micron).value))
    
    con = syn_flux['F444W']**(1 - alpha) * syn_flux['F1000W']**alpha
    
    ax_list[0].plot(con_wave, con, c = 'orange', lw = 3, zorder = 10000)
    
    # # add the F560W powerlaw
    # con_wave = np.linspace(4.5, 12, 1000)
    # alpha = (np.log10(con_wave) - np.log10(filt_dict['F560W'].pivot().to(u.micron).value))/(np.log10(filt_dict['F1000W'].pivot().to(u.micron).value)  - np.log10(filt_dict['F560W'].pivot().to(u.micron).value))
    
    # con = 0.91* syn_flux['F560W']**(1 - alpha) * syn_flux['F1000W']**alpha
    
    # ax_list[0].plot(con_wave, con, c = 'red', lw = 3, zorder = 10000)

    
    # title and labels
    ax_list[0].set_title('{:s}'.format(spec_val))
    ax_list[0].set_xlim(3, 15)
    
    # # also plot it on a normal spectrum
    # plt.figure()
    # plt.scatter(spec['wavelength'], spec['flux'])
    # plt.plot(con_wave, con, c = 'm', lw = 3)
    

    savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/other_consub/'
    savename = 'PAHFIT_PDRs4ALL_{}_spec_consub_with_F444W.pdf'.format(spec_val)
    plt.savefig(savepath + savename)

    
    
    
    
    
#     filt_Jy = []
#     for filt in filt_list:
#         filt_Jy.append(filt + '_Jy')
    
# final_table = np.hstack((Flam_master, Jy_master))
# names = filt_list + filt_Jy

# print(filt_Jy)
# print(filt_through_Jy)

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/calc_PAH_contamination/'
# final_table_name = f'{filt_contam}_PAHFIT_PDRs4ALL_spec_PAH_contamination_v4.txt'
# # ascii.write(final_table, savepath + final_table_name, names = names, overwrite = True)            
                

            
            
            




            
        
