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
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import synphot

# from pahfit.model import Model
# from specutils import Spectrum1D
from astropy.units import Quantity
from astropy.nddata import NDData
from astropy.nddata import StdDevUncertainty

plt.ion()

def split_con(wave, con1, con2, ind):
    ind[np.where((wave < con2) & ((wave > con1)))[0]] = True
    return ind

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

filt_contam = 'F335M'

if filt_contam == 'F770W':
    wave_low = 3.7    
    wave_high = 10
    
    pah_wave1 = 5
    pah_wave2 = 9.4
    
    filt_list = ['F770W', 'F560W']
    
    p_deg = 1
    
elif filt_contam == 'F335M':
    wave_low = 2
    wave_high = 5
    
    pah_wave1 = 3
    pah_wave2 = 4
    
    filt_list = ['F360M', 'F335M']
    
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

U_list = [0, 1, 2, 3, 4]
size_list = ['sma', 'std', 'lrg', 'sma', 'std']
ion_list = ['lo', 'st', 'hi', 'hi', 'lo']

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

###########################################
######## Synthetic photometry #############
###########################################

def get_synphot(angstrom, pah, filt_list):
   
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
    
    return filt_through_Jy


def remove_PAH(pah1, pah2, ind, wave):
    ind[np.where((wave > pah1) & ((wave < pah2)))[0]] = False
    return ind

Flam_master =  []
Jy_master = []

method = 'new'

for i, PDR_spec in enumerate(PDR_spec_arr):

    # # plot the raw spectrum just to see it 
    # plt.figure(1)
    # plt.clf()
    # plt.plot(spec['wavelength'], spec[PDR_spec])
    # plt.xlim(3,4)
    # plt.ylim(0, 2100)
    
    # fit and subtract the continuum similar to how was done with the D21 models
    if method == 'new':
        data =  ascii.read(pdrs_filepath + orig_spec + '.ecsv')
    
        # only use the MIR
        cutoff = np.where( (data['wavelength'] > wave_low) &  (data['wavelength'] < wave_high))[0]
        wave = data['wavelength'][cutoff]
        spec = data[PDR_spec][cutoff]
        
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
        
        # data =  ascii.read(pdrs_filepath + orig_spec + '.ecsv')

        # # only use the MIR
        # cutoff = np.where( (data['wavelength'] > wave_low) &  (data['wavelength'] < wave_high))[0]
        # wave = data['wavelength'][cutoff]
        # spec = data[PDR_spec][cutoff]
    
        # plt.figure(1)
        # plt.clf()
        # plt.plot(wave, spec)
        # # plt.ylim(0, 2100)
        
        # # convert to Flambda 
        # angstrom = wave.to(u.AA)
        # spec_val = spec.value * 1e6 * u.Jy
        # wave_val = wave.value * u.micron
        # Flambda = spec_val.to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(wave_val))
        
        # # fit continuum and subtract to isolate PAH features
        # ind = np.ones(len(Flambda), dtype = bool)
        # ind = remove_PAH(pah_wave1, pah_wave2, ind, wave)
            
        # if filt_contam == 'F1500W':
        #     ind = remove_PAH(pah_wave3, pah_wave4, ind, wave)
            
        # clipped = stats.sigma_clip(data = Flambda[ind].value, sigma = 3)
        # clipped_mask = clipped.mask
        
        # clipped_data = Flambda[ind].value[~clipped_mask]
        # clipped_angstrom = angstrom[ind].value[~clipped_mask]
        
        # # fit continuum
        # pfit = np.polyfit(clipped_angstrom, clipped_data, deg = p_deg)
        # p = np.poly1d(pfit)
        # xx = np.linspace(min(clipped_angstrom), max(clipped_angstrom))
        
        # # subtract continuum
        # pah = Flambda.value - p(angstrom.value)
        
        
    else:
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
        
        
        def split_con(con1, con2, ind):
            ind[np.where((draine_wave < con2) & ((draine_wave > con1)))[0]] = True
            return ind
    
        # only use the MIR
        draine_wave = wave_val.value           
        cutoff = np.where(draine_wave < 20)[0]
        draine_wave = draine_wave[cutoff]
        spec_val = spec_val[cutoff]
        
        # convert to Flambda 
        angstrom = wave_val[cutoff].value*1e4
        Flambda = Flambda[cutoff]
        
        # fit continuum and subtract to isolate PAH features
        ind = np.zeros(len(spec_val), dtype = bool)
        
        # isolate continuum indices
        ind = split_con(9.4, 10.2, ind)
        ind = split_con(14.5, 16, ind)
        ind = split_con(18, 20, ind)
        
        # fit continuum
        pfit = np.polyfit(angstrom[ind], Flambda[ind], deg = 4)
        p = np.poly1d(pfit)
        xx = np.linspace(min(angstrom), max(angstrom))
        
        # subtract continuum
        pah = Flambda.value - p(angstrom)
        angstrom = angstrom * u.AA
    
    # load a D21 model and do the continuum subtraction 
    # load Draine model
    U = U_list[i]
    ion = ion_list[i]
    size = size_list[i]
    drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
    model = 'BC03_Z0.0004_10Myr'
    # drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
    drainename = 'pahspec.out_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}'.format(U, ion, size)
    header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']
    
    # D21 models are given in 4pi nu jnu -> erg s^-1 H^-1
    draine_data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)

    # convert from nu * P_nu to P_nu
    Hz = c/(draine_data['wave'])
    
    # extract and convert from erg/s/H to erg/s/cm^2/Hz/sr
    pah_pos = (draine_data['PAH^+'].value/(Hz * 4 * np.pi)) * N_HI
    pah_neut = (draine_data['PAH^0'].value/(Hz * 4 * np.pi)) * N_HI
    total = pah_pos + pah_neut
    
    # only use the MIR
    draine_wave = draine_data['wave']            
    cutoff = np.where( (draine_wave > wave_low) &  (draine_wave < wave_high))[0]
    draine_wave = draine_wave[cutoff]
    total = total[cutoff]
    
    # convert to Flambda 
    D21_angstrom = draine_data['wave'][cutoff].value*1e4
    D21_Flambda = total * c_A/(D21_angstrom)**2
    
    # fit continuum and subtract to isolate PAH features
    ind = np.ones(len(total), dtype = bool)
    ind = remove_PAH(pah_wave1, pah_wave2, ind, draine_wave)
    
    if filt_contam == 'F1500W':
        ind = remove_PAH(pah_wave3, pah_wave4, ind, draine_wave)
    
    # fit continuum
    D21_pfit = np.polyfit(D21_angstrom[ind], D21_Flambda[ind], deg = p_deg)
    D21_p = np.poly1d(D21_pfit)
    xx = np.linspace(min(D21_angstrom), max(D21_angstrom))
    
    # subtract continuum
    D21_pah = D21_Flambda - D21_p(D21_angstrom)
    
    
    # # multiply by R-J tail (lambda**2) to see if it helps
    # Flambda = Flambda * angstrom**2
    # D21_Flambda = D21_Flambda * D21_angstrom**2
    
    # normalize each by the PAH peak
    feat_11 = np.where(((angstrom) > (3.25 * u.micron)) & ((angstrom) < (3.35 * u.micron)))[0]
    Flambda_norm = Flambda/np.nanmax(Flambda[feat_11])
    
    D21_feat_11 = np.where(((D21_angstrom * u.AA) > (3.25 * u.micron)) & ((D21_angstrom * u.AA) < (3.35 * u.micron)))[0]
    D21_Flambda_norm = D21_Flambda/np.nanmax(D21_Flambda[D21_feat_11])
    
    norm_pah = pah/np.nanmax(pah[feat_11])
    norm_D21_pah = D21_pah/np.nanmax(D21_pah[D21_feat_11])


    # getting the synphot values
    synphot_Jy = get_synphot(angstrom.value, pah, filt_list)
    print(synphot_Jy)
    D21_synphot_Jy = get_synphot(D21_angstrom, D21_pah, filt_list)
    print(D21_synphot_Jy)
    


    ###########################################
    ######## Plotting to check con sub ########
    ###########################################
    
    plt.figure(2, figsize = (10, 5))
    plt.clf()
    plt.plot(angstrom, Flambda_norm, c = 'orange', lw = 2, label = 'PDRs4ALL', alpha = 0.6, zorder = 10)
    plt.plot(xx, p(xx)/np.nanmax(Flambda[feat_11]), ls = '--', c  = 'salmon', label = 'PDRs4All Continuum')
    plt.plot(D21_angstrom, D21_Flambda_norm, c = 'b', lw = 2, label = 'Draine+ 2021', alpha = 0.6, zorder = 10)
    plt.plot(xx, D21_p(xx)/np.nanmax(D21_Flambda[D21_feat_11]), ls = '--', c  = 'purple', label = 'Draine+ 2021 Continuum')

    plt.xlabel('Angstrom')
    plt.ylabel('Normalized F_lambda')
    plt.legend(loc = 'best')
    plt.title('No Continuum Subtraction')
    plt.ylim(0, 1.2)
    
    # plot transmission curves
    plt.plot(angstrom, 2*filt_dict[filt_list[0]](angstrom), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    plt.plot(angstrom, 2*filt_dict[filt_list[1]](angstrom), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    
    savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/{}_con_fit_tests_v4/'.format(filt_contam)
    savename = '{}_PDRs4ALL_{}_spec_confit_compare_D21_{:s}.pdf'.format(filt_contam, PDR_spec, method)
    plt.savefig(savepath + savename)
    
    plt.figure(3, figsize = (10, 5))
    plt.clf()
    plt.plot(angstrom, norm_pah, c = 'orange', lw = 2, label = 'PDRs4ALL', alpha = 0.6)
    plt.plot(D21_angstrom, norm_D21_pah, c = 'blue', lw = 2, label = 'Draine+ 2021', alpha = 0.6)
    plt.axhline(0, c = 'k', lw = 0.5)
    plt.ylim(-0.2, 1.2)
    
    # plot transmission curves
    plt.plot(angstrom, 2*filt_dict[filt_list[0]](angstrom), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    plt.plot(angstrom, 2*filt_dict[filt_list[1]](angstrom), c = 'gray', lw = 1, alpha = 0.5, ls = '-', zorder = 0.1)    
    
    # add synphot data
    ax = plt.gca()
    text_vals = 'k = {:s}/{:s} = {:2.2f} = {:3.1f}/{:3.1f}'.format(filt_list[0], filt_list[1], synphot_Jy[0]/synphot_Jy[1], synphot_Jy[0], synphot_Jy[1])
    plt.text(0.99, 0.95, s = text_vals, horizontalalignment='right', color = 'orange', transform=ax.transAxes)
    
    D21_text_vals = 'k = {:s}/{:s} = {:2.2f} = {:3.1f}/{:3.1f}'.format(filt_list[0], filt_list[1],D21_synphot_Jy[0]/D21_synphot_Jy[1], D21_synphot_Jy[0], D21_synphot_Jy[1])
    plt.text(0.99, 0.90, s = D21_text_vals, horizontalalignment='right', color = 'blue', transform=ax.transAxes)
    
    
    # plt.xlim(0, 15)
    plt.xlabel('Angstrom')
    plt.ylabel('F_lambda')
    plt.title('{:s} Method Continuum subtracted with deg = {} poly'.format(method, p_deg))
    savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/{}_con_fit_tests_v4/'.format(filt_contam)
    savename = '{}_PDRs4ALL_{}_spec_consub_compare_D21_{:s}.pdf'.format(filt_contam, PDR_spec, method)
    plt.savefig(savepath + savename)
    
    
                

    
# final_table = np.hstack((Flam_master, Jy_master))
# names = filt_list + filt_Jy

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/calc_PAH_contamination/'
# final_table_name = '{}_PDRs4ALL_spec_PAH_contamination_v4_test.txt'.format(filt_contam)
# ascii.write(final_table, savepath + final_table_name, names = names, overwrite = True)            
                

            
            
            




            
        
