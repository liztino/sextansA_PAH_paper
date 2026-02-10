#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Isolates the PAH contribution to each filter with D21 models

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

erg_MJy = 1e-23

# calculate the column density 
N_HI = (1.105e21 * HI_flux)/(bmaj * bmin)

def split_con(con1, con2, ind):
    ind[np.where((draine_wave < con2) & ((draine_wave > con1)))[0]] = True
    return ind

# filters to apply synthetic photometry
# use only the pairs with contamination

filt_list = ['F335M', 'F360M', 'F560W', 'F770W']
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


U_list = np.arange(0, 7, 0.5)
ion_list = ['lo', 'st', 'hi']
size_list = ['sma', 'std', 'lrg']

# U_list = [7,3]
# ion_list = ['hi','lo']
# size_list = ['lrg','sma']

U_master = []
ion_master = []
size_master = []

Flam_master =  []
Jy_master = []

for i in range(len(U_list)):
    for j in range(len(size_list)):
        for k in range(len(ion_list)):
            
            U = U_list[i]
            size = size_list[j]
            ion = ion_list[k]
            
            # load Draine model
            drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
            model = 'BC03_Z0.0004_10Myr'
            # drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
            drainename = 'pahspec.out_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}'.format(U, ion, size)
            header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']
            
            # D21 models are given in 4pi nu jnu -> erg s^-1 H^-1
            draine_data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)

            ind = np.where(draine_data['wave'] > 40)[0][0]
            
            print(ind)
            print('wave end', draine_data['wave'][ind])
            
            # convert to useable units 
            micron = draine_data['wave'][:ind]
            Hz = c/(micron)
            total = ((draine_data['total'][:ind]/(Hz * 4 * np.pi)) * N_HI)/erg_MJy
            
            # prepare spectrum for pahfit
            unc = StdDevUncertainty(0.01 * total)
            # flux = NDData(data = total, uncertainty = unc, unit = u.Jy / u.sr )
            flux_q = Quantity(total * u.Jy / u.sr )
            wave = Quantity(micron * u.micron)
            fit_spec = Spectrum1D(flux = flux_q, spectral_axis = wave, redshift = 0, uncertainty = unc)

            # apply pahfit
            fit_spec.meta['instrument'] = 'draine21*'
            model = Model.from_yaml('classic_nir.yaml')
            model.guess(fit_spec)
            model.fit(fit_spec, redshift = 0, verbose = True)     
            
            tot = model.tabulate(wavelengths=wave, instrumentname='draine21*', redshift=0)
            pah = model.tabulate(wavelengths=wave, instrumentname='draine21*', redshift=0, feature_mask = model.features['kind']=='dust_feature') 
            con = model.tabulate(wavelengths=wave, instrumentname='draine21*', redshift=0, feature_mask = model.features['kind']=='dust_continuum')

        
            ###########################################
            ######## Plotting to check con sub ########
            ###########################################
            plt.figure(1)
            plt.clf()
            plt.plot(micron, Hz * total, c = 'k', lw =2, label = 'Original')
            plt.plot(tot.wavelength.to(u.micron), tot.frequency.to(u.Hz) * tot.flux, c = 'b', lw =1.5, label = 'PAHFIT', alpha = 0.7)
            plt.plot(pah.wavelength.to(u.micron), pah.frequency.to(u.Hz) * pah.flux, c = 'orange', lw =1.5, label = 'PAH', alpha = 0.7)
            plt.plot(con.wavelength.to(u.micron), con.frequency.to(u.Hz) * con.flux, c = 'purple', lw =1.5, label = 'Continuum', alpha = 0.7)

            # plt.plot(draine_wave[ind], total[ind], '.', c = 'orange', label = 'Continuum indices')
            plt.xlim(1,20)
            # # plt.ylim(1e-18,1e-13)
            # plt.semilogy()
            plt.xlabel('Wavelength (Micron)')
            plt.ylabel('$\\nu$ F$_{\\nu}$')
            plt.legend(loc = 'best')
            savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/pahfit_continuum/'
            savename = 'pahfit_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}.pdf'.format(U, ion, size)
            plt.savefig(savepath + savename)

            
            ###########################################
            ######## Synthetic photometry #############
            ###########################################
            
            
            # convert to angstroms and F_lambda units
            tot_flux_lam = (tot.flux).to(u.erg / u.s / u.AA / u.cm**2 / u.sr, equivalencies=u.spectral_density(tot.wavelength))
            pah_flux_lam = (pah.flux).to(u.erg / u.s / u.AA / u.cm**2 / u.sr, equivalencies=u.spectral_density(pah.wavelength))
            con_flux_lam = (con.flux).to(u.erg / u.s / u.AA / u.cm**2 / u.sr, equivalencies=u.spectral_density(con.wavelength))

            tot_lam = (tot.wavelength).to(u.AA, equivalencies=u.spectral_density(tot.wavelength))
            pah_lam = (pah.wavelength).to(u.AA, equivalencies=u.spectral_density(pah.wavelength))
            con_lam = (con.wavelength).to(u.AA, equivalencies=u.spectral_density(con.wavelength))
            
            # turn the total spectrum into a source spectrum object in synphot 
            sp = synphot.SourceSpectrum(Empirical1D, points = pah_lam, lookup_table=pah_flux_lam.value * units.FLAM, keep_neg = True)
            
            # loop through all filters
            filt_through = np.zeros(len(filt_list))
            filt_through_Jy = np.zeros(len(filt_list))

            for f, filt_name in enumerate(filt_list):
                # create an observation and calculate the effective stimulation
                obs = synphot.Observation(sp, filt_dict[filt_name])
                filt_through[f] = obs.effstim(flux_unit='flam').value
                
                filt_through_Jy[f] = (filt_through[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
                filt_through_Jy[f] = filt_through_Jy[f]/1e6
                
            # save the filter information
            Flam_master.append(filt_through)
            Jy_master.append(filt_through_Jy)

            U_master.append(U)
            ion_master.append(ion)
            size_master.append(size)
            
param_table = np.array([U_master, ion_master, size_master]).T
final_table = np.hstack((param_table, Flam_master, Jy_master))

name1 = ['U', 'ion', 'size']

filt_Jy = []
for filt in filt_list:
    filt_Jy.append(filt + '_Jy')
    
names = name1 + filt_list + filt_Jy

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/calc_PAH_contamination/'
final_table_name = 'pah_consub_calc_contamination_pahfit.txt'
ascii.write(final_table, savepath + final_table_name, names = names, overwrite = True)            
                

            
            
            




            
        
