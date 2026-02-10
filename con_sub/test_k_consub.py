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

plt.ion()


# constants
c = 3e14        # in microns/sec
c_A = 3e18      # in Angstrom/sec
h = 6.62e-27    # in erg s 

N_HI = 1e21


def remove_PAH(pah1, pah2, ind, wave):
    ind[np.where((wave > pah1) & ((wave < pah2)))[0]] = False
    return ind

# filt_contam = 'F335M'
# k_val = 20

filt_contam = 'F770W'
k_val =5.87

# filt_contam = 'F1500W'
# k_val = 10.17


if filt_contam == 'F770W':
    wave_low = 4.5    
    wave_high = 10
    
    pah_wave1 = 5
    pah_wave2 = 9.4
    
    filt_list = ['F560W', 'F770W', 'F1000W']
    filt_list2 = ['F770W', 'F560W']
    
    p_deg = 1
    
elif filt_contam == 'F335M':
    wave_low = 2
    wave_high = 5
    
    pah_wave1 = 3
    pah_wave2 = 4
    
    filt_list = ['F300M', 'F335M', 'F360M']
    filt_list2 = ['F335M', 'F360M']
    
    p_deg = 1
    
elif filt_contam == 'F1500W':
    wave_low = 10
    wave_high = 20
    
    pah_wave1 = 10.4
    pah_wave2 = 14.8
    
    pah_wave3 = 15.8
    pah_wave4 = 18.2
    
    filt_list = ['F1000W', 'F1130W', 'F1500W']
    filt_list2 = ['F1130W', 'F1500W']
    
    p_deg = 6
    
    
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

D21_sp = 'BC03_Z0.02_100Myr'
time = '1e8'
z = '0.02'

# U_list = [7,3]
# ion_list = ['hi','lo']
# size_list = ['lrg','sma']

U_master = []
ion_master = []
size_master = []

Flam_master =  []
Jy_master = []

# U_list = [3.0]
# size_list = ['std']
# ion_list = ['st']

for i in range(len(U_list)):
    for j in range(len(size_list)):
        for k in range(len(ion_list)):
            
            U = U_list[i]
            size = size_list[j]
            ion = ion_list[k]
            
            ### step 1 - calculate with the k method 
            
            # load Draine model
            drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
            # drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
            drainename = 'pahspec.out_bc03_z{:s}_{:s}_{:01.2f}_{:s}_{:s}.gz'.format(z, time, U, ion, size)
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
            
            # only use the MIR
            draine_wave = draine_data['wave']            
            cutoff = np.where( (draine_wave < 20))[0]
            draine_wave = draine_wave[cutoff]
            total = total[cutoff]
            
            # convert to Flambda 
            angstrom = draine_data['wave'][cutoff].value*1e4
            Flambda = total * c_A/(angstrom)**2
        
            # turn the total spectrum into a source spectrum object in synphot 
            sp = synphot.SourceSpectrum(Empirical1D, points = angstrom, lookup_table=Flambda * units.FLAM, keep_neg = True)
            
            # loop through all filters
            filt_through = np.zeros(len(filt_list))
            filt_through_Jy = np.zeros(len(filt_list))

            for f, filt_name in enumerate(filt_list):
                # create an observation and calculate the effective stimulation
                obs = synphot.Observation(sp, filt_dict[filt_name])
                filt_through[f] = obs.effstim(flux_unit='flam').value
                
                filt_through_Jy[f] = (filt_through[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
                filt_through_Jy[f] = filt_through_Jy[f]/1e6
                
            
            f1 = filt_through_Jy[0]
            f2 = filt_through_Jy[1]
            f3 = filt_through_Jy[2]
            
            lam1 = get_pivot_wave(filt_list[0])
            lam2 = get_pivot_wave(filt_list[1])
            lam3 = get_pivot_wave(filt_list[2])
            
            def fp2_v2(f1, f2, f3, lam1, lam2, lam3, k):
                # equation calculated through Julia's math
                fp1 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(1 - k - ((lam2 - lam1)/(lam3 - lam1)))
                fp2 = fp1*k
                
                print('fp1', fp1)
                print('fp2', fp2)
                
                return fp2
            
            def fc2_v2(f1, f2, f3, lam1, lam2, lam3, k):
                # equation calculated through Julia's math
                fp1 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(1 - k - ((lam2 - lam1)/(lam3 - lam1)))
                fc1 = f1 - fp1    
                fc2 = (((f3 - fc1)/(lam3 - lam1)) * (lam2 - lam1)) + fc1
                
                print('fc1', fc1)
                print('fc2', fc2)
                print('slope', (((f3 - fc1)/(lam3 - lam1))))
                # print('delta_lam', lam2 - lam1)
                
                return fc2
            
            
            def fp2(f1, f2, f3, lam1, lam2, lam3, k):
                # equation calculated through Julia's math
                fp1 = (f3 * ((lam2 - lam1)/(lam3-lam1)) + f1 * ((lam3 - lam2)/(lam3 - lam1)) - f2)/(((lam3 - lam2)/(lam3 - lam1)) - k)
                fp2 = fp1*k
                
                print('fp1', fp1)
                print('fp2', fp2)
                
                return fp2
            
            def fc2(f1, f2, f3, lam1, lam2, lam3, k):
                # equation calculated through Julia's math
                fp1 = (f3 * ((lam2 - lam1)/(lam3-lam1)) + f1 * ((lam3 - lam2)/(lam3 - lam1)) - f2)/(((lam3 - lam2)/(lam3 - lam1)) - k)
                fc1 = f1 - fp1    
                fc2 = (((f3 - fc1)/(lam3 - lam1)) * (lam2 - lam1)) + fc1
                
                print('fc1', fc1)
                print('fc2', fc2)
                
                return fc2
            
            def fp2_v3(f1, f2, f3, lam1, lam2, lam3, k):
                # equation calculated through Julia's math
                fp3 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(((lam2 - lam1)/(lam3 - lam1)) - k)
                fp2 = fp3*k
                
                print('fp3', fp3)
                print('fp2', fp2)
                
                return fp2
            
            # calculate the synthetic 
            syn_F560W = filt_through_Jy[0]
            syn_F770W = filt_through_Jy[1]
            syn_F1000W = filt_through_Jy[2]
            
            k = k_val
            
            # pah_7_sub = fp2(syn_F560W, syn_F770W, syn_F1000W, lam1, lam2, lam3, k)

            # con_7 = fc2_v2(syn_F560W, syn_F770W, syn_F1000W, F560W_pivot, F770W_pivot, F1000W_pivot, k)

            if filt_contam == 'F770W':
                print(f1, f2, f3)
                print(lam1, lam2, lam3)
                consub = k_eq.get_pah_low(f1, f2, f3, lam1, lam2, lam3, k_val)
            else:
                consub = k_eq.get_pah_up(f1, f2, f3, lam1, lam2, lam3, k_val)
                
            plt.figure(3)
            plt.clf()
            plt.plot(angstrom, Flambda)
            
                
            #### step 2 - calculate photometry with continuum subtraction
            
            # now we will use just the PAH portion of the spectrum 
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
            angstrom = draine_data['wave'][cutoff].value*1e4
            Flambda = total * c_A/(angstrom)**2
            
            # fit continuum and subtract to isolate PAH features
            ind = np.ones(len(total), dtype = bool)
            ind = remove_PAH(pah_wave1, pah_wave2, ind, draine_wave)
            
            if filt_contam == 'F1500W':
                ind = remove_PAH(pah_wave3, pah_wave4, ind, draine_wave)
            
            # convert to Flambda before consub 
            angstrom = draine_wave.value*1e4
            Flambda = total * c_A/(angstrom)**2
            
            # fit continuum
            pfit = np.polyfit(angstrom[ind], Flambda[ind], deg = p_deg)
            p = np.poly1d(pfit)
            xx = np.linspace(min(angstrom), max(angstrom))
            
            # subtract continuum
            pah = Flambda - p(angstrom)
            
            # perform synthetic photometry on the continuum subtracted model spectrum 
            sp = synphot.SourceSpectrum(Empirical1D, points = angstrom, lookup_table=pah * units.FLAM, keep_neg = True)
            
            eff_stim = np.zeros(len(filt_list))
            eff_stim_Jy = np.zeros(len(filt_list))
            
            for f, filt_name in enumerate(filt_list):
                # create an observation and calculate the effective stimulation
                obs = synphot.Observation(sp, filt_dict[filt_name], force='extrap')
                eff_stim[f] = obs.effstim(flux_unit='flam').value
                
                eff_stim_Jy[f] = (eff_stim[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
                eff_stim_Jy[f] = eff_stim_Jy[f]/1e6
                 
            
            if filt_contam == 'F770W':
                D21_pah = eff_stim_Jy[1]
                D21_contam = eff_stim_Jy[0]
                model_k = D21_pah/D21_contam
            else:
                D21_pah = eff_stim_Jy[1]
                D21_contam = eff_stim_Jy[2]
                model_k = D21_pah/D21_contam
                        
            # print('PAH from k method:', pah_7_sub)
            # print('CON from k method', consub['con'])
            # print('Slope from k method', consub['slope'])


            # print('PAH calc with con sub:', D21_pah)
            # print('k assumed:', k_val)
            # print('k from model', model_k)
            print(consub)
            
            pah_k_method.append(consub['pah'])
            pah_synphot.append(D21_pah)
            contam_synphot.append(D21_contam)
            
            U_master.append(U)
            ion_master.append(ion)
            size_master.append(size)
            
            f1_array.append(f1)
            f2_array.append(f2)
            f3_array.append(f3)
            
            ###########################################

            ######## Plotting to check con sub ########
            ###########################################
            plt.figure(4)
            plt.clf()
            plt.plot(angstrom, Flambda, c = 'k', lw = 2, label = 'D21 PAH')
            plt.plot(angstrom[ind], Flambda[ind], '.', c = 'orange', label = 'Continuum indices')
            plt.plot(xx, p(xx), ls = '--', c  = 'purple', label = 'Continuum Fit')
            # plt.xlim(0, 15)
            # # plt.ylim(1e-18,1e-13)
            # plt.semilogy()
            plt.xlabel('Angstrom')
            plt.ylabel('F_lambda')
            plt.legend(loc = 'best')
            savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/{}_con_fit_tests/'.format(filt_contam)
            savename = 'Flam_pahspec_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}.pdf'.format(U, ion, size)
            #plt.savefig(savepath + savename)

            plt.figure(5)
            plt.clf()
            plt.plot(angstrom, pah, c = 'b', lw = 2)
            plt.axhline(0, c = 'k', lw = 0.5)
            # plt.xlim(0, 15)
            plt.xlabel('Angstrom')
            plt.ylabel('F_lambda')
            plt.title('Continuum subtracted with deg = {} poly'.format(p_deg))
            savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/{}_con_fit_tests/'.format(filt_contam)
            savename = 'Flam_consub_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}.pdf'.format(U, ion, size)
            #plt.savefig(savepath + savename)
            
            
            

# save the files
pah_k_method = np.array(pah_k_method)
pah_synphot = np.array(pah_synphot)
f1_array = np.array(f1_array)
f2_array = np.array(f2_array)
f3_array = np.array(f3_array)


final_table = np.array([U_master, ion_master, size_master, pah_k_method, pah_synphot, f1_array, f2_array, f3_array]).T
# final_table = np.hstack((param_table, pah_k_method, pah_synphot, f1_array, f2_array, f3_array))

name1 = ['U', 'ion', 'size']
name2 = ['pah_k_method', 'pah_synphot']

# filt_Jy = []
# for filt in filt_list:
#     filt_Jy.append(filt + '_Jy')
    
names = name1 + name2 + filt_list 

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/D21_models/'
final_table_name = D21_sp +  '_' + filt_contam + '_D21_models_with_k_consub_none.txt'
# ascii.write(final_table, savepath + final_table_name, names = names, overwrite = True)            
                



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

