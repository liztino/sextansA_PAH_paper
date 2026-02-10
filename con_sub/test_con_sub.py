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
import pysynphot as S
import stsynphot as stsyn
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import synphot
from specutils import Spectrum1D
from pahfit.model import Model
from astropy.nddata import StdDevUncertainty
import seaborn as sns
import matplotlib.lines as mlines


plt.ion()


# constants
c = 3e14        # in microns/sec
c_A = 3e18      # in Angstrom/sec
h = 6.62e-27    # in erg s 

N_HI = 1e21

# gets the filters 
filt_list = ['F560W', 'F770W', 'F1000W']
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

def split_con(con1, con2, ind):
    ind[np.where((draine_wave < con2) & ((draine_wave > con1)))[0]] = True
    return ind

# start arrays that we will populate
pah_math = []
pah_synphot = []


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

U_list = [1.0]
size_list = ['std']
ion_list = ['st']

for i in range(len(U_list)):
    for j in range(len(size_list)):
        for k in range(len(ion_list)):

    
            # # open a model from the D21 models
            # U = 1.0
            # size = 'std'
            # ion = 'st'
            
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
            
            # convert from nu * P_nu to P_nu
            Hz = c/(draine_data['wave'])
            
            # extract and convert from erg/s/H to erg/s/cm^2/Hz/sr# extract and convert from erg/s/H to erg/s/cm^2/Hz/sr
            pah_pos = (draine_data['PAH^+'].value/(Hz * 4 * np.pi)) * N_HI
            pah_neut = (draine_data['PAH^0'].value/(Hz * 4 * np.pi)) * N_HI
            total = pah_pos + pah_neut
            # total = (draine_data['total'].value/(Hz * 4 * np.pi)) * N_HI
            
            # only use the MIR
            draine_wave = draine_data['wave']            
            cutoff = np.where(draine_wave < 12)[0]
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
            
            # calculate the synthetic 
            syn_F560W = filt_through_Jy[0]
            syn_F770W = filt_through_Jy[1]
            syn_F1000W = filt_through_Jy[2]
            
            print('Filter values', filt_through_Jy)
            
            k = 5.564
            F560W_pivot = 5.635
            F770W_pivot = 7.639
            F1000W_pivot = 9.953
            
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
            
            pah_7_sub = fp2(syn_F560W, syn_F770W, syn_F1000W, F560W_pivot, F770W_pivot, F1000W_pivot, k)
            pah_7_sub_v2 = fp2_v2(syn_F560W, syn_F770W, syn_F1000W, F560W_pivot, F770W_pivot, F1000W_pivot, k)
            pah_7_sub_v3 = fp2_v3(syn_F560W, syn_F770W, syn_F1000W, F560W_pivot, F770W_pivot, F1000W_pivot, k)

            con_7 = fc2_v2(syn_F560W, syn_F770W, syn_F1000W, F560W_pivot, F770W_pivot, F1000W_pivot, k)
            
            # plot the synthetic spectrum
            plt.figure(1)
            plt.clf()
            plt.plot(draine_wave, total)
            
            ind = np.zeros(len(total), dtype = bool)
            ind = split_con(4,5, ind)
            ind = split_con(9.4, 10.2, ind)
            
            # fit continuum
            p_deg = 3
            pfit = np.polyfit(angstrom[ind], Flambda[ind], deg = p_deg)
            p = np.poly1d(pfit)
            xx = np.linspace(min(angstrom), max(angstrom))
            
            plt.figure(3)
            plt.clf()
            plt.plot(angstrom, Flambda)
            plt.plot(xx, p(xx), ls = '--', c  = 'purple')
            
            # subtract continuum
            pah = Flambda - p(angstrom)
            
            plt.figure(2)
            plt.clf()
            plt.plot(draine_wave, pah)
            plt.axhline(0, c = 'k', lw = 0.5)
            
            spec_ind = np.zeros(len(total), dtype = bool)
            b1 = 6.96
            b2 = 8.34
            spec_ind = split_con(b1, b2, spec_ind)
            
            plt.plot(draine_wave[spec_ind], pah[spec_ind])
            
                        
            ###########################################

            ######## Plotting to check con sub ########
            ###########################################
            plt.figure(1)
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

            plt.figure(2)
            plt.clf()
            plt.plot(angstrom, pah, c = 'b', lw = 2)
            plt.axhline(0, c = 'k', lw = 0.5)
            # plt.xlim(0, 15)
            plt.xlabel('Angstrom')
            plt.ylabel('F_lambda')
            plt.title('Continuum subtracted with deg = {} poly'.format(p_deg))

            
            # perform synthetic photometry on the continuum subtracted model spectrum 
            sp = synphot.SourceSpectrum(Empirical1D, points = angstrom, lookup_table=pah * units.FLAM, keep_neg = True)
            
            eff_stim = np.zeros(len(filt_list))
            eff_stim_Jy = np.zeros(len(filt_list))
            
            for f, filt_name in enumerate(filt_list):
                # create an observation and calculate the effective stimulation
                obs = synphot.Observation(sp, filt_dict[filt_name])
                eff_stim[f] = obs.effstim(flux_unit='flam').value
                
                eff_stim_Jy[f] = (eff_stim[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
                eff_stim_Jy[f] = eff_stim_Jy[f]/1e6
            
            model_k = eff_stim_Jy[1]/eff_stim_Jy[0]
            
            
            print('7.7 PAH orig v:', pah_7_sub)
            print('7.7 PAH v2:', pah_7_sub_v2)
            print('7.7 PAH v3:', pah_7_sub_v3)

            # print('7.7 PAH with filter math:', pah_7_sub)

            print('7.7 PAH calc with con sub:', eff_stim_Jy[1])
            print('Continuum at F770W', con_7)
            print('k assumed:', k)
            print('k from model', model_k)
            
            pah_math.append(pah_7_sub)
            pah_synphot.append(eff_stim_Jy[1])
            
            U_master.append(U)
            ion_master.append(ion)
            size_master.append(size)
            
            
            

# pah_math = np.array(pah_math)
# pah_synphot = np.array(pah_synphot)


# colors = sns.color_palette("hls", 3)

# # translate size text into actual size 
# d = {'sma': 20, 'lrg': 80, 'std': 40}
# size_arr = [d[x] for x in size_master] 

# # translate ionization into colors
# d = {'hi': colors[0], 'lo': colors[2], 'st': colors[1]}
# ion_arr = [d[x] for x in ion_master] 

# plt.figure(1)
# plt.clf()
# plt.scatter(pah_synphot, pah_math, s = size_arr, alpha = 0.7, c = ion_arr)
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

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/plots/new_method/'
# savename = 'test_newmethod_synphot_compare.pdf'
# plt.savefig(savepath + savename)

# perc = (pah_synphot -  pah_math)/pah_synphot

# plt.figure(2)
# plt.clf()
# plt.scatter(U_master, perc, s = size_arr, alpha = 0.7, c = ion_arr)
# plt.axhline(0, c = 'k', lw = 0.5)
# plt.xlabel('Ionization Parameter')
# plt.ylabel('Percent Difference')

# plt.legend(handles=[ hi, st, lo], loc = 'best')

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/plots/new_method/'
# savename = 'test_newmethod_synphot_percdiff.pdf'
# plt.savefig(savepath + savename)

