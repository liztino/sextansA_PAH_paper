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

def fitFunc(xy,a,b,c,d):
	'''
	functional form that takes in x = F1130W/F1000W and y = F770W/F1000W
	returns z = percentage of filter from PAHs
	'''
	x,y = xy
	return a + b*x + c*y + d*x*y

def colorcolorFunc770(pop,x,y):
	'''
	returns the percentage of the F770W filter from PAHs
	x = F1130W/F1000W
	y = F770W/F1000W
	'''
	a = pop[0] 
	b = pop[1] 
	c = pop[2] 
	d = pop[3] 
	

	return fitFunc([x,y], a,b,c,d)
	
def colorcolorFunc1130(pop,x,y):
	'''
	returns the percentage of the F1130W filter from PAHs
	x = F1130W/F1000W
	y = F770W/F1000W
	'''
	a = pop[0] 
	b = pop[1] 
	c = pop[2] 
	d = pop[3] 
	

	return fitFunc([x,y], a,b,c,d)


coeff_1130 = [-14.1122251,32.7828093,20.1457013,-7.84496217]
coeff_770 = [-68.3960087,37.7059632,46.5365895,-11.0003656]


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
            
            # apply Lindsay's subtraction
            ### plug MJy/sr values here
            here = [f3/f2,f1/f2] 
            
            ### calculate PAH percentage from function form 
            pahperc_770 = colorcolorFunc770(coeff_770,here[0],here[1])
            pahperc_1130 = colorcolorFunc1130(coeff_1130,here[0],here[1])
            
            ### make sure the percentages don't go out of physical bounds (this is optional but I argue for it)
            ### calculate F770W and F1130W PAH emission 
            if pahperc_770>100 or pahperc_770<0:
                F770W_PAH = np.nan
            else:
                F770W_PAH = pahperc_770/100*f1
                
            if pahperc_1130>100 or pahperc_1130<0:
                F1130W_PAH = np.nan
            else:
                F1130W_PAH = pahperc_1130/100*f3
            
                
            # plt.figure()
            # plt.plot(angstrom, Flambda)
            # plt.title('pah {:5.4f} con {:5.4f}'.format(consub['pah'], consub['con']))


            pah_k_method.append(F1130W_PAH)
            
            print(F770W_PAH)
            print(F1130W_PAH)
            
            U_master.append(U)
            ion_master.append(ion)
            size_master.append(size)
            
            f1_array.append(f1)
            f2_array.append(f2)
            f3_array.append(f3)
            

            
            

# save the files
pah_k_method = np.array(pah_k_method)
pah_synphot = np.array(pah_synphot)
f1_array = np.array(f1_array)
f2_array = np.array(f2_array)
f3_array = np.array(f3_array)


final_table = np.array([U_master, ion_master, size_master, pah_k_method,  f1_array, f2_array, f3_array]).T
# final_table = np.hstack((param_table, pah_k_method, pah_synphot, f1_array, f2_array, f3_array))

name1 = ['U', 'ion', 'size']
name2 = ['pah_k_method']

# filt_Jy = []
# for filt in filt_list:
#     filt_Jy.append(filt + '_Jy')
    
names = name1 + name2 + filt_list 

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/D21_models/'
final_table_name = D21_sp +  '_' + filt_pah + '_D21_models_Lindsay_Hands.txt'
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

