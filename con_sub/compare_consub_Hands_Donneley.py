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


# values for our filters
F560W = 0.36
F770W = 1.05
F1000W = 0.97
F1130W = 2.20
F1500W = 1.67

# get wavelengths
F560W_wave = get_pivot_wave('F560W')
F770W_wave = get_pivot_wave('F770W')
F1000W_wave = get_pivot_wave('F1000W')
F1130W_wave = get_pivot_wave('F1130W')
F1500W_wave = get_pivot_wave('F1500W')

# run through k method
pah_7 = k_eq.get_pah_low(F560W, F770W, F1000W, F560W_wave, F770W_wave, F1000W_wave, 4.33)
pah_11 = k_eq.get_pah_up(F1000W, F1130W, F1500W, F1000W_wave, F1130W_wave, F1500W_wave,  k=7.21)

# now test Grant's
# for 7.7
alpha_7 = 0.53
g_7 = 0.91
con_7 = (g_7*F560W**(1-alpha_7) * F1000W**(alpha_7))

print('My method 7 con:', pah_7['con'])
print('Grants method 7 con:', con_7)

# for 11.3
alpha_11 = 0.31
g_11 = 1
con_11 = (g_11 * F1000W**(1-alpha_11) * F1500W**(alpha_11))

print('My method 11 con:', pah_11['con'])
print('Grants method 11 con:', con_11)



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




