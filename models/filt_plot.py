#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Makes a quick plot of JWST Filters on D21 Model

@author: etarantino
"""


from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import scipy.interpolate as inter
# import seaborn as sns

# os.environ['PYSYN_CDBS']
# import pysynphot as S
# import stsynphot as stsyn
# from synphot import SourceSpectrum, SpectralElement, units
# from synphot.models import Empirical1D
# import synphot

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

# calculate the column density 
N_HI = (1.105e21 * HI_flux)/(bmaj * bmin)

filt_list = [ 'F560W', 'F770W', 'F1000W']
miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
filt_dict = {}
filt_wave = np.zeros(len(filt_list))



def get_filt(filt_name):
    # load filter info
    if filt_name in nircam_filts:
        filts_dir = '/Users/etarantino/Documents/JWST/filts/nircam/'
    else:
        filts_dir = '/Users/etarantino/Documents/JWST/filts/miri/'
        
    filt_file = f'{filts_dir}/{filt_name}_mean_system_throughput.txt'
    filt = ascii.read(filt_file)
    
    return filt
    


U = 3
ion = 'lo'
size = 'sma'
qpah = 0.5

# load Draine model
drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
model = 'BC03_Z0.0004_10Myr'
# drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
drainename = 'pahspec.out_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}'.format(U, ion, size)
header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']

draine_data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)


orig_qpah = 3.51
factor = qpah/orig_qpah

# convert from nu * P_nu to P_nu
Hz = c/(draine_data['wave'])

# extract and convert from erg/s/H to erg/s/cm^2/Hz/sr
astrodust = (draine_data['Astrodust'].value/(Hz * 4 * np.pi)) * N_HI
pah_pos = (draine_data['PAH^+'].value/(Hz * 4 * np.pi)) * N_HI
pah_neut = (draine_data['PAH^0'].value/(Hz * 4 * np.pi)) * N_HI
total = (draine_data['total'].value/(Hz * 4 * np.pi)) * N_HI
    

# only use the MIR
draine_wave = draine_data['wave']            
cutoff = np.where(draine_wave < 15)[0]
draine_wave = draine_wave[cutoff]
total = total[cutoff]


plt.figure(1)
plt.clf()
plt.plot(draine_wave, total, c = 'k', lw = 3)
plt.xlim(4,13)
plt.ylim(0, 5e-14)

ax = plt.gca()

F560W = get_filt('F560W')
y_coord = np.zeros(len(F560W['Throughput']))
ax.fill_between(F560W['Microns'], F560W['Throughput']*1.5*np.max(total), y_coord, alpha = 0.3, color = 'magenta', zorder = 0.0001)

F770W = get_filt('F770W')
y_coord = np.zeros(len(F770W['Throughput']))
ax.fill_between(F770W['Microns'], F770W['Throughput']*1.5*np.max(total), y_coord, alpha = 0.3, color = 'cyan', zorder = 0.0001)

F1000W = get_filt('F1000W')
y_coord = np.zeros(len(F1000W['Throughput']))
ax.fill_between(F1000W['Microns'], F1000W['Throughput']*1.5*np.max(total), y_coord, alpha = 0.3, color = 'yellow', zorder = 0.0001)

plt.xlabel('Wavelength (Microns)', size = 'large')
plt.ylabel('Flux', size = 'large')

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/'
savename = 'PAH_contamination_filter.pdf'
plt.savefig(savepath + savename)


