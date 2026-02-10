#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculates the contribution of the hot dust in the F770W and F2100W filters from D21 models

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import scipy.interpolate as inter
from astropy.table import Table

os.environ['PYSYN_CDBS']
import pysynphot as S
import stsynphot as stsyn
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import synphot

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

# filters to apply synthetic photometry
# use only the pairs with contamination

filt_list = ['F770W', 'F2100W']

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W', 'F2100W'}
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

wave_low = 4
wave_high = 30

U_master = []
ion_master = []
size_master = []

Flam_master =  []
Jy_master = []

plt.figure(2)
plt.clf()
fig, axs = plt.subplots(nrows = 4, ncols = 4, num =2)
axs_list = axs.ravel()

q_pah_list = np.linspace(0, 10, 200)

for i in range(len(q_pah_list)):
            
    # U = U_list[i]
    U = 0
    ion = 'st'
    size = 'std'
    
    # load Draine model
    drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
    model = 'BC03_Z0.0004_10Myr'
    # drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
    drainename = 'pahspec.out_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}'.format(U, ion, size)
    header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']
    
    # D21 models are given in 4pi nu jnu -> erg s^-1 H^-1
    draine_data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)
    
    qpah = q_pah_list[i]
    # qpah_text = '0p2'
    orig_qpah = 3.8
    factor = qpah/orig_qpah
    

    # convert from nu * P_nu to P_nu
    Hz = c/(draine_data['wave'])
    
    # astrodust = (draine_data['total'].value/(Hz * 4 * np.pi)) * N_HI
    
    # extract and convert from erg/s/H to erg/s/cm^2/Hz/sr
    astrodust = (draine_data['Astrodust'].value/(Hz * 4 * np.pi)) * N_HI
    pah_pos = (draine_data['PAH^+'].value/(Hz * 4 * np.pi)) * N_HI
    pah_neut = (draine_data['PAH^0'].value/(Hz * 4 * np.pi)) * N_HI
    total = (draine_data['total'].value/(Hz * 4 * np.pi)) * N_HI
    
    # convert to the given qpah
    new_pah = factor * pah_pos + factor * pah_neut
    new_tot = new_pah + astrodust
    
    
    # only use the MIR
    draine_wave = draine_data['wave']            
    cutoff = np.where( (draine_wave > wave_low) &  (draine_wave < wave_high))[0]
    draine_wave = draine_wave[cutoff]
    total = new_tot[cutoff]
    
    # convert to Flambda 
    angstrom = draine_data['wave'][cutoff].value*1e4
    Flambda = total * c_A/(angstrom)**2

    # axs_list[i].plot(draine_wave, total, c = 'k', lw = 2)
    # axs_list[i].text(0.95, 0.05, 'U = {:2.1f}'.format(U), transform = axs_list[i].transAxes, ha = 'right')

    
    ###########################################
    ######## Synthetic photometry #############
    ###########################################
    
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
        
    # save the filter information
    Flam_master.append(filt_through)
    Jy_master.append(filt_through_Jy)
    
plt.tight_layout()
fig.delaxes(axs_list[15])
fig.delaxes(axs_list[14])  

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/hot_dust_D21/figures/'
# savename = 'D21_hot_dust_F770W_F2100W_ratio_spec_withpah_{:s}.pdf'.format(qpah_text)
# plt.savefig(savepath + savename)


Jy_master = np.array(Jy_master)
# ratio = Jy_master[:,0]/Jy_master[:,1]

# plt.figure(1)
# plt.clf()
# plt.scatter(U_list, np.log10(ratio), s = 20)
# plt.xlabel('Ionization Parameter, logU')
# plt.ylabel('Log(F770W/F2100W)')
# plt.title('qPAH = {:2.1f}'.format(qpah))
# # plt.semilogy()

# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/hot_dust_D21/figures/'
# savename = 'D21_hot_dust_F770W_F2100W_ratio_withpah_{:s}.pdf'.format(qpah_text)
# plt.savefig(savepath + savename)

# save the text file 
savedata = Table(Jy_master, names = ['F770W', 'F2100W'])
savedata.add_column(q_pah_list, name = 'q_pah')
reorder = ['q_pah', 'F770W', 'F2100W']
savedata = savedata[reorder]

filename = 'D21_ascii_hot_dust_F770W_F2100W_ratio_vary_q_pah.txt'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/hot_dust_D21/'

ascii.write(savedata, savepath + filename, overwrite = True)


            
            
            