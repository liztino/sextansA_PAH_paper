#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script to make figure for paper that shows the contamination on the filters with D21 models and PDRs4All

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter
from astropy.table import Table
import os, glob
import seaborn as sns
import astropy.units as u
from matplotlib import ticker

os.environ['PYSYN_CDBS']
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


filt_list = ['F300M', 'F335M', 'F360M','F560W',  'F770W', 'F1000W', 'F1130W', 'F1500W']

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W', 'F2100W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
filt_dict = {}
filt_wave = np.zeros(len(filt_list))

# load the filter transmission curves and info 
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
    
# load in a D21 model
drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
model = 'BC03_Z0.0004_10Myr'

U = 1
ion = 'st'
size = 'std'
drainename = 'pahspec.out_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}'.format(U, ion, size)
header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']

# D21 models are given in 4pi nu jnu -> erg s^-1 H^-1
draine_data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)

# convert from nu * P_nu to P_nu
Hz = c/(draine_data['wave'])

# extract and convert from erg/s/H to erg/s/cm^2/Hz/sr
total = (draine_data['total'].value/(Hz * 4 * np.pi)) * N_HI
# convert to MJy to be in line with JWST units
total = (total * u.erg / u.s / u.cm**2 / u.Hz).to(u.Jy)/1e6

# load a PDRs4All spectrum 
pdrs_filepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/PDRs4ALL_data/orion_template_fits_and_clean_specs/'
spec_val = 'DF1'
orig_spec = f'spec_orion_{spec_val}'
pdrs_data = ascii.read(pdrs_filepath + orig_spec + '.ecsv')

PDR_spec_arr = ['HII_flux', 'Atomic_flux', 'DF1_flux', 'DF2_flux', 'DF3_flux']


###########################
####### Now plot ##########
###########################

plt.figure(1, figsize = (8.5, 11))
plt.clf()
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True, num = 1)

plt.subplots_adjust(hspace = 0)
xlim = (2.5, 19)
ylim = (0, 45)

ax1.plot(draine_data['wave'], total * c/draine_data['wave'] * 1e-15 , c ='k', lw = 1.5, alpha =0.8)
ax1.set_xlim(xlim)
# ax1.set_ylim(ylim)
ax1.semilogx()

props = dict(boxstyle='round', facecolor='white', alpha=1)

ax1.text(0.015, 0.9, 'Draine+ 2021', transform = ax1.transAxes, size = 16, bbox = props)


ax2.plot(pdrs_data['wavelength'], pdrs_data['flux'] * c/pdrs_data['wavelength'], c ='k', lw = 1.5,  alpha =0.8)
ax2.set_xlim(xlim)
ax2.set_ylim(-2e16, 5.75e17)
ax2.semilogx()

ax2.text(0.015, 0.9, 'PDRs4All', transform = ax2.transAxes, size = 16, bbox = props)

ax2.get_yaxis().get_offset_text().set_visible(False)



colors = sns.color_palette("hls", len(filt_list))[::-1]


for f, filt_name in enumerate(filt_list):
    if filt_name in nircam_filts:
        filts_dir = '/Users/etarantino/Documents/JWST/filts/nircam/'
    else:
        filts_dir = '/Users/etarantino/Documents/JWST/filts/miri/'
        
    filt_file = f'{filts_dir}/{filt_name}_mean_system_throughput.txt'
    th = ascii.read(filt_file)
    
    cutoff = np.nanmax(th['Throughput'])*0.1
    ind = np.where(th['Throughput'] >  cutoff)[0]
    ind1 = ind[0]
    ind2 = ind[-1]
    mic_low = th['Microns'][ind1]
    mic_high = th['Microns'][ind2]
    
    print(filt_name, mic_low, mic_high)
    
    # ax1.scatter(phot[row][f'{filt}_pivot'].value, phot[row][f'{filt}_jy'].value *1e3, s = 200, marker = 'o', color = colors[f], zorder = 10000, alpha = 0.7, edgecolor = 'k', linewidth = 1.5, label = 'MIRI Photometry')
    # ax1.errorbar(phot[row][f'{filt}_pivot'].value, phot[row][f'{filt}_jy'].value *1e3, yerr = phot[row][f'{filt}_err_jy'].value *1e3, 
    #               fmt = '.', c = 'k', alpha = 0.5, elinewidth = 0.7, capsize = 10 )
    # print(phot[row][f'{filt}_jy'].value *1e3)
    # print(phot[row][f'{filt}_err_jy'].value *1e3)
    
        
    pah_filts = {'F335M', 'F770W', 'F1130W'}
    con_filts = {'F300M', 'F1000W'}
    pah_con_filts = {'F360M', 'F560M', 'F1500W'}
        
    
    if filt_name in pah_filts:
        # gold
        color = '#009E73'
        # color = '#CC79A7'
        # color = 'limegreen'
    elif filt_name in con_filts:
        # teal
        color = '#CC79A7'
        # color = 'darkorchid'
    elif filt_name in pah_con_filts:
        # purple
        color = '#E69F00'
        color = '#F7B67C'
        # color = 'darkorchid'
    
    ax3.plot(th['Microns'], th['Throughput'], lw = 1.5, color = color)
    ax3.fill_between(th['Microns'], th['Throughput'], np.zeros(len(th['Microns'])), color = color, alpha = 0.5)
    ax3.set_ylim(0, 0.65)
    

    y_label_pos = np.max(th['Throughput']) + 0.025
    ax3.text((mic_high - mic_low)/2 + mic_low, y_label_pos, filt_name, c = color, va = 'center', ha = 'center', fontweight = 'bold')
    ax3.semilogx()
    
def calc_syn_phot(wave, flux, filt_list):
    # conversions
    # convert microns to angstroms
    angstrom = (wave).to(u.AA)
    
    # convert MJy to F_lambda
    Jy = (flux.value * 1e6) * u.Jy
    Flambda = Jy.to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(angstrom))
    Flambda = Flambda.value
    
    ###########################################
    ######## Synthetic photometry #############
    ###########################################
    
    # turn the total spectrum into a source spectrum object in synphot 
    sp = synphot.SourceSpectrum(Empirical1D, points = angstrom, lookup_table=Flambda * units.FLAM, keep_neg = True)
    
    # loop through all filters
    filt_through = np.zeros(len(filt_list))
    filt_through_Jy = np.zeros(len(filt_list))
    
    # Flam_master = []
    # Jy_master = []
    for f, filt_name in enumerate(filt_list):
        # create an observation and calculate the effective stimulation
        obs = synphot.Observation(sp, filt_dict[filt_name], force = 'taper')
        filt_through[f] = obs.effstim(flux_unit='flam').value
        
        filt_through_Jy[f] = (filt_through[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
        filt_through_Jy[f] = filt_through_Jy[f]/1e6
        
    
    return filt_through_Jy

d21_synphot = calc_syn_phot(draine_data['wave'] * u.micron, total *1e-15, filt_list)
pdrs_synphot = calc_syn_phot(pdrs_data['wavelength'], pdrs_data['flux'], filt_list)


filt_Hz = c/((filt_wave * u.AA).to(u.micron)).value

# for color and symbols
for f, filt_name in enumerate(filt_list):
    if filt_name in pah_filts:
        # gold
        color = '#CC79A7'
        color = '#009E73'
        # color = 'limegreen'
        symb = '^'
        s = 200
    elif filt_name in con_filts:
        # teal
        # color = 'darkorchid'
        color = '#CC79A7'
        symb = 'o'
        s=150
    elif filt_name in pah_con_filts:
        # purple
        # color = '#E69F00'
        color = '#F7B67C'
        # color = 'darkorchid'
        symb = 's'
        s = 150
        
    #     color = '#B156A3'
    # elif filt_name in con_filts:
    #     # teal
    #     color = '#A2D7F1'
    # elif filt_name in pah_con_filts:
    #     # purple
    #     color = '#DECC7C'

    ax1.scatter((filt_wave * u.AA).to(u.micron)[f], d21_synphot[f] * filt_Hz[f], c = color, marker = symb, s = s, edgecolor = 'k', alpha = 0.7, zorder = 100000)
    ax2.scatter((filt_wave * u.AA).to(u.micron)[f], pdrs_synphot[f] * filt_Hz[f], c = color, marker = symb, s = s, edgecolor = 'k', alpha = 0.7, zorder = 100000)
    
# ax3.xaxis.set_major_formatter(ScalarFormatter())
# ax3.xaxis.get_major_formatter().set_scientific(False)
# ax3.xaxis.get_major_formatter().set_powerlimits((0, 1))

ticks = np.linspace(3, xlim[1], 9)
print(ticks)
# ticks = [ 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18]
plt.xticks(ticks, labels=[f'{int(t)}' for t in ticks])

ax3.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax3.xaxis.set_minor_formatter(ticker.NullFormatter())

ax3.set_xlabel('Wavelength ($\mu$m)', size = 22)
ax3.set_ylabel('Transmission', size = 22)

ax3.tick_params(axis='x', which='major', labelsize=14)

# ax1.set_ylabel('$\mathrm{\\nu F_{\\nu}}$ (Hz MJy sr$^{-1}$)', size = 22)
# ax2.set_ylabel('$\mathrm{\\nu F_{\\nu}}$ (Hz MJy sr$^{-1}$)', size = 22)
ax1.set_ylabel('$\mathrm{\\nu F_{\\nu}}$', size = 22)
ax2.set_ylabel('$\mathrm{\\nu F_{\\nu}}$', size = 22)

plt.minorticks_on()

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/plots/'
savename = 'filters_on_data_D21_PDRs4All_v2.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)
