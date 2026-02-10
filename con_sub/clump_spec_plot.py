#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates proposal plot

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import matplotlib
import pandas as pd
import seaborn as sns
import matplotlib.colors
from astropy.constants import c

import synphot
from synphot.models import Empirical1D

grid_name = 'clump_phoenix_D21_metallicity'

# load clump file that has the fluxes for each clump
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_file = 'clump_props_pah_filt_flux_corrected.txt'
clumps = ascii.read(clump_path + clump_file)

clump_num = clumps['clump_num']

clump_to_model = 1

clump_num = clumps['clump_num']


clumps = clumps[clump_num == clump_to_model]

plt.ion()
griddir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/' + grid_name + '/clump_{:d}/'.format(clump_to_model)
filename = 'D21_combo_clump{:d}_models_master_table.txt'.format(clump_to_model)
data = ascii.read(griddir + filename)

regname = 'SexA_pah_small_box'
regdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/obs/'
regfile = regdir + f'sexa_SED_region_{regname}.txt'
jwst_data = ascii.read(regfile)
jwst_data.sort('mu_low')

modeldir = griddir + 'spec/'
throughdir = griddir + 'throughput/'


miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

filt_list = ['F115W', 'F150W', 'F200W','F300M', 'F335M', 'F360M', 'F560W', 'F770W', 'F1000W', 'F1130W', 'F1500W']
Jy_flux = np.zeros(len(filt_list))
Jy_flux_err = np.zeros(len(filt_list))
Flam_flux = np.zeros(len(filt_list))
micron_mid = np.zeros(len(filt_list))

ind = np.where(data['chisq_Jy'] <10000)[0]

plt.figure(1)
plt.clf()
plt.hist(data['chisq_Jy'][ind], bins = 50)
                

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

    micron_mid[i] = filt_feature.pivot().to(u.micron).value
    
    Jy_flux[i] = clumps[filt_name + '_flux']
    
    Jy_flux_err[i] = clumps[filt_name + '_err']
    
    # convert the flux to Flam units 
    Flam_flux[i] = (Jy_flux[i] * 1e6 * u.Jy).to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(micron_mid[i] * u.micron)).value

# checking k method math
print('PAH', clumps['clump_flux_3_k1'])
print('CON', clumps['clump_con_3_k1'])
print('PAH + CON', clumps['clump_flux_3_k1'] + clumps['clump_con_3_k1'])
print('F770W', clumps['F335M_flux'])


cutoff = 8000
ind = np.where(data['chisq_Jy'] <cutoff)[0]


# plot all the parameters relative to chisquared
plt.figure(3)
plt.clf()
im = plt.scatter(data['U'], data['chisq_Jy'], s = 2, c = data['const'], norm=matplotlib.colors.LogNorm(), cmap = 'rainbow')
plt.xlabel('U')
plt.ylabel('Chisq')
plt.semilogy()
cb = plt.colorbar(im)
cb.set_label('Const')
plt.axhline(cutoff, c = 'k', lw = 0.5)

plt.figure(4)
plt.clf()
im = plt.scatter(data['size'], data['chisq_Jy'], s = 2, c = data['qpah'], cmap = 'rainbow')
plt.xlabel('size')
plt.ylabel('Chisq')
plt.semilogy()
cb = plt.colorbar(im)
cb.set_label('qpah')
plt.axhline(cutoff, c = 'k', lw = 0.5)


# load an initial model
i = 1
U = data['U'][ind[i]]
qpah = data['qpah'][ind[i]]
const = np.log10(data['const'][ind[i]])
ion = data['ion'][ind[i]]
size = data['size'][ind[i]]
chisq =  data['chisq_Jy'][ind[i]]

model_name = 'D21_combo_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.txt'.format(U, qpah, const, ion, size)
model = ascii.read(modeldir + model_name)

master_mod = np.zeros(len(model))

all_mods_flux = []
all_mods_wave = []
all_mods_7_flux = []
all_through_Jy = []
all_chi = []
all_ion = []
all_size = []
all_qpah = []
all_U = []
all_const = []

plt.figure(7, figsize = (10,5))
plt.clf()

for i in range(len(ind)):

    print(data[ind[i]])
    
    U = data['U'][ind[i]]
    qpah = data['qpah'][ind[i]]
    const = np.log10(data['const'][ind[i]])
    ion = data['ion'][ind[i]]
    size = data['size'][ind[i]]
    chisq =  data['chisq_Jy'][ind[i]]
    
    model_name = 'D21_combo_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.txt'.format(U, qpah, const, ion, size)
    model = ascii.read(modeldir + model_name)
    
    through_name = 'D21_combo_throughput_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.txt'.format(U, qpah, const, ion, size)
    through = ascii.read(throughdir + through_name)
    through_Jy = (through['throughput'] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(micron_mid[3:] * u.micron))
    through_MJy = through_Jy/1e6
    
    model_Jy = (model['total'] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(model['microns'] * u.micron))
    
    master_mod = master_mod + model_Jy  

    
    large_7_val = np.where((model['microns'] > 7) & (model['microns'] < 8))[0]
    norm_val = np.nanmax(model_Jy[large_7_val])
    
    all_mods_7_flux.append(norm_val)
    
    
    all_through_Jy.append(through_Jy)
    all_mods_wave.append(model['microns'])
    all_mods_flux.append(model_Jy)
    all_chi.append(chisq)
    all_ion.append(ion)
    all_size.append(size)
    all_qpah.append(qpah)
    all_U.append(U)
    all_const.append(const)
    

# all_mods_7_flux.sort()
index = np.arange(0, len(ind), 1)
ind_sort = [i[0] for i in sorted(zip(index, all_mods_7_flux), key=lambda l: l[1], reverse=False)]


# sort_all_mods_wave = [i[0] for i in sorted(zip(all_mods_wave, all_mods_7_flux), key=lambda l: l[1], reverse=False)]
# sort_all_mods_flux = [i[0] for i in sorted(zip(all_mods_flux, all_mods_7_flux), key=lambda l: l[1], reverse=False)]

# U_arr = np.arange(0, 4, 0.5)
# colors = sns.color_palette("hls", len(U_arr))
# colors = sns.color_palette("hls", len(ind))

qpah_arr = np.arange(0, 1.5, 0.25)
colors = sns.color_palette("hls", len(qpah_arr))

    
sma_count = 0
std_count = 0
lrg_count = 0

qpah_05_count = 0
qpah_1_count = 0
qpah_15_count = 0

sort_mu_low = np.copy(jwst_data['mu_low'])
sort_mu_high = np.copy(jwst_data['mu_high'])
sort_mu_mid = np.copy(jwst_data['mu_mid'])
sort_mu_low.sort()
sort_mu_high.sort()
sort_mu_mid.sort()

for i in range(len(ind)):
    
    this_mods_flux = all_mods_flux[ind_sort[i]]
    this_mods_wave = all_mods_wave[ind_sort[i]]
    this_through_Jy = all_through_Jy[ind_sort[i]]
    this_chi = all_chi[ind_sort[i]]
    this_ion = all_ion[ind_sort[i]]
    this_size = all_size[ind_sort[i]]
    this_qpah = all_qpah[ind_sort[i]]
    this_U = all_U[ind_sort[i]]
    
    # if this_size == 'sma':
    #     if this_ion == 'lo':
    #         color = colors[0]
    #     elif this_ion == 'st':
    #         color = colors[1]
    #     elif this_ion == 'hi':
    #         color = colors[2]
    # elif this_size == 'std':
    #     if this_ion == 'lo':
    #         color = colors[3]
    #     elif this_ion == 'st':
    #         color = colors[4]
    #     elif this_ion == 'hi':
    #         color = colors[5]
    # elif this_size == 'lrg':
    #     if this_ion == 'lo':
    #         color = colors[6]
    #     elif this_ion == 'st':
    #         color = colors[7]
    #     elif this_ion == 'hi':
    #         color = colors[8]

    if this_size == 'sma':
        color = 'purple'
        if this_ion == 'lo':
            color = 'blue'
        elif this_ion == 'st':
            color = 'blueviolet'
        elif this_ion == 'hi':
            color = 'cornflowerblue'
            
        sma_count +=1
    elif this_size == 'std':
        color = 'green'
        if this_ion == 'lo':
            color = 'green'
        elif this_ion == 'st':
            color = 'lime'
        elif this_ion == 'hi':
            color = 'aquamarine'
        std_count +=1
        
    elif this_size == 'lrg':
        color = 'orange'
        if this_ion == 'lo':
            color = 'orange'
        elif this_ion == 'st':
            color = 'goldenrod'
        elif this_ion == 'hi':
            color = 'red'
        lrg_count +=1
        
    # color_ind = int(np.floor(this_chi*10))
    # color = colors[color_ind]    
    
    # this_ind = np.where(this_U == U_arr)[0][0]
    # color = colors[this_ind]
    
    if this_qpah == 0.5:
        color = 'blue'
        qpah_05_count +=1
    elif this_qpah == 1.0:
        color = 'orange'
        qpah_1_count +=1
    elif this_qpah == 1.5:
        color = 'black'
        qpah_15_count +=1
        
    # this_ind = np.where(this_qpah == qpah_arr)[0][0]
    # color = colors[this_ind]
    
    
    # if this_U == 0:
    #     color = colors[0]
    # elif this_U == 0.5:
    #     color = colors[1]
    # elif this_U == 1:
    #     color = colors[2]
    # elif this_U == 1.5:
    #     color = colors[3]
    # elif this_U == 2:
    #     color = colors[4]
    # elif this_U == 2.5:
    #     color = colors[5]
    # elif this_U == 3:
    #     color = colors[6]
    # elif this_U == 3.5:
    #     color = colors[7]
    
    
    max_chi = np.max(data['chisq_Jy'][ind])
    min_chi = np.min(data['chisq_Jy'][ind])
    norm = 1/((min_chi/max_chi))
    alpha_val = ((1/((this_chi/max_chi)))/norm)/2
    
    
    # alpha_val = 0.5
    
    # alpha_val = 0.3
    
        
    plt.figure(7)
    # plt.plot(this_mods_wave, this_mods_flux*1e6, c = 'gray', alpha = alpha_val, zorder = 0.0001, lw = 0.3)
    # plt.plot(model['microns'], D21_MJy, c = 'green', alpha = 0.7)
    # plt.plot(model['microns'], ph_MJy, c = 'purple', alpha = 0.7)
    
    # plot throughput
    # im = plt.scatter(micron_mid[3:], this_through_Mjy, marker = '.', alpha = 0.6, c = colors[i])
    # plt.hlines(this_through_Jy*1e6, sort_mu_low[3:], sort_mu_high[3:], lw = 1, zorder = 10, alpha = alpha_val, color = 'cornflowerblue')
    # plt.hlines(this_through_Jy, sort_mu_low[3:], sort_mu_high[3:], lw = 1, zorder = 10, alpha = alpha_val, color = color)


    # convert to nu * Fnu
    nu_this_mods_wave = (c/(this_mods_wave * u.micron)).to(u.Hz)
    plt.plot(this_mods_wave, this_mods_flux*1e-23*nu_this_mods_wave*1e15, c = 'gray', alpha = alpha_val, zorder = 0.0001, lw = 0.3)

    
    # plot throughput
    nu_mu_mid = (c/(sort_mu_mid[3:] * u.micron)).to(u.Hz) 
    plt.hlines(this_through_Jy*1e-23*nu_mu_mid*1e15, sort_mu_low[3:], sort_mu_high[3:], lw = 1, zorder = 10, alpha = alpha_val, color = '#B4A0FF')
    # plt.hlines(this_through_Jy, sort_mu_low[3:], sort_mu_high[3:], lw = 1, zorder = 10, alpha = alpha_val, color = color)


# nonsense to plot trends in the SED plot
# plt.colorbar(im)
    
    # plt.title('U = {:01.2f} qpah = {:01.2f} logC = {:3.1f} ion = {:s} size = {:s} chisq = {:2.1e}'.format(U, qpah, const, ion, size, chisq))
    # plt.ylim(0, np.max(MJy_flux) + 1.5*np.max(MJy_flux))
        # plt.scatter(micron_mid[3:], through['throughput'], s = 10, c = 'r')
    # plt.scatter(micron_mid[3:], through_MJy, s = 10, c = 'r', label = 'Model Throughput')

    
    # save name = 'D21_visualize_MJy_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.pdf'.format(U, qpah, const, ion, size)
    # plt.savefig(Usavepath + savename)
    
print('SMALL', sma_count)
print('STD', std_count)
print('LARGE', lrg_count)

print('qpah=0.5', qpah_05_count)
print('qpah=1.0', qpah_1_count)
print('qpah=1.5', qpah_15_count)


master_mod = master_mod/len(ind)

c_list = []
for i in range(len(jwst_data)):
    filt = jwst_data['filt'][i]

    if filt == 'F335M':
        # c_list.append('#ff5fa2')
        c_list.append('#009E73')
        lw = 3
    elif filt == 'F770W':
        # c_list.append('#ff5fa2')
        c_list.append('#009E73')
        lw = 3
    elif filt == 'F1130W':
        # c_list.append('#ff5fa2')
        c_list.append('#009E73')
        lw = 3
        
    else:
        # c_list.append('springgreen')
        c_list.append('#CC79A7') # magenta
        # c_list.append('#E69F00') # goldenrod
        # c_list.append('#0072B2') # deep blue
        # c_list.append('#CC79A7') # rose
        # c_list.append('#F7B67C') # gold
        
    
        # E69F00, CC79A7
        lw = 3
        
    color = matplotlib.colors.to_rgb(c_list[i])
        
    
    nu_mu_mid = (c/(jwst_data['mu_mid'][i] * u.micron)).to(u.Hz).value
    # plt.hlines(Jy_flux[i]*1e6, jwst_data['mu_low'][i], jwst_data['mu_high'][i], lw = lw, zorder = 10000, color = color, label = 'Observed Data')
    plt.hlines(Jy_flux[i]*1e-23*nu_mu_mid*1e15, jwst_data['mu_low'][i], jwst_data['mu_high'][i], lw = lw, zorder = 10000, color = color, label = 'Observed Data')
    plt.errorbar(jwst_data['mu_mid'][i], Jy_flux[i]*1e-23*nu_mu_mid*1e15, yerr = Jy_flux_err[i]*1e-23*nu_mu_mid*1e15, fmt = '.', ecolor = color, ms = 0.01, elinewidth = 3, capsize = 3, zorder = 1000000000)

# plt.scatter(micron_mid, Jy_flux)
    
nu_model_micron = (c/(model['microns'] * u.micron)).to(u.Hz).value
# plt.plot(model['microns'], master_mod*1e6, c = 'k', alpha = 1, lw = 2, zorder  = 10, label = 'Average Model')
plt.plot(model['microns'], master_mod*1e-23 *nu_model_micron*1e15, c = 'k', alpha = 1, lw = 2, zorder  = 10, label = 'Average Model')



plt.minorticks_on()

# plt.legend(loc = 'best')
plt.xlabel('Wavelength ($\mu$m)', size = 16)
# plt.ylabel('$\mathrm{F_{\\nu}}$ ($\mathrm{\mu}$Jy)', size = 18)
plt.ylabel('$\mathrm{\\nu F_{\\nu}}$ ($\mathrm{10^{-15} \ erg \, s^{-1} \, cm^{-2}}$)', size = 16)

plt.xlim(0.5,20)
# plt.ylim(0, 55)
plt.ylim(1.2, 15)

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/SEDs/' + grid_name + '/'

# # report the range of parameters preferred
# plt.figure(8)
# plt.clf()
# all_qpah = np.round(all_qpah, decimals = 1)
# plt.hist(all_qpah, align = 'mid', bins = 7)
# plt.xlabel('qpah')
# savename = 'D21_phoenix_clump_{:d}_qpah.png'.format(clump_to_model)
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 300)

# plt.figure(9)
# plt.clf()
# all_U = np.round(all_U, decimals = 1)
# plt.hist(all_U, align = 'mid', bins = 8)
# plt.xlabel('U')
# savename = 'D21_phoenix_clump_{:d}_U.png'.format(clump_to_model)
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 300)

# plt.figure(10)
# plt.clf()
# all_size = np.array(all_size)
# print('Small models', len(all_size[all_size=='lrg'])/len(all_size))
# plt.hist(all_size, align = 'mid', bins = 3)
# savename = 'D21_phoenix_clump_{:d}_size.png'.format(clump_to_model)
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 300)

# plt.figure(11)
# plt.clf()
# all_ion = np.array(all_ion)
# print('High models', len(all_ion[all_ion=='lo'])/len(all_ion))
# plt.hist(all_ion, align = 'mid', bins = 3)
# savename = 'D21_phoenix_clump_{:d}_ion.png'.format(clump_to_model)
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 300)


# plt.figure(12)
# plt.clf()
# plt.hist(all_const, align = 'mid', bins = 10)


savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/SEDs/' + grid_name + '/'
if not os.path.exists(savepath):
    os.makedirs(savepath)    
savename = 'D21_phoenix_clump_{:d}_SED_nu_Fnu_metallicity_corr_purple.pdf'.format(clump_to_model)
# savename = 'D21_phoenix_clump_{:d}_SED_Fnu.png'.format(clump_to_model)
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 600)

ETC = np.array([model['microns'], master_mod * 1e9]).T
ETC_name = 'ETC_SexA_model_average_spec_clump_{:d}.txt'.format(clump_to_model)
names = ['microns', 'mJy']
# ascii.write(ETC, savepath + ETC_name, names = names, overwrite = True)