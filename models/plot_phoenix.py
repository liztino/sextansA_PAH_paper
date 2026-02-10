#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plots the results from comparing the grid of phoenix models 

@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import matplotlib

plotdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/phoenix_more_filts/'


regname = 'SexA_pah_small_box'
regdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/obs/'
regfile = regdir + f'sexa_SED_region_{regname}.txt'
jwst_data = ascii.read(regfile)


savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/phoenix_more_filts/'
filename = f'phoenix_models_{regname}_master_table.txt'

modeldir = savedir + 'spec/'
throughdir = savedir + 'throughput/'

data = ascii.read(savedir + filename)

cutoff = 5e-9


plt.figure(1)
plt.clf()
im = plt.scatter(data['T'], data['Chisq'], s = 2, c = data['Const'], norm=matplotlib.colors.LogNorm(), cmap = 'rainbow')
plt.semilogy()
plt.colorbar(im)
plt.axhline(cutoff, c = 'k', lw = 0.5)

plt.figure(2)
plt.clf()
im = plt.scatter(data['logg'], data['Chisq'], s = 2, c = data['T'], norm=matplotlib.colors.LogNorm(), cmap = 'rainbow')
plt.semilogy()
plt.colorbar(im)
plt.axhline(cutoff, c = 'k', lw = 0.5)

filt_list = ['F115W', 'F150W', 'F200W', 'F300M', 'F360M']
MJy_flux = np.zeros(len(filt_list))
Flam_flux = np.zeros(len(filt_list))
micron_mid = np.zeros(len(filt_list))

                
for i, filt_name in enumerate(filt_list):
    
    # find the JWST for the filter
    arg = np.where(jwst_data['filt'] == filt_name)[0]
    MJy_flux[i] = jwst_data['med'][arg].value
    micron_mid[i] = jwst_data['mu_mid'][arg].value
    
    # convert the flux to Flam units 
    Flam_flux[i] = (MJy_flux[i] * 1e6 * u.Jy).to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(micron_mid[i] * u.micron)).value


ind = np.argmin(data['Chisq'])
print(data[ind])

ind = np.where(data['Chisq'] < cutoff)[0]


T_search = 5500
ind = np.where((data['T'] == T_search))[0]

Tsavepath = plotdir + '{:04.0f}/'.format(T_search)
if not os.path.exists(Tsavepath):
    os.makedirs(Tsavepath)    


Jy_chisq = np.ones(len(data['Chisq']))

for i in range(len(ind)):

    print(data[ind[i]])
    
    T = data['T'][ind[i]]
    logg = data['logg'][ind[i]]
    const = np.log10(data['Const'][ind[i]])
    chisq =  data['Chisq'][ind[i]]
    
    ph_name = 'phoenix_{:4.0f}_{:2.1f}_{:3.2f}.txt'.format(T, logg, const)
    ph_model = ascii.read(modeldir + ph_name)
    
    through_name = 'phoenix_throughput_{:4.0f}_{:2.1f}_{:3.2f}.txt'.format(T, logg, const)
    through = ascii.read(throughdir + through_name)
    
    ph_Jy = (ph_model['Flambda'] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(ph_model['microns'] * u.micron))
    ph_MJy = ph_Jy/1e6
    
    through_Jy = (through['throughput'] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(micron_mid * u.micron))
    through_MJy = through_Jy.value/1e6
    
    plt.figure(3)
    plt.clf()
    # plt.plot(ph_model['microns'], ph_model['Flambda'], alpha = 0.5)
    plt.plot(ph_model['microns'], ph_MJy, alpha = 0.5)
    plt.title('T = {:4.0f} logg = {:2.1f} logC = {:3.2f}'.format(T, logg, const))
    plt.xlim(0,10)
    
    plt.scatter(micron_mid, MJy_flux, s = 10, c = 'k', label = 'data')
    plt.scatter(micron_mid, through_MJy, s = 10, c = 'r', label = 'model through')
    
    # plt.scatter(micron_mid, Flam_flux, s = 10, c = 'k', label = 'data')
    # plt.scatter(micron_mid, through['throughput'], s = 10, c = 'r', label = 'model through')
    
    Jy_chisq[ind[i]] = np.nansum((MJy_flux[0:3] - through_MJy[0:3])**2/(through_MJy[0:3]))

    plt.title('T = {:4.0f} logg = {:2.1f} const = {:3.2f} chisq={:3.2e}.txt'.format(T, logg, const, Jy_chisq[ind[i]]))
              
    plt.legend(loc = 'best')
    plt.xlabel('Microns')
    plt.ylabel('MJy/sr')
    
    savename = 'phoenix_visualize_MJy_{:s}_{:4.0f}_{:2.1f}_{:3.2f}.pdf'.format(regname, T, logg, const)
    plt.savefig(Tsavepath + savename)

# temp_list = np.arange(2000, 6000, 500)

# for i in range(len(temp_list)):
#     ind1 = np.where(data['T'] == temp_list[i])
    
#     ind2 = np.argmin(data['Chisq'][ind1])
    
#     print(data['T'][ind1][ind2], data['Chisq'][ind1][ind2])
    
#     T = data['T'][ind1][ind2]
#     logg = data['logg'][ind1][ind2]
#     const = np.log10(data['Const'][ind1][ind2])
    
#     ph_name = 'phoenix_{:4.0f}_{:2.1f}_{:3.1f}.txt'.format(T, logg, const)
#     ph_model = ascii.read(modeldir + ph_name)
    
#     through_name = 'phoenix_throughput_{:4.0f}_{:2.1f}_{:3.1f}.txt'.format(T, logg, const)
#     through = ascii.read(throughdir + through_name)
    
#     plt.figure()
#     plt.clf()
#     plt.plot(ph_model['microns'], ph_model['Flambda'], alpha = 0.5)
#     plt.title('T = {:4.0f} logg = {:2.1f} logC = {:3.1f}'.format(T, logg, const))
#     plt.xlim(0,10)
    
#     plt.scatter(micron_mid, Flam_flux, s = 10, c = 'k')
#     plt.scatter(micron_mid, through['throughput'], s = 10, c = 'r')

    


    
    