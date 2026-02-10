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
import pandas as pd
import seaborn as sns

plt.ion()
savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/D21_new_phoenix/'
filename = 'D21_combo_SexA_pah_small_box_models_master_table_v2.txt'
plotdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/D21_new_phoenix/'

def mscatter(x,y,ax=None, m=None, **kw):
    import matplotlib.markers as mmarkers
    if not ax: ax=plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                        marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc

modeldir = savedir + 'spec/'
throughdir = savedir + 'throughput/'

if not os.path.exists(throughdir):
    os.makedirs(throughdir)    

if not os.path.exists(modeldir):
    os.makedirs(modeldir)  

data = ascii.read(savedir + filename)

cutoff = 6e-9

plt.figure(1)
plt.clf()
im = plt.scatter(data['U'], data['chisq_Jy'], s = 5, c = data['const'], norm=matplotlib.colors.LogNorm(), cmap = 'rainbow', edgecolor = 'k')
plt.semilogy()
cb = plt.colorbar(im)
cb.set_label('Constant')
plt.axhline(cutoff, c = 'k', lw = 0.5)
plt.xlabel('U')
plt.ylabel('Chisq')

plt.figure(2)
plt.clf()
im = plt.scatter(data['qpah'], data['chisq_Jy'], s = 5, c = data['U'], norm=matplotlib.colors.LogNorm(), cmap = 'rainbow', edgecolor = 'k')
plt.semilogy()
cb = plt.colorbar(im)
cb.set_label('U')
plt.axhline(cutoff, c = 'k', lw = 0.5)
plt.xlabel('qpah')
plt.ylabel('Chisq')

# cutoff = 5e-2


# plt.figure(3)
# plt.clf()
# im = plt.scatter(data['U'], data['chisq_Jy'], s = 5, c = data['const'], norm=matplotlib.colors.LogNorm(), cmap = 'rainbow', edgecolor = 'k')
# plt.semilogy()
# cb = plt.colorbar(im)
# cb.set_label('Constant')
# plt.axhline(cutoff, c = 'k', lw = 0.5)
# plt.xlabel('U')
# plt.ylabel('Chisq')

# plt.figure(4)
# plt.clf()
# im = plt.scatter(data['qpah'], data['chisq_Jy'], s = 5, c = data['U'], norm=matplotlib.colors.LogNorm(), cmap = 'rainbow', edgecolor = 'k')
# plt.semilogy()
# cb = plt.colorbar(im)
# cb.set_label('U')
# plt.axhline(cutoff, c = 'k', lw = 0.5)
# plt.xlabel('qpah')
# plt.ylabel('Chisq')

ind = np.where((data['qpah'] == 1.0) & (data['ion'] == 'lo') & (data['size'] == 'sma'))[0]
U = data['U'][ind]
const = data['const'][ind]
chisq = data['chisq_Jy'][ind]

# U = data['U']
# const = data['const']
# chisq = data['chisq_Jy']

df = pd.DataFrame(dict(x=U, y=const, z=chisq))
xcol, ycol, zcol = "x", "y", "z"
df = df.sort_values(by=[xcol, ycol])
xvals = df[xcol].unique()
yvals = df[ycol].unique()
zvals = df[zcol].values.reshape(len(xvals), len(yvals)).T

plt.figure(5)
plt.clf()
extent = [np.min(xvals), np.max(xvals), np.min(np.log10(yvals)),np.max(np.log10(yvals))]
im = plt.imshow(zvals, cmap = 'rainbow', norm=matplotlib.colors.LogNorm(), extent = extent)
plt.xlabel('U')
plt.ylabel('log C')
cb = plt.colorbar(label = 'log chisq')
plt.title('qpah = 1.0')


cut = np.where(data['chisq_Jy'] < 10)

ion_arr = []
size_arr = []
for i in range(len(data['size'][cut])):
    
    if data['size'][cut][i] == 'lrg':
        size_arr.append(15)
    elif data['size'][cut][i] == 'std':
        size_arr.append(10)
    elif data['size'][cut][i] == 'sma':
        size_arr.append(5)
        
    if data['ion'][cut][i] == 'hi':
        ion_arr.append('P')
    elif data['ion'][cut][i] == 'st':
        ion_arr.append('o')
    elif data['ion'][cut][i] == 'lo':
        ion_arr.append('s')
    

plt.figure(4)
plt.clf()
im = mscatter(data['U'][cut], data['qpah'][cut], s = 30, c = data['chisq_Jy'][cut], cmap = 'rainbow')

# try likelihoods
l = np.exp(-data['chisq_Jy']/2)
plt.figure(7)
plt.clf()
plt.hist(l, bins = 50)

plt.figure(8)
plt.clf()
plt.plot(data['U'], l,'.')


def get_likilihood(var, text, data):
    
    like_array = np.zeros(len(var))
    for i in range(len(var)):
        ind = np.where(data[text] == var[i])
        
        like_array[i] = np.nansum(np.exp(-data['chisq_Jy'][ind]/2))
        
    return like_array
        
    
        
plt.figure(10)
plt.clf()
U_arr = np.unique(data['U'])
like_array = get_likilihood(U_arr, 'U', data)
sumval = np.sum(like_array)
plt.plot(U_arr, like_array/sumval,marker = '.', lw = 4, linestyle = '-', ms = 25)
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('U, radiation field', size = 24)
plt.ylabel('PDF', size = 24)

savepath = '/Users/etarantino/Documents/proposals/JWST_C3/JWST_C3_SexA/plots/'
savename = 'U_PDF_degeneracy.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)

plt.figure(10)
plt.clf()
U_arr = np.unique(data['size'])
like_array = get_likilihood(U_arr, 'size', data)
sumval = np.sum(like_array)
plt.plot(U_arr, like_array/sumval,marker = '.', lw = 4, linestyle = '-', ms = 25)
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('PAH size', size = 24)
plt.ylabel('PDF', size = 24)
plt.ylim(0, 0.4)

savepath = '/Users/etarantino/Documents/proposals/JWST_C3/JWST_C3_SexA/plots/'
savename = 'size_PDF_degeneracy.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)


plt.figure(10)
plt.clf()
U_arr = np.unique(data['ion'])
like_array = get_likilihood(U_arr, 'ion', data)
sumval = np.sum(like_array)
plt.plot(U_arr, like_array/sumval,marker = '.', lw = 4, linestyle = '-', ms = 25)
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('PAH ionization', size = 24)
plt.ylabel('PDF', size = 24)
plt.ylim(0, 0.4)

savepath = '/Users/etarantino/Documents/proposals/JWST_C3/JWST_C3_SexA/plots/'
savename = 'ion_PDF_degeneracy.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)

plt.figure(10)
plt.clf()
U_arr = np.unique(data['qpah'])
like_array = get_likilihood(U_arr, 'qpah', data)
sumval = np.sum(like_array)
plt.plot(U_arr, like_array/sumval,marker = '.', lw = 4, linestyle = '-', ms = 25)
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('q$_{\mathrm{pah}}$', size = 24)
plt.ylabel('PDF', size = 24)
plt.ylim(0, 0.4)

savepath = '/Users/etarantino/Documents/proposals/JWST_C3/JWST_C3_SexA/plots/'
savename = 'qpah_PDF_degeneracy.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)

plt.figure(10)
plt.clf()
U_arr = np.unique(data['const'])
like_array = get_likilihood(U_arr, 'const', data)
sumval = np.sum(like_array)
plt.plot(U_arr/1e22, like_array/sumval,marker = '.', lw = 3, linestyle = '-', ms = 20)
ax = plt.gca()
ax.tick_params(axis='both', which='major', labelsize=18)
plt.xlabel('N$_{\mathrm{HI}}$ (10$^{23}$ cm$^{-2}$)', size = 24)
plt.ylabel('PDF', size = 24)
plt.ylim(0, 0.04)

savepath = '/Users/etarantino/Documents/proposals/JWST_C3/JWST_C3_SexA/plots/'
savename = 'NHI_PDF_degeneracy.pdf'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)

# plt.semilogy()
# cb = plt.colorbar(im)
# cb.set_label('Constant')


ind = np.argmin(data['chisq'])
print(data[ind])


regname = 'SexA_pah_small_box'
regdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/obs/'
regfile = regdir + f'sexa_SED_region_{regname}.txt'
jwst_data = ascii.read(regfile)


filt_list = ['F115W', 'F150W', 'F200W','F300M', 'F335M', 'F360M', 'F560W', 'F770W', 'F1000W', 'F1130W', 'F1500W']
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


# cutoff = 0.5
# ind = np.where(data['chisq_Jy'] < cutoff)[0]

# # U_search = 2.5
# # ind = np.where((data['U'] == U_search)  & (data['qpah'] == 0.5) & (data['ion'] == 'lo') & (data['size'] == 'lrg'))[0]

# # Usavepath = plotdir + '{:01.2f}/'.format(U_search)
# Usavepath = plotdir + 'low_chisq_Jy2/'
# if not os.path.exists(Usavepath):
#     os.makedirs(Usavepath)    

# for i in range(len(ind)):

#     print(data[ind[i]])
    
#     U = data['U'][ind[i]]
#     qpah = data['qpah'][ind[i]]
#     const = np.log10(data['const'][ind[i]])
#     ion = data['ion'][ind[i]]
#     size = data['size'][ind[i]]
#     chisq =  data['chisq_Jy'][ind[i]]
    
    
#     model_name = 'D21_combo_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.txt'.format(U, qpah, const, ion, size)
#     model = ascii.read(modeldir + model_name)
    
#     through_name = 'D21_combo_throughput_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.txt'.format(U, qpah, const, ion, size)
#     through = ascii.read(throughdir + through_name)
#     through_Jy = (through['throughput'] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(micron_mid[3:] * u.micron))
#     through_MJy = through_Jy/1e6
    
#     model_Jy = (model['total'] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(model['microns'] * u.micron))
#     model_MJy = model_Jy/1e6
    
#     D21_Jy = (model['D21_Flam'] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(model['microns'] * u.micron))
#     D21_MJy = D21_Jy/1e6
    
#     ph_Jy = (model['phoenix'] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(model['microns'] * u.micron))
#     ph_MJy = ph_Jy/1e6
    
#     plt.figure(6)
#     plt.clf()
#     plt.plot(model['microns'],  model['total'], alpha = 0.5)
#     # plt.plot(model['microns'], model_MJy, c = 'b', alpha = 0.7)
#     # plt.plot(model['microns'], D21_MJy, c = 'green', alpha = 0.7)
#     # plt.plot(model['microns'], ph_MJy, c = 'purple', alpha = 0.7)
    
#     plt.title('U = {:01.2f} qpah = {:01.2f} logC = {:3.1f} ion = {:s} size = {:s} chisq = {:2.1e}'.format(U, qpah, const, ion, size, chisq))
#     plt.xlim(0,20)
#     # plt.ylim(0, np.max(MJy_flux) + 1.5*np.max(MJy_flux))
    
#     plt.scatter(micron_mid, Flam_flux, s = 10, c = 'k', label = 'Observed')
#     plt.scatter(micron_mid[3:], through['throughput'], s = 10, c = 'r')
#     # plt.scatter(micron_mid[3:], through_MJy, s = 10, c = 'r', label = 'Model Throughput')
#     plt.legend(loc = 'best')
#     plt.xlabel('Microns')
#     plt.ylabel('MJy/sr')
    
#     savename = 'D21_visualize_MJy_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.pdf'.format(U, qpah, const, ion, size)
#     # plt.savefig(Usavepath + savename)
    

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

    


    
    