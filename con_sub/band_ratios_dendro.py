#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates band ratio plots and images from the data


@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
# import reproject
from astropy.wcs import WCS
import numpy.ma as ma
import astropy.visualization.wcsaxes as vis
import astropy.units as u 
import matplotlib
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
import seaborn as sns
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import StrMethodFormatter, NullFormatter, LogLocator, ScalarFormatter
from matplotlib.colors import LogNorm

plt.ion()

# custom functions 
from get_pivot_wave import get_pivot_wave
import k_eq

clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'

clump_file = 'clump_props_pah_filt_flux_corrected.txt'

clumps = ascii.read(clump_path + clump_file)

# remove clump 14 just for a quick test with the Halpha data
# clumps.remove_row(13)


D21_sp = 'BC03_Z0.0004_10Myr_'


def plot_D21_mods(U_cut=3.5, color = 'k'):
    
    # load the band ratios from the D21 k analysis for PAHs
    
    inpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/D21_models/'
    D21_F335M_table = D21_sp + 'F335M_D21_models_with_k_consub_new_D21_k_D21.txt'
    D21_F335M = ascii.read(inpath + D21_F335M_table )  
    print(D21_F335M.colnames)

    D21_F770W_table = D21_sp + 'F770W_D21_models_with_k_consub_new_D21_k_D21.txt'
    D21_F770W = ascii.read(inpath + D21_F770W_table )  

    D21_F1130W_table = D21_sp + 'F1130W_D21_models_with_k_consub_new_D21_k_D21.txt'
    D21_F1130W = ascii.read(inpath + D21_F1130W_table )  

    # cut each table to have U<3
    D21_F335M = D21_F335M[D21_F335M['U'] < U_cut]
    D21_F770W = D21_F770W[D21_F770W['U'] < U_cut]
    D21_F1130W = D21_F1130W[D21_F1130W['U'] < U_cut]
    
    D21_arr = [D21_F335M, D21_F770W, D21_F1130W]
    
    size_arr = ['sma','std', 'lrg']
    ion_arr = ['lo', 'st', 'hi']
    
    prop_F335M = np.zeros((3,3))
    prop_F770W = np.zeros((3,3))
    prop_F1130W = np.zeros((3,3))
    
    prop_arr = [prop_F335M, prop_F770W, prop_F1130W]
    
    prop_label = np.zeros((3,3), dtype = '<U11')
    
    
    i = 0
    for size in size_arr:
        j = 0
        for ion in ion_arr:
    
            ind_F335M = np.where((D21_F335M['size'] == size_arr[i]) & (D21_F335M['ion'] == ion_arr[j]) )
            prop_F335M[i][j] = np.nanmean(D21_F335M[ind_F335M]['pah_k_method'])
            
            ind_F770W = np.where((D21_F770W['size'] == size_arr[i]) & (D21_F770W['ion'] == ion_arr[j]) )
            prop_F770W[i][j] = np.nanmean(D21_F770W[ind_F770W]['pah_k_method'])
                            
            ind_F1130W = np.where((D21_F1130W['size'] == size_arr[i]) & (D21_F1130W['ion'] == ion_arr[j]) )
            prop_F1130W[i][j] = np.nanmean(D21_F1130W[ind_F1130W]['pah_k_method'])
            
            prop_label[i][j] = size_arr[i] + '-' + ion_arr[j]
            
            
            j = j+1
            
        i = i + 1

    
    # plot the models on top
    plt.scatter(prop_F335M/prop_F770W, prop_F335M/prop_F1130W, color = color, s = 30, alpha = 0.5)
    plt.plot(prop_F335M/prop_F770W, prop_F335M/prop_F1130W, color = color, alpha = 0.5)
    for i in range(3):
        plt.plot(prop_F335M[i][:]/prop_F770W[i][:], prop_F335M[i][:]/prop_F1130W[i][:], color = color, alpha = 0.5)
        
    # # plot the models on top
    # plt.scatter(prop_F335M/prop_F770W, prop_F335M/prop_F1130W, color = color, s = 30, alpha = 0.5)
    # plt.plot(prop_F335M/prop_F770W, prop_F335M/prop_F1130W, color = color, alpha = 0.5)
    # for i in range(3):
    #     plt.plot(prop_F335M[i][:]/prop_F770W[i][:], prop_F335M[i][:]/prop_F1130W[i][:], color = color, alpha = 0.5)

    label = False
    
    if label:
        i = 0
        for size in size_arr:
            j = 0
            for ion in ion_arr:
                plt.text(prop_F335M[i][j]/prop_F770W[i][j], prop_F335M[i][j]/prop_F1130W[i][j], s = prop_label[i][j])
                
                j = j+1
                
            i = i + 1



pah3_sigma = clumps['clump_flux_3_k1_err_k']
pah7_sigma = clumps['clump_flux_7_k1_err_k']
pah11_sigma = clumps['clump_flux_11_k1_err_k']

pah3 = clumps['clump_flux_3_k1']
pah7 = clumps['clump_flux_7_k1']
pah11 = clumps['clump_flux_11_k1']

err_x = np.sqrt((pah3_sigma/pah3)**2 + (pah7_sigma/pah7)**2) * (pah3/pah7)
err_y = np.sqrt((pah3_sigma/pah3)**2 + (pah11_sigma/pah11)**2) * (pah3/pah11)

# load the uv data to color points
anc_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
uv_name = 'sextans-a_m2_cr.fits'
uv_name = 'sextans-a_w1_cr.fits'
# uv_name = 'sextans-a_w2_cr.fits'



uv_hdr = fits.open(anc_path + uv_name)
uv_data = uv_hdr[0].data
uv_head = uv_hdr[0].header
uv_w = WCS(uv_head)

# find world coordinates for all of the clumps
coord = SkyCoord(ra = clumps['ra_cen'] * u.deg, dec = clumps['dec_cen']* u.deg)

# find pixel coordinates for the uv data
uv_x, uv_y = uv_w.world_to_pixel(coord)
uv_x = np.int32(uv_x)
uv_y = np.int32(uv_y)

# load the Halpha data to color points
anc_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
ha_name = 'D09_SextansA_Halpha_sub_units.fits'

ha_hdr = fits.open(anc_path + ha_name)
ha_data = ha_hdr[0].data
ha_head = ha_hdr[0].header
ha_w = WCS(ha_head)

# convert to something useable
ha_data = ha_data*1e15

# find world coordinates for all of the clumps
coord = SkyCoord(ra = clumps['ra_cen'] * u.deg, dec = clumps['dec_cen']* u.deg)

# find pixel coordinates for the uv data
ha_x, ha_y = ha_w.world_to_pixel(coord)
ha_x = np.int32(ha_x)
ha_y = np.int32(ha_y)

# load the HI
anc_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
hi_name = 'DDO75_R_X0_P_R_FLATTEN.FITS'

hi_hdr = fits.open(anc_path + hi_name)
hi_data = hi_hdr[0].data
hi_head = hi_hdr[0].header
hi_w = WCS(hi_head)


# find world coordinates for all of the clumps
coord = SkyCoord(ra = clumps['ra_cen'] * u.deg, dec = clumps['dec_cen']* u.deg)

# find pixel coordinates for the uv data
hi_x, hi_y = hi_w.world_to_pixel(coord)
hi_x = np.int32(hi_x)
hi_y = np.int32(hi_y)


# plt.figure(5, figsize = (8, 6))
# plt.clf()
# plt.errorbar(F335M_k_2['flux']/F770W_k_4['flux'], F335M_k_2['flux']/F1130W_k_10['flux'],
#              yerr = err_y, xerr = err_x, fmt = '.', ecolor = 'k', capsize = 2, zorder = 1, markersize = 0.01)
# im = plt.scatter(F335M_k_2['flux']/F770W_k_4['flux'], F335M_k_2['flux']/F1130W_k_10['flux'], 
#             c = np.log10(uv_data[uv_y, uv_x]), s = 100)
# for i in range(len(master_clump)):
#     plt.text(F335M_k_2['flux'][i]/F770W_k_4['flux'][i], F335M_k_2['flux'][i]/F1130W_k_7['flux'][i],
#              s = '{:d}'.format(master_clump['Clump'][i]))

plt.figure(1, figsize = (8, 6))
plt.clf()
plt.errorbar(pah3/pah7, pah3/pah11,
             yerr = err_y, xerr = err_x, fmt = '.', ecolor = 'k', capsize = 2, zorder = 1, markersize = 0.01)

# color by UV
# im = plt.scatter(pah3/pah7, pah3/pah11, c = np.log10(uv_data[uv_y, uv_x]), s = 100)

# color by Halpha
im = plt.scatter(pah3/pah7, pah3/pah11, c = ha_data[ha_y, ha_x], s = 100, cmap = 'viridis', label = 'Sextans A (this paper)')


# # color by Sigma PAH/F1500W in the clump 
# sigma_pah_F1500W = (pah3 + pah7 + pah11)/clumps['F1500W_flux']
# im = plt.scatter(pah3/pah7, pah3/pah11, c = sigma_pah_F1500W, s = 80, cmap = 'plasma', edgecolors = 'k', linewidth = 0.3)

# colorbar label
cb = plt.colorbar(im)
# cb.set_label('log(UV)')
cb.set_label('H$\\alpha$ Flux (10$^{-15}$ erg s$^{-1}$ cm$^{-2}$)', size = 'large')
# cb.set_label('$\Sigma$PAH/F1500W')
    
plt.minorticks_on()
plt.xlabel('3.3/7.7 PAH', size = 'x-large')
plt.ylabel('3.3/11.3 PAH', size = 'x-large')

# add clump labels
ax = plt.gca()
# for i in range(len(clumps)):
#     ax.text(pah3[i]/pah7[i], pah3[i]/pah11[i], s = clumps['clump_num'][i])

# # add multiple U models
# U_list = np.arange(0, 1, 1)
# c_list = sns.color_palette("viridis", len(U_list))

# for i in range(len(U_list)):
#     plot_D21_mods(U_cut=U_list[i], color = 'k') 

plot_D21_mods() 
    
    
# # plot PDRs4all
inpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/D21_models/'
pdrs_F335M  = ascii.read(inpath + 'PDRs4All_F335M_D21_models_with_k_consub_new_pdrs4all_k.txt' )  
pdrs_F770W  = ascii.read(inpath + 'PDRs4All_F770W_D21_models_with_k_consub_new_pdrs4all_k.txt' )  
pdrs_F1130W  = ascii.read(inpath + 'PDRs4All_F1130W_D21_models_with_k_consub_new_pdrs4all_k.txt' )  

plt.scatter(pdrs_F335M['pah_k_method'][0:4]/pdrs_F770W['pah_k_method'][0:4], pdrs_F335M['pah_k_method'][0:4]/pdrs_F1130W['pah_k_method'][0:4], 
            s = 150, marker = '^', alpha = 0.5, zorder = 1000, c = 'tomato', label = 'Orion Bar (PDRs4All Team)')


ax.set_xlabel('3.3 $\mathrm{\mu}$m/7.7 $\mathrm{\mu}$m PAH', size = 'x-large')
ax.set_ylabel('3.3 $\mathrm{\mu}$m/11.3 $\mathrm{\mu}$m PAH', size = 'x-large')
ax.minorticks_on()

ax.semilogx()
ax.semilogy()

ax.set_xlim(0.03, 0.92)
ax.set_ylim(0.03, 0.22)

# --- Format ticks ---
# Major ticks as numbers, minor ticks hidden (to avoid overlap)
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(plt.NullFormatter())
ax.yaxis.set_minor_formatter(plt.NullFormatter())

# --- Control locations ---
# X-axis: fewer major ticks (prevents overlap)
# ax.xaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0, 2.0, 5.0], numticks=5))
ax.set_xticks([0.05, 0.1, 0.2, 0.9])


# Y-axis: more major ticks
# ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=np.arange(1.0, 10.0)*0.1, numticks=10))
# ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=np.arange(1.0, 10.0)*0.1, numticks=10))
ax.set_yticks([0.05, 0.1, 0.15, 0.2])


ax = plt.gca()
ax.grid(True, which='major', axis='y', linestyle='-', linewidth=0.2, color='gray', alpha=0.4)
ax.grid(True, which='minor', axis='y', linestyle='-', linewidth=0.2, color='gray', alpha=0.4)

ax.grid(True, which='major', axis='x', linestyle='-', linewidth=0.2, color='gray', alpha=0.4)
ax.grid(True, which='minor', axis='x', linestyle='-', linewidth=0.2, color='gray', alpha=0.4)

plt.legend(loc = 'upper left')


savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/'
savename = '{:s}D21_clumps_SigmaPAH_band_ratios_PDRs4All_k_D21_Halpha_K08.png'.format(D21_sp)
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 300)

# plot individual ratios against properties
plt.figure(6, figsize = (10,10))
plt.clf()
fig, axs = plt.subplots(nrows=3, ncols=1, num = 6, sharex = True)

# prop = np.log10(uv_data[uv_y, uv_x])
prop = ha_data[ha_y, ha_x]/hi_data[hi_y, hi_x]


axs[0].scatter(prop, pah3/pah7)
axs[0].set_ylabel('3.3/7.7')
for i in range(len(clumps)):
    axs[0].text(prop[i], pah3[i]/pah7[i], s = clumps['clump_num'][i])
    
res = stats.spearmanr(prop, pah3/pah7)    
print('3.3/7.7 sp coeff:', res.statistic)


axs[1].scatter(prop, pah3/pah11)
axs[1].set_ylabel('3.3/11.3')
for i in range(len(clumps)):
    axs[1].text(prop[i], pah3[i]/pah11[i], s = clumps['clump_num'][i])
    
res = stats.spearmanr(prop, pah3/pah11)    
print('3.3/11.3 sp coeff:', res.statistic)

axs[2].scatter(prop, pah7/pah11)
axs[2].set_ylabel('7.7/11.3')
for i in range(len(clumps)):
    axs[2].text(prop[i], pah7[i]/pah11[i], s = clumps['clump_num'][i])

plt.xlabel('H$\\alpha$/HI')
# plt.xlabel('log(UV)')

res = stats.spearmanr(prop, pah7/pah11)    
print('7.7/11.3 sp coeff:', res.statistic)

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/'
savename = 'band_ratios_prop_corr_Halpha_K08_HI.png'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 300)

# quick sigma PAH check 
sigma_pah_F1500W = (pah3 + pah7 + pah11)/clumps['F1500W_flux']

plt.figure(4)
plt.clf()
plt.scatter(prop, sigma_pah_F1500W)
plt.ylabel('$\Sigma$PAH/F1500W')

for i in range(len(clumps)):
    plt.text(prop[i], sigma_pah_F1500W[i], s = clumps['clump_num'][i])

plt.xlabel('log(H$\\alpha$)')
# plt.xlabel('log(UV)')

res = stats.spearmanr(prop, sigma_pah_F1500W)    
print('\$Sigma$PAH/F1500W sp coeff:', res.statistic)


plt.figure(5)
plt.clf()
plt.scatter(prop, clumps['F1500W_flux'])
plt.ylabel('F1500W')

for i in range(len(clumps)):
    plt.text(prop[i], clumps['F1500W_flux'][i], s = clumps['clump_num'][i])

plt.xlabel('log(H$\\alpha$)')
# plt.xlabel('log(UV)')

res = stats.spearmanr(prop, clumps['F1500W_flux'])    
print('F1500W sp coeff:', res.statistic)



