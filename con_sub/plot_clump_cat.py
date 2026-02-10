#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plots the results from the dendrogram PAH clumps

@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS

clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/clump_prop_cats/'

F335M_k_2 = ascii.read(clump_path + 'F335M_k_2.07_dendro_clump_props.txt')
F335M_k_4 = ascii.read(clump_path + 'F335M_k_4.45_dendro_clump_props.txt')
F1130W_k_10 = ascii.read(clump_path + 'F1130W_k_10.17_dendro_clump_props.txt')
F770W_k_4 = ascii.read(clump_path + 'F770W_k_4.33_dendro_clump_props.txt')
F770W_k_5 = ascii.read(clump_path + 'F770W_k_5.84_dendro_clump_props.txt')
F1130W_k_7 = ascii.read(clump_path + 'F1130W_k_7.21_dendro_clump_props.txt')
F1130W_k_10 = ascii.read(clump_path + 'F1130W_k_10.17_dendro_clump_props.txt')

pah3_sigma = 0.0036
pah7_sigma =  0.0066
pah11_sigma = 0.0138

savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/band_ratios/'
master_clump = ascii.read('/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/master_clump_cat.csv')

plt.figure(1)
plt.clf()
plt.scatter(master_clump['Clump'], F770W_k_4['flux'], label = 'F770W k=4.33', alpha = 0.5)
plt.scatter(master_clump['Clump'], F770W_k_5['flux'], label = 'F770W k=5.84', alpha = 0.5)
plt.legend(loc='best')
plt.xlabel('Clump number')
plt.ylabel('Flux')
plt.title('F770W')
savename = 'dendro_clump_F770W_compare.pdf'
plt.savefig(savedir + savename)

plt.figure(2)
plt.clf()
plt.scatter(master_clump['Clump'], F1130W_k_7['flux'], label = 'F1130W k=7.21', alpha = 0.5)
plt.scatter(master_clump['Clump'], F1130W_k_10['flux'], label = 'F1130W k=10.17', alpha = 0.5)
plt.legend(loc='best')
plt.xlabel('Clump number')
plt.ylabel('Flux')
plt.title('F1130W')
savename = 'dendro_clump_F1130W_compare.pdf'
plt.savefig(savedir + savename)

plt.figure(3)
plt.clf()
plt.scatter(master_clump['Clump'], F335M_k_2['flux'], label = 'F335M k=2', alpha = 0.5)
plt.scatter(master_clump['Clump'], F335M_k_4['flux'], label = 'F335M k=4', alpha = 0.5)
plt.legend(loc='best')
plt.xlabel('Clump number')
plt.ylabel('Flux')
plt.title('F335M')
savename = 'dendro_clump_F335M_compare.pdf'
plt.savefig(savedir + savename)


plt.figure(4)
plt.clf()
# plt.scatter(master_clump['Clump'], F770W_k_4['flux']/F1130W_k_7['flux'])
# plt.scatter(master_clump['Clump'], F770W_k_5['flux']/F1130W_k_10['flux'])
plt.scatter(master_clump['Clump'], F335M_k_2['flux']/F1130W_k_7['flux'], label = 'F335M k = 2/F1130W k = 7')
plt.scatter(master_clump['Clump'], F335M_k_4['flux']/F1130W_k_7['flux'], label = 'F335M k = 4/F1130W k = 7')
plt.legend(loc='best')
plt.xlabel('Clump number')
plt.ylabel('3.3/7.7')
savename = 'dendro_clump_F335M_F1130W_ratio_compare.pdf'
plt.savefig(savedir + savename)

err_x = np.sqrt(pah3_sigma**2 + pah7_sigma**2)
err_y = np.sqrt(pah3_sigma**2 + pah11_sigma**2)

# load the uv data to color points
anc_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
uv_name = 'sextans-a_m2_cr.fits'

uv_hdr = fits.open(anc_path + uv_name)
uv_data = uv_hdr[0].data
uv_head = uv_hdr[0].header
uv_w = WCS(uv_head)

# load in a continuum subtraction to get wcs
filt_mid = 'F770W'
k = 4.33
datapath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
dataname = 'SextansA_{:s}_k_method_pah_k_{:3.2f}.fits'.format(filt_mid, k)

pah_hdr = fits.open(datapath + dataname)
pah_head = pah_hdr[0].header
pah_w = WCS(pah_head)

# find world coordinates for all of the clumps
coord = pah_w.pixel_to_world(master_clump['x'], master_clump['y'])

# find pixel coordinates for the uv data
uv_x, uv_y = uv_w.world_to_pixel(coord)
uv_x = np.int32(uv_x)
uv_y = np.int32(uv_y)

plt.figure(5)
plt.clf()
plt.errorbar(F335M_k_2['flux']/F770W_k_4['flux'], F335M_k_2['flux']/F1130W_k_7['flux'],
             yerr = err_y, xerr = err_x, fmt = '.', ecolor = 'k', capsize = 2, zorder = 1, markersize = 0.01)
im = plt.scatter(F335M_k_2['flux']/F770W_k_4['flux'], F335M_k_2['flux']/F1130W_k_7['flux'], 
            c = np.log10(uv_data[uv_y, uv_x]), s = 100)
for i in range(len(master_clump)):
    plt.text(F335M_k_2['flux'][i]/F770W_k_4['flux'][i], F335M_k_2['flux'][i]/F1130W_k_7['flux'][i],
             s = '{:d}'.format(master_clump['Clump'][i]))
    
plt.minorticks_on()
plt.xlabel('3.3/7.7 PAH', size = 'x-large')
plt.ylabel('3.3/11.3 PAH', size = 'x-large')
plt.title('F335M k = 2, F770W k = 4, F1130W k = 7')
cb = plt.colorbar(im)
cb.set_label('log(UV)')
savename = 'dendro_clump_F335M_2_F770W_4_F1130W_7_ratio_UV_color.pdf'
plt.savefig(savedir + savename)

plt.figure(6)
plt.clf()
plt.errorbar(np.log10(uv_data[uv_y, uv_x]), F335M_k_2['flux']/F1130W_k_7['flux'],
             yerr = err_y, xerr = err_x, fmt = '.', ecolor = 'k', capsize = 2, zorder = 1, markersize = 0.01)
im = plt.scatter(np.log10(uv_data[uv_y, uv_x]), F335M_k_2['flux']/F1130W_k_7['flux'], 
            c = np.log10(uv_data[uv_y, uv_x]), s = 100)

plt.minorticks_on()
plt.title('F335M k = 2, F1130W k = 7')
plt.xlabel('log(UV)', size = 'x-large')
plt.ylabel('3.3/11.3 PAH', size = 'x-large')
cb = plt.colorbar(im)
cb.set_label('log(UV)')

savename = 'dendro_clump_F335M_2_F1130W_7_ratio_UV_color.pdf'
plt.savefig(savedir + savename)

plt.figure(7)
plt.clf()
plt.errorbar(np.log10(uv_data[uv_y, uv_x]), F770W_k_4['flux']/F1130W_k_7['flux'],
             yerr = err_y, xerr = err_x, fmt = '.', ecolor = 'k', capsize = 2, zorder = 1, markersize = 0.01)
im = plt.scatter(np.log10(uv_data[uv_y, uv_x]), F770W_k_4['flux']/F1130W_k_7['flux'], 
            c = np.log10(uv_data[uv_y, uv_x]), s = 100)

plt.minorticks_on()
plt.title('F770W k = 4, F1130W k = 7')
plt.xlabel('log(UV)', size = 'x-large')
plt.ylabel('7.7/11.3 PAH', size = 'x-large')
cb = plt.colorbar(im)
cb.set_label('log(UV)')
savename = 'dendro_clump_F770W_F1130W_7_ratio_UV_color.pdf'
plt.savefig(savedir + savename)


plt.figure(8)
plt.clf()
plt.errorbar(np.log10(uv_data[uv_y, uv_x]), F770W_k_4['flux']/F1130W_k_7['flux'],
             yerr = err_y, xerr = err_x, fmt = '.', ecolor = 'k', capsize = 2, zorder = 1, markersize = 0.01)
im = plt.scatter(np.log10(uv_data[uv_y, uv_x]), F770W_k_4['flux']/F1130W_k_7['flux'], 
            c = np.log10(uv_data[uv_y, uv_x]), s = 100)

plt.minorticks_on()
plt.xlabel('log(UV)', size = 'x-large')
plt.ylabel('7.7/11.3 PAH', size = 'x-large')
cb = plt.colorbar(im)
cb.set_label('log(UV)')

tot_PAH = F335M_k_2['flux'] + F770W_k_4['flux'] + F1130W_k_7['flux'] 
plt.figure(9)
plt.clf()

plt.scatter(master_clump['Clump'], F335M_k_2['flux']/tot_PAH, marker = '*', c = 'red', label = '3.3')
plt.scatter(master_clump['Clump'], F770W_k_4['flux']/tot_PAH, marker = 'o', c = 'orange', label = '7.7')
plt.scatter(master_clump['Clump'], F1130W_k_7['flux']/tot_PAH, marker = 's', c = 'purple', label = '11.3')
plt.legend(loc = 'best')
plt.ylabel('PAH/$\Sigma$PAH')
plt.xlabel('Clump')
savename = 'dendro_clump_PAH_SigmaPAH.pdf'
plt.savefig(savedir + savename)


