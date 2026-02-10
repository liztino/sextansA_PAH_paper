#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plots the results of calc_PAH_contamination.py

@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import seaborn as sns
import matplotlib.lines as mlines

plt.ion()

inpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/calc_PAH_contamination/'
final_table_name = 'pah_consub_calc_contamination_v1.txt'
data = ascii.read(inpath + final_table_name )  

# load the empirical table 
pdrs_table_name = 'PDRs4ALL_DF1_spec_PAH_contamination_v1.txt'
pdrs = ascii.read(inpath + pdrs_table_name )  

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/comparison/'

colors = sns.color_palette("hls", 3)

# translate size text into actual size 
d = {'sma': 20, 'lrg': 80, 'std': 40}
size_arr = [d[x] for x in data['size']] 

# translate ionization into colors
d = {'hi': colors[0], 'lo': colors[2], 'st': colors[1]}
ion_arr = [d[x] for x in data['ion']] 

# translate ionization into symbols
d = {'hi': '^', 'st': 'o', 'lo': 's'}
shape_arr = [d[x] for x in data['ion']] 


pfit = np.polyfit(data['F335M'], data['F360M'], deg = 1)
xx = np.linspace(0, 0.06, 1000)
p = np.poly1d(pfit)

plt.figure(1)
plt.clf()
im = plt.scatter(data['F335M'], data['F360M'], c = data['U'], cmap = 'magma', alpha = 0.7, edgecolor = 'k', s = size_arr)
plt.plot(xx, p(xx), c = 'k', label = 'Linear\na={:4.3e}\nb={:4.3e}'.format(pfit[0], pfit[1]))
cb = plt.colorbar(im)
cb.set_label('U')
plt.xlabel('F335M')
plt.ylabel('F360M')
plt.legend(loc='best')

savename = 'PAH_contamination_F335M_F360M_linear_fit.pdf'
# plt.savefig(savepath + savename)


pfit = np.polyfit(data['F770W'], data['F560W'], deg = 1)
xx = np.linspace(0, 0.04, 1000)
p = np.poly1d(pfit)

plt.figure(2)
plt.clf()
im = plt.scatter(data['F770W'], data['F560W'], c = data['U'], cmap = 'magma', alpha = 0.7, edgecolor = 'k', s = size_arr)
plt.plot(xx, p(xx), c = 'k', label = 'Linear\na={:4.3e}\nb={:4.3e}'.format(pfit[0], pfit[1]))
cb = plt.colorbar(im)
cb.set_label('U')
plt.xlabel('F770W')
plt.ylabel('F560W')
plt.legend(loc='best')

savename = 'PAH_contamination_F770W_F560W_linear_fit.pdf'
# plt.savefig(savepath + savename)

pfit = np.polyfit(np.log10(data['F335M']), np.log10(data['F360M']), deg = 1)
xx = np.logspace(min(np.log10(data['F335M'])), max(np.log10(data['F335M'])), 1000)
p = np.poly1d(pfit)

plt.figure(1)
plt.clf()
im = plt.scatter(np.log10(data['F335M']), np.log10(data['F360M']), c = data['U'], cmap = 'magma', alpha = 0.7, edgecolor = 'k', s = size_arr)
plt.plot(np.log10(xx), p(np.log10(xx)), c = 'k', label = 'Log\na={:4.3f}\nb={:4.3f}'.format(pfit[0], pfit[1]))
cb = plt.colorbar(im)
cb.set_label('U')
plt.xlabel('Log F335M')
plt.ylabel('Log F360M')
# plt.semilogx()
# plt.semilogy()
plt.legend(loc='best')

savename = 'PAH_contamination_F335M_F360M_log_fit.pdf'
# plt.savefig(savepath + savename)

pfit = np.polyfit(np.log10(data['F770W']), np.log10(data['F560W']), deg = 1)
xx = np.logspace(min(np.log10(data['F770W'])), max(np.log10(data['F770W'])), 1000)
p = np.poly1d(pfit)

plt.figure(2)
plt.clf()
im = plt.scatter(np.log10(data['F770W']), np.log10(data['F560W']), c = data['U'], cmap = 'magma', alpha = 0.7, edgecolor = 'k', s = size_arr)
plt.plot(np.log10(xx), p(np.log10(xx)), c = 'k', label = 'Log\na={:4.3f}\nb={:4.3f}'.format(pfit[0], pfit[1]))
cb = plt.colorbar(im)
cb.set_label('U')
plt.xlabel('Log F770W')
plt.ylabel('Log F560W')
# plt.semilogx()
# plt.semilogy()
plt.legend(loc='best')

savename = 'PAH_contamination_F770W_F560W_log_fit.pdf'
# plt.savefig(savepath + savename)


# find average value from U=0 to U=4
ind = np.where(data['U'] < 4)[0]
F360M_mean = np.nanmean(data['F360M'][ind]/data['F335M'][ind])
F560W_mean = np.nanmean(data['F770W_Jy'][ind]/data['F560W_Jy'][ind])

# now plot ratios vs U
plt.figure(3)
plt.clf()
im = plt.scatter(data['U'], data['F360M']/data['F335M'], s = size_arr, c = ion_arr, alpha = 0.7, edgecolor = 'k')
plt.axhline(pdrs['F360M'][0]/pdrs['F335M'][0], c= 'k', lw = 2)
plt.axhline(F360M_mean, c= 'blue', lw = 2)
plt.xlabel('U')
plt.ylabel('F360M/F335M')

hi = mlines.Line2D([], [], color=colors[0], marker='^', ls='', label='Hi', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
st = mlines.Line2D([], [], color=colors[1], marker='o', ls='', label='St', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
lo = mlines.Line2D([], [], color=colors[2], marker='s', ls='', label='Lo', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
pdrs_label = mlines.Line2D([], [], color='k',  ls='-', label='Empirical: PDRS4ALL', alpha = 1.0)
mean_label = mlines.Line2D([], [], color='blue',  ls='-', label='Mean U<4', alpha = 1.0)

# etc etc
plt.legend(handles=[mean_label, pdrs_label, hi, st, lo], loc = 'best')
plt.legend(handles=[mean_label, pdrs_label, hi, st, lo], loc = 'best')
plt.semilogy()


savename = 'PAH_contamination_F335M_F360M_U_with_PDRS.pdf'
# plt.savefig(savepath + savename)

plt.figure(4)
plt.clf()
for i in range(len(size_arr)):
    im = plt.scatter(data['U'][i], data['F770W_Jy'][i]/data['F560W_Jy'][i], s = size_arr[i], c = ion_arr[i], marker = shape_arr[i], alpha = 0.7, edgecolor = 'k')
plt.axhline(pdrs['F770W_Jy'][0]/pdrs['F560W_Jy'][0], c= 'k', lw = 2)
plt.axhline(F560W_mean, c= 'blue', lw = 2)
plt.xlabel('U', size = 'x-large')
plt.ylabel('$\mathrm{k_{F770W}= F770W_{PAH}/F560W_{PAH}}$', size = 'x-large')
plt.xlim(-0.2,4.2)

plt.legend(handles=[mean_label, pdrs_label, hi, st, lo], loc = 'best')

savename = 'PAH_contamination_F770W_F560W_U_with_PDRS_Jy.pdf'
plt.savefig(savepath + savename)

