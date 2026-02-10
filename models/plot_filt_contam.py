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

filt_contam = 'F1500W'
filt_con = 'F1130W'

inpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/calc_PAH_contamination/'
final_table_name = filt_contam + '_pah_consub_calc_contamination_v2.txt'
data = ascii.read(inpath + final_table_name )  



# load the empirical table 
PDR4all = True

if PDR4all:

    pdrs_table_name = filt_contam + '_PDRs4ALL_spec_PAH_contamination_v2.txt'
    pdrs = ascii.read(inpath + pdrs_table_name )  
    
    pdrs_table_name2 = filt_contam + '_PDRs4ALL_spec_PAH_contamination_v5_emission_line.txt'
    pdrs2 = ascii.read(inpath + pdrs_table_name2 )  
    
    # final_table_name = 'PAHFIT_' + filt_contam + '_PDRs4ALL_spec_PAH_contamination_v2.txt'
    # pahfit_contam = ascii.read(inpath + final_table_name )  

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/plots/calc_PAH_contamination/comparison_filt/'

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


flux_contam = data[filt_contam + '_Jy']
flux_con = data[filt_con + '_Jy']


pfit = np.polyfit(flux_contam, flux_con, deg = 1)
xx = np.linspace(0, 0.06, 1000)
p = np.poly1d(pfit)

plt.figure(1)
plt.clf()
im = plt.scatter(flux_contam, flux_con, c = data['U'], cmap = 'magma', alpha = 0.7, edgecolor = 'k', s = size_arr)
plt.plot(xx, p(xx), c = 'k', label = 'Linear\na={:4.3e}\nb={:4.3e}'.format(pfit[0], pfit[1]))
cb = plt.colorbar(im)
cb.set_label('U')
plt.xlabel(filt_contam)
plt.ylabel(filt_con)
plt.legend(loc='best')

savename = 'PAH_contamination_{}_{}_linear_fit_Jy.pdf'.format(filt_contam, filt_con)
# plt.savefig(savepath + savename)


pfit = np.polyfit(np.log10(flux_contam), np.log10(flux_con), deg = 1)
xx = np.logspace(min(np.log10(flux_contam)), max(np.log10(flux_con)), 1000)
p = np.poly1d(pfit)

plt.figure(2)
plt.clf()
im = plt.scatter(np.log10(flux_contam), np.log10(flux_con), c = data['U'], cmap = 'magma', alpha = 0.7, edgecolor = 'k', s = size_arr)
plt.plot(np.log10(xx), p(np.log10(xx)), c = 'k', label = 'Log\na={:4.3f}\nb={:4.3f}'.format(pfit[0], pfit[1]))
cb = plt.colorbar(im)
cb.set_label('U')
plt.xlabel('Log ' + filt_contam)
plt.ylabel('Log ' + filt_con)
# plt.semilogx()
# plt.semilogy()
plt.legend(loc='best')

savename = 'PAH_contamination_{}_{}_log_fit_Jy.pdf'.format(filt_contam, filt_con)
# plt.savefig(savepath + savename)


# find average value from U=0 to U=4
ind = np.where(data['U'] < 4)[0]
mean_ratio1 = np.nanmean(flux_con[ind]/flux_contam[ind])
mean_ratio2 = np.nanmean(flux_contam[ind]/flux_con[ind])


# now plot ratios vs U
plt.figure(3)
plt.clf()
# im = plt.scatter(data['U'], flux_con/flux_contam, s = size_arr, c = ion_arr, alpha = 0.7, edgecolor = 'k')
plt.clf()
# plt.axhline(np.mean(pdrs[filt_con + '_Jy'][1:]/pdrs[filt_contam + '_Jy'][1:]), c= 'k', lw = 2)

for i in range(len(size_arr)):
    im = plt.scatter(data['U'][i], flux_con[i]/flux_contam[i], s = size_arr[i], c = ion_arr[i], marker = shape_arr[i], alpha = 0.7, edgecolor = 'k')
# plt.axhline(pdrs['F770W_Jy'][0]/pdrs['F560W_Jy'][0], c= 'k', lw = 2)
plt.axhline(mean_ratio1, c= 'blue', lw = 2, ls = '--')
plt.xlabel('U', size = 'x-large')
# plt.ylabel('$\mathrm{k_{F770W}= F770W_{PAH}/F560W_{PAH}}$', size = 'x-large')
plt.xlim(-0.2,4.2)
plt.ylim(5,15)

hi = mlines.Line2D([], [], color=colors[0], marker='^', ls='', label='Hi', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
st = mlines.Line2D([], [], color=colors[1], marker='o', ls='', label='St', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
lo = mlines.Line2D([], [], color=colors[2], marker='s', ls='', label='Lo', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
pdrs_label = mlines.Line2D([], [], color='k', ls='-', label='PDRs4All', alpha = 1.0)
plt.ylabel('$\mathrm{k_{F1130W}= F1130W_{PAH}/F1500W_{PAH}}$', size = 'x-large')

avg_k = 7.21
std_k = 0.92
plt.axhline(avg_k, c = 'k', lw = 1.0)
ax = plt.gca()
xlim = ax.get_xlim()
x_vals = np.arange(xlim[0], xlim[1], 10)
plt.fill_between(x_vals, avg_k - std_k, avg_k + std_k)

mean_label = mlines.Line2D([], [], color='blue',  ls='-', label='D21 Models', alpha = 1.0)

ax = plt.gca()
ax2 = ax.twiny()

xlabel = ['T1', 'T2','T3','T4','T5']
if PDR4all:
    ax2.scatter(xlabel, pdrs[filt_con + '_Jy']/pdrs[filt_contam + '_Jy'], c= 'k', marker = '*', s = 50)
    ax2.scatter(xlabel, pdrs2[filt_con + '_Jy']/pdrs2[filt_contam + '_Jy'], c= 'red', marker = '*', s = 50)

    # ax2.scatter(xlabel, pahfit_contam[filt_con + '_Jy']/pahfit_contam[filt_contam + '_Jy'], c= 'blue', marker = '*', s = 50, alpha = 0.7)

    # plt.ylim(8, 26)

# etc etc
# plt.legend(handles=[mean_label, pdrs_label,  hi, st, lo], loc = 'best')
# plt.legend(handles=[mean_label, pdrs_label, pdrs_label2, hi, st, lo], loc = 'best')

# plt.semilogy()

# plt.ylim(4, 15)

if PDR4all: 
    savename = 'PAH_contamination_{}_{}_U_Jy_final.pdf'.format(filt_contam, filt_con)
else:
    savename = 'PAH_contamination_{}_{}_U_Jy.pdf'.format(filt_contam, filt_con)
plt.savefig(savepath + savename)

# now plot ratios vs U
plt.figure(4)
plt.clf()
# im = plt.scatter(data['U'], flux_contam/flux_con, s = size_arr, c = ion_arr, alpha = 0.7, edgecolor = 'k')
plt.axhline(np.mean(pdrs[filt_contam + '_Jy'][1:]/pdrs[filt_con + '_Jy'][1:]), c= 'k', lw = 2)
for i in range(len(size_arr)):
    im = plt.scatter(data['U'][i], flux_contam[i]/flux_con[i], s = size_arr[i], c = ion_arr[i], marker = shape_arr[i], alpha = 0.7, edgecolor = 'k')
plt.axhline(mean_ratio2, c= 'blue', lw = 2, ls = '--')
plt.xlabel('U', size = 'x-large')
plt.ylabel('$\mathrm{k_{F335M}= F335M_{PAH}/F360M_{PAH}}$', size = 'x-large')

print(filt_contam, mean_ratio2)
if PDR4all:
    print('PDRs4all',  np.mean(pdrs[filt_contam + '_Jy'][1:]/pdrs[filt_con + '_Jy'][1:]))
    # print('PDRs4all',  np.mean(pahfit_contam[filt_contam + '_Jy']/pahfit_contam[filt_con + '_Jy']))

hi = mlines.Line2D([], [], color=colors[0], marker='^', ls='', label='Hi', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
st = mlines.Line2D([], [], color=colors[1], marker='o', ls='', label='St', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
lo = mlines.Line2D([], [], color=colors[2], marker='s', ls='', label='Lo', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
pdrs_label = mlines.Line2D([], [], color='k', ls='-', label='PDRs4All', alpha = 1.0)
# pdrs_label = mlines.Line2D([], [], color='b', marker='*', ls='', label='PDRS4ALL w/ PAHFIT', alpha = 1.0)

mean_label = mlines.Line2D([], [], color='blue',  ls='--', label='D21 Models', alpha = 1.0)

plt.xlim(-0.2, 4.3)
# plt.ylim(1,  32)


# ax = plt.gca()
# ax2 = ax.twiny()

# xlabel = ['T1', 'T2','T3','T4','T5']
# if PDR4all:
#     ax2.scatter(xlabel, pdrs[filt_contam + '_Jy']/pdrs[filt_con + '_Jy'], c= 'k', marker = '*', s = 50)
    # ax2.scatter(xlabel, pahfit_contam[filt_contam + '_Jy']/pahfit_contam[filt_con + '_Jy'], c= 'blue', marker = '*', s = 50, alpha = 0.7)

    
    # plt.ylim(0, 6)


# plt.legend(handles=[mean_label, pdrs_label, hi, st, lo], loc = 'upper left')
# plt.legend(handles=[mean_label, pdrs_label, pdrs_label2, hi, st, lo], loc = 'best')

# plt.semilogy()

if PDR4all:
    savename = 'PAH_contamination_{}_{}_ratio2_U_Jy_pdrs4all_zoom_v2.pdf'.format(filt_con, filt_contam)
else:
    savename = 'PAH_contamination_{}_{}_ratio2_U_Jy_v2.pdf'.format(filt_con, filt_contam)
# plt.savefig(savepath + savename)



