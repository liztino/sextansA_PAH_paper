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

filt_contam = 'F335M'
filt_con = 'F360M'

inpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/calc_PAH_contamination/'
final_table_name = filt_contam + '_pah_consub_calc_contamination_v2.txt'
data = ascii.read(inpath + final_table_name )  


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

# find average value from U=0 to U=4
ind = np.where(data['U'] < 4)[0]
mean_ratio1 = np.nanmean(flux_con[ind]/flux_contam[ind])
mean_ratio2 = np.nanmean(flux_contam[ind]/flux_con[ind])

std_ratio1 = np.nanstd(flux_con[ind]/flux_contam[ind])
std_ratio2 = np.nanstd(flux_contam[ind]/flux_con[ind])


if filt_contam == 'F1500W':
    # now plot ratios vs U
    plt.figure(1)
    plt.clf()
    
    # results from D21
    for i in range(len(size_arr)):
        im = plt.scatter(data['U'][i], flux_con[i]/flux_contam[i], s = size_arr[i], c = ion_arr[i], marker = shape_arr[i], alpha = 0.7, edgecolor = 'k')
    plt.xlabel('U', size = 'x-large')
    plt.ylabel('$\mathrm{k_{F1130W}= F1130W_{PAH}/F1500W_{PAH}}$', size = 'x-large')
    plt.xlim(-0.2,4.2)
    plt.ylim(5,15)
    
    
    # bands from empirical results
    avg_k = 7.21
    std_k = 0.92
    plt.axhline(avg_k, c = 'cornflowerblue', lw = 1.5, ls = '--')
    ax = plt.gca()
    xlim = ax.get_xlim()
    x_vals = np.linspace(xlim[0], xlim[1], 10)
    plt.fill_between(x_vals, avg_k - std_k, avg_k + std_k, alpha = 0.1, color = 'cornflowerblue' )
    
    plt.axhline(mean_ratio1, c = 'k', lw = 1.5, ls = '--')
    plt.fill_between(x_vals, mean_ratio1 - std_ratio1, mean_ratio1 + std_ratio1, color = 'k', alpha = 0.1)
    
    
    
    
    # legend things
    hi = mlines.Line2D([], [], color=colors[0], marker='^', ls='', label='Hi', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
    st = mlines.Line2D([], [], color=colors[1], marker='o', ls='', label='St', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
    lo = mlines.Line2D([], [], color=colors[2], marker='s', ls='', label='Lo', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
    pdrs_label = mlines.Line2D([], [], color='k', ls='-', label='PDRs4All', alpha = 1.0)
    # plt.legend(handles=[mean_label, pdrs_label,  hi, st, lo], loc = 'best')
    # plt.legend(handles=[mean_label, pdrs_label, pdrs_label2, hi, st, lo], loc = 'best')
    mean_label = mlines.Line2D([], [], color='blue',  ls='-', label='D21 Models', alpha = 1.0)
    
    
    savename = 'PAH_contamination_{}_{}_U_Jy_final.pdf'.format(filt_contam, filt_con)
    plt.savefig(savepath + savename)
    
elif filt_contam == 'F770W':
    # now plot ratios vs U
    plt.figure(1)
    plt.clf()
    
    # results from D21
    for i in range(len(size_arr)):
        im = plt.scatter(data['U'][i], flux_contam[i]/flux_con[i], s = size_arr[i], c = ion_arr[i], marker = shape_arr[i], alpha = 0.7, edgecolor = 'k')
    plt.xlabel('U', size = 'x-large')
    plt.ylabel('$\mathrm{k_{F770W}= F770W_{PAH}/F560W_{PAH}}$', size = 'x-large')
    plt.xlim(-0.2,4.2)
    # plt.ylim(5,15)
    
    
    # bands from empirical results
    avg_k = 4.33
    std_k = 0.35
    plt.axhline(avg_k, c = 'tomato', lw = 1.5, ls = '--')
    ax = plt.gca()
    xlim = ax.get_xlim()
    x_vals = np.linspace(xlim[0], xlim[1], 10)
    plt.fill_between(x_vals, avg_k - std_k, avg_k + std_k, alpha = 0.1, color = 'tomato' )
    
    plt.axhline(mean_ratio2, c = 'k', lw = 1.5, ls = '--')
    plt.fill_between(x_vals, mean_ratio2 - std_ratio2, mean_ratio2 + std_ratio2, color = 'k', alpha = 0.1)
    
    
    
    
    # legend things
    hi = mlines.Line2D([], [], color=colors[0], marker='^', ls='', label='Hi', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
    st = mlines.Line2D([], [], color=colors[1], marker='o', ls='', label='St', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
    lo = mlines.Line2D([], [], color=colors[2], marker='s', ls='', label='Lo', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
    pdrs_label = mlines.Line2D([], [], color='k', ls='-', label='PDRs4All', alpha = 1.0)
    # plt.legend(handles=[mean_label, pdrs_label,  hi, st, lo], loc = 'best')
    # plt.legend(handles=[mean_label, pdrs_label, pdrs_label2, hi, st, lo], loc = 'best')
    mean_label = mlines.Line2D([], [], color='blue',  ls='-', label='D21 Models', alpha = 1.0)
    
    
    savename = 'PAH_contamination_{}_{}_U_Jy_final.pdf'.format(filt_contam, filt_con)
    plt.savefig(savepath + savename)
    
elif filt_contam == 'F335M':
    # now plot ratios vs U
    plt.figure(1)
    plt.clf()
    
    # results from D21
    for i in range(len(size_arr)):
        im = plt.scatter(data['U'][i], flux_contam[i]/flux_con[i], s = size_arr[i], c = ion_arr[i], marker = shape_arr[i], alpha = 0.7, edgecolor = 'k')
    plt.xlabel('U', size = 'x-large')
    plt.ylabel('$\mathrm{k_{F335M}= F335M_{PAH}/F360M_{PAH}}$', size = 'x-large')
    plt.xlim(-0.2,4.2)
    plt.ylim(0,31)
    
    
    # bands from empirical results
    avg_k = 4.45
    std_k = 0.39
    plt.axhline(avg_k, c = 'cornflowerblue', lw = 1.5, ls = '--')
    ax = plt.gca()
    xlim = ax.get_xlim()
    x_vals = np.linspace(xlim[0], xlim[1], 10)
    plt.fill_between(x_vals, avg_k - std_k, avg_k + std_k, alpha = 0.1, color = 'cornflowerblue' )
    
    avg_k = 2.07
    std_k = 0.3
    plt.axhline(avg_k, c = 'tomato', lw = 1.5, ls = '--')
    ax = plt.gca()
    xlim = ax.get_xlim()
    x_vals = np.linspace(xlim[0], xlim[1], 10)
    plt.fill_between(x_vals, avg_k - std_k, avg_k + std_k, alpha = 0.1, color = 'tomato' )
    
    plt.axhline(mean_ratio2, c = 'k', lw = 1.5, ls = '--')
    plt.fill_between(x_vals, mean_ratio2 - std_ratio2, mean_ratio2 + std_ratio2, color = 'k', alpha = 0.1)
    
    
    
    
    # legend things
    hi = mlines.Line2D([], [], color=colors[0], marker='^', ls='', label='Hi', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
    st = mlines.Line2D([], [], color=colors[1], marker='o', ls='', label='St', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
    lo = mlines.Line2D([], [], color=colors[2], marker='s', ls='', label='Lo', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
    # pdrs_label = mlines.Line2D([], [], color='k', ls='-', label='PDRs4All', alpha = 1.0)
    # plt.legend(handles=[mean_label, pdrs_label,  hi, st, lo], loc = 'best')
    # plt.legend(handles=[mean_label, pdrs_label, pdrs_label2, hi, st, lo], loc = 'best')
    # mean_label = mlines.Line2D([], [], color='blue',  ls='-', label='D21 Models', alpha = 1.0)
    
    plt.legend(handles = [hi, st, lo], loc = (0.8, 0.22), title = 'D21 Models')

    
    
    savename = 'PAH_contamination_{}_{}_U_Jy_final.pdf'.format(filt_contam, filt_con)
    plt.savefig(savepath + savename)    

# # now plot ratios vs U
# plt.figure(4)
# plt.clf()
# # im = plt.scatter(data['U'], flux_contam/flux_con, s = size_arr, c = ion_arr, alpha = 0.7, edgecolor = 'k')
# plt.axhline(np.mean(pdrs[filt_contam + '_Jy'][1:]/pdrs[filt_con + '_Jy'][1:]), c= 'k', lw = 2)
# for i in range(len(size_arr)):
#     im = plt.scatter(data['U'][i], flux_contam[i]/flux_con[i], s = size_arr[i], c = ion_arr[i], marker = shape_arr[i], alpha = 0.7, edgecolor = 'k')
# plt.axhline(mean_ratio2, c= 'blue', lw = 2, ls = '--')
# plt.xlabel('U', size = 'x-large')
# plt.ylabel('$\mathrm{k_{F335M}= F335M_{PAH}/F360M_{PAH}}$', size = 'x-large')

# print(filt_contam, mean_ratio2)
# if PDR4all:
#     print('PDRs4all',  np.mean(pdrs[filt_contam + '_Jy'][1:]/pdrs[filt_con + '_Jy'][1:]))
#     # print('PDRs4all',  np.mean(pahfit_contam[filt_contam + '_Jy']/pahfit_contam[filt_con + '_Jy']))

# hi = mlines.Line2D([], [], color=colors[0], marker='^', ls='', label='Hi', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
# st = mlines.Line2D([], [], color=colors[1], marker='o', ls='', label='St', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
# lo = mlines.Line2D([], [], color=colors[2], marker='s', ls='', label='Lo', alpha = 0.7, markeredgecolor  = 'k', ms = 7)
# pdrs_label = mlines.Line2D([], [], color='k', ls='-', label='PDRs4All', alpha = 1.0)
# # pdrs_label = mlines.Line2D([], [], color='b', marker='*', ls='', label='PDRS4ALL w/ PAHFIT', alpha = 1.0)

# mean_label = mlines.Line2D([], [], color='blue',  ls='--', label='D21 Models', alpha = 1.0)

# plt.xlim(-0.2, 4.3)
# # plt.ylim(1,  32)


# # ax = plt.gca()
# # ax2 = ax.twiny()

# # xlabel = ['T1', 'T2','T3','T4','T5']
# # if PDR4all:
# #     ax2.scatter(xlabel, pdrs[filt_contam + '_Jy']/pdrs[filt_con + '_Jy'], c= 'k', marker = '*', s = 50)
#     # ax2.scatter(xlabel, pahfit_contam[filt_contam + '_Jy']/pahfit_contam[filt_con + '_Jy'], c= 'blue', marker = '*', s = 50, alpha = 0.7)

    
#     # plt.ylim(0, 6)


# # plt.legend(handles=[mean_label, pdrs_label, hi, st, lo], loc = 'upper left')
# # plt.legend(handles=[mean_label, pdrs_label, pdrs_label2, hi, st, lo], loc = 'best')

# # plt.semilogy()

# if PDR4all:
#     savename = 'PAH_contamination_{}_{}_ratio2_U_Jy_pdrs4all_zoom_v2.pdf'.format(filt_con, filt_contam)
# else:
#     savename = 'PAH_contamination_{}_{}_ratio2_U_Jy_v2.pdf'.format(filt_con, filt_contam)
# # plt.savefig(savepath + savename)



