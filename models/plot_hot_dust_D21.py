#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plots the results from hot_dust_D21.py 

@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import scipy.interpolate as inter
import seaborn as sns


# qpah_list = [3.8, 1.0, 0.5, 0.2, 0.0]
# qpah_text_list = ['3p8', '1p0', '0p5', '0p2', '0p0']

# plt.figure(3)
# plt.clf()

# colors = sns.color_palette("hls", len(qpah_list))

# for i in range(len(qpah_list)):

#     savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/hot_dust_D21/'
#     filename = 'D21_ascii_hot_dust_F770W_F2100W_ratio_withpah_{:s}.txt'.format(qpah_text_list[i])
    
#     data = ascii.read(savepath + filename)
    
#     ratio = data['F770W']/data['F2100W']

#     plt.scatter(data['logU'], np.log10(ratio), s = 35, alpha = 0.7, color = colors[i], 
#                 label = '$\mathrm{{q_{{PAH}} = {:1.2f}}}$'.format(qpah_list[i]))

# plt.xlabel('Radiation Field, logU', size = 'x-large')
# plt.ylabel('Log(F770W/F2100W)', size = 'x-large')
# plt.minorticks_on()
# plt.legend(loc = 'best')

# savename = 'D21_models_F770W_F2100W_ratio_vs_logU.pdf'
# savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/hot_dust_D21/'
# plt.savefig(savepath + savename)

# plot qpah

filename = 'D21_ascii_hot_dust_F770W_F2100W_ratio_vary_q_pah.txt'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/hot_dust_D21/'

ratio_q_pah = ascii.read(savepath + filename)

plt.figure(4)
plt.clf()
plt.plot(ratio_q_pah['q_pah'], np.log10(ratio_q_pah['F770W']/ratio_q_pah['F2100W']), lw = 2, label = 'D21 Models')

sutter = (0.33 * ratio_q_pah['q_pah'] + 1.3)/2.57

plt.plot(ratio_q_pah['q_pah'], sutter, lw = 2, label = 'Sutter+24 Empirical')

plt.xlabel('$\mathrm{q_{PAH}}$', size = 'x-large')
plt.ylabel('Log(F770W/F2100W)', size = 'x-large')
plt.minorticks_on()
plt.title('logU = 0.0')
plt.legend(loc = 'best')

savename = 'D21_models_F770W_F2100W_ratio_vs_qpah.pdf'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/hot_dust_D21/'
plt.savefig(savepath + savename)