#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:14:29 2023

@author: etarantino
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

filepath = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/nircam_mar/'
savepath = '/astro/dust_kg/etarantino/JWST_PAHS_2391/plots/modsub_test/'
    
short_filts = {'f115w', 'f150w', 'f200w'}
long_filts = {'f300m', 'f335m', 'f360m'}
short_det = ['nrca1', 'nrca2', 'nrca3', 'nrca4', 'nrcb1', 'nrcb2', 'nrcb3', 'nrcb4']
long_det = ['nrcalong', 'nrcblong']

filt_list = ['f115w', 'f150w', 'f200w', 'f300m', 'f335m', 'f360m']
short_det = ['nrca1', 'nrca2', 'nrca3', 'nrca4', 'nrcb1', 'nrcb2', 'nrcb3', 'nrcb4']
long_det = ['nrcalong', 'nrcblong']



detector_plots = False
mosaic_plots = True

# if mosaic_plots:
#     for filt in filt_list:     
    
    
# for filt in filt_list:
#     if filt in short_filts:
#         det_array = short_det
#     elif filt in long_filts:
#         det_array = long_det
        
#     up_filt = filt.upper()
#     workdir =  filepath + up_filt + '/'
        
#     # det_array = ['nrca3']
#     for det in det_array:
#         for i in np.arange(1,5,1):
#             orig_file = glob.glob(workdir + f"*_0000{str(i)}_{det}_a3001_crf.fits")[0]
#             modsub_file = glob.glob(workdir + f"*_0000{str(i)}_{det}_modsub_a3001_crf.fits")[0]
#             # print(orig_file)
#             # print(modsub_file)
            
#             orig = fits.open(orig_file)
#             orig_data = orig['SCI'].data
            
#             modsub = fits.open(modsub_file)
#             modsub_data = modsub['SCI'].data
            
#             shape_orig = np.shape(orig_data)
#             shape_modsub = np.shape(modsub_data)
            
#             # print(shape_orig, shape_modsub)
#             if not (shape_orig == shape_modsub):
#                 print(filt + ' ' + det)
            
        
    
filt = 'f150w'
up_filt = filt.upper()
workdir =  filepath + up_filt + '/'

orig_file = f'jw02391007001_nircam_{filt}_asn_i2d.fits'
# orig_file = 'jw02391007001_05101_00001_nrcb2_a3001_crf.fits'
orig = fits.open(workdir + orig_file)
orig_data = orig['SCI'].data

modsub_file = f'jw02391007001_nircam_{filt}_asn_modsub_i2d.fits'
# modsub_file = 'jw02391007001_05101_00001_nrcb2_modsub_a3001_crf.fits'
modsub = fits.open(workdir + modsub_file)
modsub_data = modsub['SCI'].data

fig = plt.figure(1, figsize = (15,10))
plt.clf()
gs = gridspec.GridSpec(2, 4)
gs.update(wspace=0.5)
ax1 = plt.subplot(gs[0, :2], )
ax2 = plt.subplot(gs[0, 2:])
ax3 = plt.subplot(gs[1, 1:3])
plt.show()

# fig.subplots_adjust(top=0.8)

ax_list = [ax1, ax2, ax3]

vmin = np.nanpercentile(orig_data, 1)
vmax = np.nanpercentile(orig_data, 99)
im0 = ax_list[0].imshow(orig_data, vmin = vmin, vmax = vmax, cmap = 'viridis', origin = 'lower')
ax_list[0].set_title('Orig')
ax_list[0].set_xlim(5310, 5485)
ax_list[0].set_ylim(1530, 1700)

divider = make_axes_locatable(ax_list[0])
cax = divider.append_axes('bottom', size='5%', pad=0.4)
fig.colorbar(im0, cax=cax, orientation='horizontal')

im1 = ax_list[1].imshow(modsub_data, vmin = vmin, vmax = vmax, cmap = 'viridis', origin = 'lower')
ax_list[1].set_title('Modsub')
ax_list[1].set_xlim(5310, 5485)
ax_list[1].set_ylim(1530, 1700)

divider = make_axes_locatable(ax_list[1])
cax = divider.append_axes('bottom', size='5%', pad=0.4)
fig.colorbar(im1, cax=cax, orientation='horizontal')

diff = orig_data - modsub_data

vmin = np.nanpercentile(diff, 1)
vmax = np.nanpercentile(diff, 99)
im2 = ax_list[2].imshow(diff, vmin = vmin, vmax = vmax, cmap = 'viridis', origin = 'lower')
ax_list[2].set_title('Orig - Modsub')

ax_list[2].set_xlim(5310, 5485)
ax_list[2].set_ylim(1530, 1700)

divider = make_axes_locatable(ax_list[2])
cax = divider.append_axes('bottom', size='5%', pad=0.4)
fig.colorbar(im2, cax=cax, orientation='horizontal')

fig.suptitle(f"Filter: {filt}", size = 'x-large', y = .95)


plt.draw()

    
savedir = savepath + up_filt + '/'

if not os.path.isdir(savedir): 
    os.mkdir(savedir)

savename = f'nircam_call_{filt}_modsub_mosaic_compare.pdf'

plt.savefig(savedir + savename, bbox_inches='tight',pad_inches = 0.1)

if detector_plots:
    
    for filt in filt_list:
        if filt in short_filts:
            det_array = short_det
        elif filt in long_filts:
            det_array = long_det
            
        # det_array = ['nrca3']
        for det in det_array:
            
            plt.figure(1, figsize = (15,5))
            plt.clf()
            fig, axs = plt.subplots(1,4, figsize = (15,5), num = 1)
            fig.subplots_adjust(top=0.8)
    
            ax_list = axs.flatten()
        
            up_filt = filt.upper()
            
            # pipeline working directory 
            workdir = filepath + up_filt + '/'
            
            sky_file = glob.glob(workdir + f"*{det}_sky_median.fits")[0]
            sky = fits.open(sky_file)
            sky_data = sky[0].data
            
            mod_file = glob.glob(workdir + f"*{det}_model_bkgr.fits")[0]
            mod = fits.open(mod_file)
            mod_data = mod[0].data
            
            sky_mod_diff = sky_data - mod_data
            
            vmin = np.nanpercentile(sky_data, 1)
            vmax = np.nanpercentile(sky_data, 99)
            im0 = ax_list[0].imshow(sky_data, vmin = vmin, vmax = vmax, cmap = 'viridis', origin = 'lower')
            ax_list[0].set_title('Background')
            
            divider = make_axes_locatable(ax_list[0])
            cax = divider.append_axes('bottom', size='5%', pad=0.4)
            fig.colorbar(im0, cax=cax, orientation='horizontal')
            
            vmin = np.nanpercentile(mod_data, 1)
            vmax = np.nanpercentile(mod_data, 99)
            im1 = ax_list[1].imshow(mod_data, vmin = vmin, vmax = vmax, cmap = 'magma', origin = 'lower')
            ax_list[1].set_title('Model')
            
            divider = make_axes_locatable(ax_list[1])
            cax = divider.append_axes('bottom', size='5%', pad=0.4)
            fig.colorbar(im1, cax=cax, orientation='horizontal')
            
            im2 = ax_list[2].imshow(sky_data, vmin = vmin, vmax = vmax, cmap = 'magma', origin = 'lower')
            ax_list[2].set_title('Background on Model Scale')
            
            divider = make_axes_locatable(ax_list[2])
            cax = divider.append_axes('bottom', size='5%', pad=0.4)
            fig.colorbar(im2, cax=cax, orientation='horizontal')
            
            vmax = np.nanpercentile(sky_mod_diff, 99)
            im3 = ax_list[3].imshow(sky_mod_diff, vmin = -vmax, vmax = vmax, cmap = 'seismic', origin = 'lower')
            ax_list[3].set_title('Difference')
            
            divider = make_axes_locatable(ax_list[3])
            cax = divider.append_axes('bottom', size='5%', pad=0.4)
            fig.colorbar(im3, cax=cax, orientation='horizontal')
            
            fig.suptitle(f"Filter: {filt} Detector: {det}", size = 'x-large', y = .95)
        
            plt.draw()
            
            savedir = savepath + up_filt + '/'
            
            if not os.path.isdir(savedir): 
                os.mkdir(savedir)
            
            savename = f'nircam_call_{filt}_{det}_modsub_compare.pdf'
            
            plt.savefig(savedir + savename, bbox_inches='tight',pad_inches = 0.1)
            

