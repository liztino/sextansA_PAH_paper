#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 14:14:29 2023

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from astropy.visualization import make_lupton_rgb
from regions.core import PixCoord
from regions.shapes.rectangle import RectanglePixelRegion
from astropy import wcs
from regions import Regions
import matplotlib as mpl
from astropy.modeling import models, fitting
from astropy import units as u
import seaborn as sns

os.environ['PYSYN_CDBS']
import pysynphot as S

plt.ion()

# constants
c = 3e14        # in microns/sec
h = 6.62e-27    # in erg s 

# data paths
nircam_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/nircam_may/'
miri_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/miri_may/mosaic/'

# analysis paths
reg_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/SEDs/II6_reg_tests/reg/'
savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/SEDs/II6_reg_tests/SED/'

reg_name = 'pah335_all_small'

temp = 4000
logg = 4.0 
Z = -0.2
sp = S.Icat('phoenix', temp,Z, logg)

savepath = savedir + reg_name + '/'
if not os.path.isdir(savepath): 
    os.mkdir(savepath)

filt_list = ['sp_24um','F1500W', 'F1130W', 'F1000W', 'F560W', 'F770W', 'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M']

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
spitzer_filt = {'sp_24um'}

reg_list = Regions.read(reg_path + reg_name + '.reg')

micron_low = np.zeros(len(filt_list))
micron_high = np.zeros(len(filt_list))
flux = np.zeros(len(filt_list))
avg = np.zeros(len(filt_list))
std = np.zeros(len(filt_list))

for j in range(len(reg_list)):
# for j in range(1):
    reg = reg_list[j]
    c_list = []
    reg_name = reg.meta['text']

    for i, filt in enumerate(filt_list):
        # check if it's a miri or nircam filter
        if filt in miri_filts:
            if filt == 'F1130W':
                infile = miri_path + f'miri_{filt}_stage3_asn_skysub_i2d.fits'
                
            else:
                infile = miri_path + f'miri_{filt}_stage3_asn_fixed_wcs_skysub_i2d.fits'
                
            filt_file = 'JWST_MIRI.{:s}.dat'.format(filt)
            filter_path = '/astro/dust_kg/etarantino/JWST/filts/'
            spitzer = False


        elif filt in nircam_filts:
            filt_name = filt.lower()
            infile = nircam_path + filt + '/' + f'jw02391007001_nircam_{filt_name}_asn_modsub_i2d.fits'
            filt_file = 'JWST_NIRCam.{:s}.dat'.format(filt)
            filter_path = '/astro/dust_kg/etarantino/JWST/filts/'
            spitzer = False

        elif filt in spitzer_filt:
            infile = '/user/etarantino/SEXA_ANCILLARY/Sextans_A I MIPS24 d2009.fits'
            filter_path = '/user/etarantino/SEXA_ANCILLARY/filts/'
            filt_file = 'Spitzer_MIPS.24mu.dat'
            spitzer = True
                    
        # open fits file
        hdu = fits.open(infile)

        if spitzer: 
            header = hdu[0].header
            data = hdu[0].data
        else:   
            header = hdu['SCI'].header
            data = hdu['SCI'].data
        
        # open filter file 
        # wavelength is in Angstroms, convert to microns
        filt_curve = ascii.read(filter_path + filt_file, names = ['wave', 'trans'])
        filt_curve['wave'] = filt_curve['wave']/1e4
        
        # # get the middle of the filter for plotting purposes
        # micron[i] = ((filt_curve['wave'][-1] - filt_curve['wave'][0])/2) + filt_curve['wave'][0]
        # print(micron[i])
        
        # get position of filter where it is half it's max transmission to plot
        cutoff = np.nanmax(filt_curve['trans'])*0.5
        ind = np.where(filt_curve['trans'] >  cutoff)[0]
        ind1 = ind[0]
        ind2 = ind[-1]
        micron_low[i] = filt_curve['wave'][ind1]
        micron_high[i] = filt_curve['wave'][ind2]
    
        
        # convert from pixel/degree to pixel/sr
        w = wcs.WCS(header)
        pix_area = wcs.utils.proj_plane_pixel_area(w) * 3600**2 / (4.25e10)
    
        #converts region to frame of data and creates a mask 
        pixel_reg = reg.to_pixel(w)
        pix_mask = pixel_reg.to_mask()
        shape = np.shape(data)
        mask = pix_mask.to_image(shape)
        final_mask = np.array(mask, dtype=bool)
        
        # add up the flux and multiply by 
        flux[i] = np.nansum(data[final_mask]) * pix_area
        avg[i] = np.nanmean(data[final_mask])
        std[i] = np.nanstd(data[final_mask])
        
        print(flux[i], avg[i], std[i])
            
        # # test to make sure region is correct 
        # test_data = np.zeros(np.shape(data))
        # test_data[final_mask] = data[final_mask]
        # plt.figure()
        # plt.imshow(test_data, origin = 'lower')
        # plt.show()    
            
                
        if filt == 'F335M':
            c_list.append('magenta')
        elif filt == 'F770W':
            c_list.append('cyan')
        elif filt == 'F1130W':
            c_list.append('orange')
        elif filt == 'sp_24um':
            c_list.append('gray')
        else:
            c_list.append('black')
            
    
    plt.figure(j, figsize = (10,5))
    plt.clf()
    plt.hlines(flux, micron_low, micron_high, lw = 3.0, colors = c_list, zorder = 10000)
    plt.title('Region ' + str(reg_name))
    plt.xlabel('Wavelength ($\mu$m)', size = 'large')
    plt.ylabel('Flux Sum (MJy)', size = 'large')
    plt.minorticks_on()
    
    # plot draine model 
    drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
    model = 'BC03_Z0.0004_10Myr'
    
    drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
    header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']

    data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)
    
    # convert from nu * P_nu to P_nu
    Hz = c/data['wave']
    draine_flux =  data['total']*Hz
    
    # # match the draine model to the value at the F1500W filter
    # ind_filt = np.where(micron_high>15)[0][0]
    # ind_draine = np.where(data['wave']>15)[0][0]

    # flux_ratio = flux[ind_filt]/draine_flux[ind_draine]
    
    # flux_norm = draine_flux*flux_ratio
    
    norm_val = np.nanmax(flux)
    flux_norm = draine_flux/np.nanmax(draine_flux) * norm_val
    
    plt.plot(data['wave'], flux_norm, c = 'gray', alpha = 0.5, lw = 1.5, zorder = 1, label = 'D21 PAH st sma')
        
    mic_mid = micron_low + (micron_high - micron_low)/2
    
    fit_dust = False
    if fit_dust:
    
        bb = models.BlackBody(temperature = 100*u.K, scale = 0.05)
    
        fit = fitting.LevMarLSQFitter()
        
        wave_Hz = c/mic_mid[0:4]
        fit_bb = fit(bb, wave_Hz, flux[0:4])
        print(fit_bb)
        
        wave_vals = np.linspace(0.25,27, 100)
        wave_vals_units = c/wave_vals
        
        bb_vals = fit_bb(wave_vals_units)
        
        plt.plot(wave_vals, bb_vals, c = 'r', alpha = 0.5, lw = 1.5, label = 'Dust T={:4.1f} K'.format(fit_bb.temperature.value))
    
    else:
        temp_array = np.array([200, 300, 500])
        scale_array = np.array([0.002, 0.0008, 0.0002])
        bb_array = []
        
        wave_vals = np.linspace(3,20, 100)
        wave_vals_units = c/wave_vals
        
        bb_sum = np.zeros(len(wave_vals))
        
        c_list = sns.color_palette("Reds", len(temp_array))
    
        for i in range(len(temp_array)):
            bb_mod =  models.BlackBody(temperature = temp_array[i]*u.K, scale =scale_array[i])
            bb_array.append(bb_mod)
        
            bb_vals = bb_mod(wave_vals_units)
            
            bb_sum = np.nansum([bb_sum, bb_vals], axis = 0)
            
            # plt.plot(wave_vals, bb_vals,  alpha = 0.5, lw = 1.5, label = 'Dust T={:4.1f} K'.format(bb_mod.temperature.value))
            plt.plot(wave_vals, bb_vals,  alpha = 0.5, lw = 1.5, c = c_list[i])

            
        plt.plot(wave_vals, bb_sum, alpha = 0.7, lw = 2, label = 'Mixed Dust', color = 'red')
        
        
    
    # plot phoenix models
    ph_wave = sp.wave/1e4
    ind_filt = np.where(micron_high<1.3)[0][0]
    ind_ph = np.where(ph_wave>1.15)[0][0]

    flux_ratio = flux[ind_filt]/sp.flux[ind_ph]
    ph_norm = sp.flux * flux_ratio
    
    plt.plot(sp.wave/1e4, ph_norm, label = 'Phoenix T={:5.0f} log g={:3.1f}'.format(temp, logg), alpha = 0.5, c= '#1f77b4')
    
    max_flux = np.nanmax(flux)
    plt.ylim(0.0, max_flux+ max_flux*0.1)
    
    plt.legend(loc='upper right')
    
    plt.xlim(0.25,20)

    
    plt.show()
    
    savename = f'sexa_rough_SED_region_{reg_name}_flux_sum_F1130W_fit_dust_mods.pdf'
    plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)
    
    
    # ##### plot the mean and standard deviation #####
    
    # plt.figure(j, figsize = (10,5))
    # plt.clf()
    # plt.hlines(avg, micron_low, micron_high, lw = 3.0, colors = c_list, zorder = 10000)
    # ax = plt.gca()
    # xlim = ax.get_xlim()
    # ylim = ax.get_ylim()

    # micron = micron_low + (micron_high - micron_low)/2
    # plt.vlines(micron, ymin = avg-std, ymax = avg+std, lw = 1.0, zorder = 100, colors = c_list)
    # plt.title('Region ' + str(reg_name))
    # plt.xlabel('Wavelength ($\mu$m)', size = 'large')
    # plt.ylabel('Flux Average (MJy/sr)', size = 'large')
    # plt.minorticks_on()

    
    # # plot draine model 
    # drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
    # model = 'BC03_Z0.0004_10Myr'
    
    # drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
    # header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']

    # data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)
    
    # # convert from nu * P_nu to P_nu
    # Hz = c/data['wave']
    # draine_flux =  data['total']*Hz
    
    # # # match the draine model to the value at the F1500W filter
    # # ind_filt = np.where(micron_high>15)[0][0]
    # # ind_draine = np.where(data['wave']>15)[0][0]

    # # flux_ratio = flux[ind_filt]/draine_flux[ind_draine]
    
    # # flux_norm = draine_flux*flux_ratio
    
    # norm_val = np.nanmax(avg)
    # flux_norm = draine_flux/np.nanmax(draine_flux) * norm_val
    
    # plt.plot(data['wave'], flux_norm, c = 'gray', alpha = 0.5, lw = 1.5, zorder = 1)
    
    # plt.xlim(0.25,17)
    # plt.ylim(0, ylim[1])
    
    # plt.show()
    
    # savename = f'sexa_rough_SED_region_{reg_name}_avg.pdf'
    # # plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)
        
    
    
    
    
    
        
