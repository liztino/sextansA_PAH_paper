#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculates filter flux from the PAH clumps defined by the F770W mask

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.table import Table
import glob
import os
from astropy.wcs import WCS
from scipy.ndimage import binary_dilation, binary_fill_holes, generate_binary_structure, binary_erosion
from statistics import stdev

from get_pivot_wave import get_pivot_wave

# load pah flux files
# def load_PAH(filt_mid, k):
#     datapath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
#     dataname = 'SextansA_{:s}_k_method_pah_k_{:3.2f}.fits'.format(filt_mid, k)
#     pah_hdr = fits.open(datapath + dataname)[0]
#     pah_data = pah_hdr.data
#     pah_head = pah_hdr.header
    
#     errname = 'SextansA_{:s}_k_method_err_k_{:3.2f}.fits'.format(filt_mid, k)
#     err_hdr = fits.open(datapath + errname)[0]
#     err_data = err_hdr.data
#     err_head = err_hdr.header
    
#     conname = 'SextansA_{:s}_k_method_con_k_{:3.2f}.fits'.format(filt_mid, k)
#     con_hdr = fits.open(datapath + conname)[0]
#     con_data = con_hdr.data
#     con_head = con_hdr.header    
    
#     return pah_data, pah_head, err_data, con_data

# load master clump file
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_file = 'dendro_mask_clump_flux_with_filt.txt'

clumps = ascii.read(clump_path + clump_file)
clump_num = clumps['clump_num']

# constants 
pixel_area = (3.0555555555555E-05)**2

# define filters
miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}
lw_filts = {'F300M', 'F335M', 'F360M'}


# load filters
# slightly changing this script so it loads from the original data, not the PSF matched
def load_filter(filt):
    
    # original data path not PSF matched
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
    filename = f'{filt}_reproject_to_F1500W_rot.fits'
    
    # load error file
    err_file = 'err_reg_reproject_to_F1500W_rot.txt'
    err_filt = ascii.read(filepath + err_file)
    
    # load the filter and info
    hdu = fits.open(filepath + filename)
    header = hdu[0].header
    data = hdu[0].data
    
    # get error
    err_ind = err_filt['filt'] == filt
    err_val = err_filt[err_ind]['err'][0]
    
    filt_dict = {'name': filt, 'data': data, 'header': header, 'err': err_val}
        
    return filt_dict

F115W = load_filter('F115W')
F150W = load_filter('F150W')
F200W = load_filter('F200W')
F300M = load_filter('F300M')
F335M = load_filter('F335M')
F360M = load_filter('F360M')
F560W = load_filter('F560W')
F770W = load_filter('F770W')
F1000W = load_filter('F1000W')
F1130W = load_filter('F1130W')
F1500W = load_filter('F1500W')


filt_list = [F115W, F150W, F200W, F300M, F335M, F360M, F560W, F770W, F1000W, F1130W, F1500W]

# load F1500W data to convert the mask ra and dec
path_F1500W = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
file_F1500W = 'F1500W_reproject_to_F1500W.fits'

hdu_F1500W = fits.open(path_F1500W + file_F1500W)
data_F1500W = hdu_F1500W[0].data
head_F1500W = hdu_F1500W[0].header

F1500W_pixel_scale = head_F1500W['CDELT1']**2

w_F1500W = WCS(head_F1500W)

maxval = []

# first loop through the filters
for filt in filt_list:
    
    # prep the wcs
    w_filt = WCS(filt['header'])
    
    filt['clump_flux'] = np.zeros(len(clump_num))
    filt['clump_err'] = np.zeros(len(clump_num))
    
    # then loop through clumps
    # for i in range(len(clump_num)):
    for i in range(2):
        
        if clump_num[i] == 1:
            # load clump file
            clumppath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/clump_masks/'
            clumpname = 'F770W_k_4.33_clump{:d}_mask.txt'.format(clump_num[i])
            
            mask = ascii.read(clumppath + clumpname)
            
            filt['clump_flux'][i] = np.nansum(filt['data'][mask['col0'], mask['col1']]) * F1500W_pixel_scale * 1e6 * (np.pi/180)**2   
    
            # calculate error
            filt['clump_err'][i] = np.sqrt(len(filt['data'][mask['col0'], mask['col1']])) * filt['err'] * F1500W_pixel_scale * 1e6 * (np.pi/180)**2 
    
            print(filt['name'], filt['clump_flux'][i], filt['clump_err'][i] )
            
            # calculate statistics    
            save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/clump_stats/'
            save_path = save_path + 'clump_{:d}/'.format(clump_num[i])
            if not os.path.exists(save_path):
                os.makedirs(save_path)    
                
            save_name = 'clump_{:d}_filt_{:s}_flux_hist_and_stats.pdf'.format(clump_num[i], filt['name'])
            
            plt.figure(1)
            plt.clf()
            plt.hist(filt['data'][mask['col0'], mask['col1']], bins = 50, histtype = 'step', lw = 2, color = 'k')
            plt.xlabel('Flux (MJy/sr)')
            plt.ylabel('N Pixels')
            plt.title('Clump {:d} Filter {:s}'.format(clump_num[i], filt['name']))
            
            med = np.nanmedian(filt['data'][mask['col0'], mask['col1']])
            mean = np.nanmean(filt['data'][mask['col0'], mask['col1']])
            sumval = np.nansum(filt['data'][mask['col0'], mask['col1']])
            stdev_val = stdev(filt['data'][mask['col0'], mask['col1']])
            
            
            maxval.append(np.nanmax(filt['data'][688, 833]))

            
            plt.axvline(med, c = 'r')
            plt.axvline(mean, c = 'g')
            
            ax = plt.gca()
            s = 'Med = {:3.2f}\nMean = {:3.2f}\nSum = {:3.2f}\nstdev = {:3.2f}'.format(med, mean, sumval, stdev_val)
            plt.text(0.99, 0.99, s, ha = 'right', va = 'top', transform = ax.transAxes)
            
            plt.savefig(save_path + save_name)
            
# just a quick and dirty plot of the maxval to test my theory 
wave_arr = np.zeros(len(filt_list))
for i, filt in enumerate(filt_list):
    wave_arr[i] = get_pivot_wave(filt['name'])
    
flux_arr = np.zeros(len(filt_list))
err_arr = np.zeros(len(filt_list))
for j, filt in enumerate(filt_list):
    flux_arr[j] = filt['clump_flux'][0]
    err_arr[j] = filt['clump_err'][0]

maxval = np.array(maxval)

norm = maxval[7]/flux_arr[7]
    
plt.figure(3, figsize = (10,5))
plt.clf()
plt.scatter(wave_arr, maxval,  alpha = 0.5, label = 'Max')
plt.scatter(wave_arr, flux_arr* norm, alpha = 0.5, label = 'Sum')
# plt.errorbar(wave_arr, flux_arr*1e6, yerr = err_arr*1e6, alpha = 0.5, fmt = 'none', capsize=5, capthick=2)
plt.title('Clump 1')
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Flux (MJy/sr)')
plt.legend(loc = 'best')
            
        
      
plots = False        
if plots:        
    # now make a figure for the LW filters that we can test with
    plt.figure(1)
    plt.clf()
    plt.scatter(clumps['clump_num'], clumps['F300M_flux'], alpha = 0.5, label = 'PSF-matched')
    plt.scatter(clumps['clump_num'], F300M['clump_flux'], alpha = 0.5, label = 'regridded')
    plt.legend(loc = 'best')
    plt.xlabel('Clump Number')
    plt.ylabel('Flump Flux (Jy)')
    plt.title('F300M')
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    figname = 'SextansA_F300M_clump_PSF-matched_regrid.png'
    plt.savefig(figpath + figname, dpi = 300)
    
    F300M_perc = np.mean(abs(clumps['F300M_flux'] - F300M['clump_flux'])/(F300M['clump_flux']))
    print('F300M', F300M_perc)
    
    plt.figure(2)
    plt.clf()
    plt.scatter(clumps['clump_num'], clumps['F335M_flux'], alpha = 0.5, label = 'PSF-matched')
    plt.scatter(clumps['clump_num'], F335M['clump_flux'], alpha = 0.5, label = 'regridded')
    plt.legend(loc = 'best')
    plt.xlabel('Clump Number')
    plt.ylabel('Flump Flux (Jy)')
    plt.title('F335M')
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    figname = 'SextansA_F335M_clump_PSF-matched_regrid.png'
    plt.savefig(figpath + figname, dpi = 300)
    
    F335M_perc = np.mean(abs(clumps['F335M_flux'] - F335M['clump_flux'])/(F335M['clump_flux']))
    print('F300M', F335M_perc)
    
    
    plt.figure(3)
    plt.clf()
    plt.scatter(clumps['clump_num'], clumps['F360M_flux'], alpha = 0.5, label = 'PSF-matched')
    plt.scatter(clumps['clump_num'], F335M['clump_flux'], alpha = 0.5, label = 'regridded')
    plt.legend(loc = 'best')
    plt.xlabel('Clump Number')
    plt.ylabel('Flump Flux (Jy)')
    plt.title('F360M')
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    figname = 'SextansA_F360M_clump_PSF-matched_regrid.png'
    plt.savefig(figpath + figname, dpi = 300)
    
    F360M_perc = np.mean(abs(clumps['F360M_flux'] - F360M['clump_flux'])/(F360M['clump_flux']))
    print('F360M', F360M_perc)

    # plt.figure(4)
    # plt.clf()
    # plt.scatter(clumps['clump_num'], clumps['F1500W_flux'], alpha = 0.5, label = 'PSF-matched')
    # plt.scatter(clumps['clump_num'], F1500W['clump_flux'], alpha = 0.5, label = 'regridded')
    # plt.legend(loc = 'best')
    # plt.xlabel('Clump Number')
    # plt.ylabel('Flump Flux (Jy)')
    # plt.title('F360M')
    # figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    # figname = 'SextansA_F360M_clump_PSF-matched_regrid.png'
    # plt.savefig(figpath + figname, dpi = 300)
    
    # F360M_perc = np.mean(abs(clumps['F360M_flux'] - F360M['clump_flux'])/(F360M['clump_flux']))
    # print('F360M', F360M_perc)

make_table = False
if make_table:
        
    # save EVERYTHING 
    table = Table()
    table['clump_num'] = clump_num
    
    for filt in filt_list:
        table[filt['name'] + '_flux'] = filt['clump_flux']
        table[filt['name'] + '_err'] = filt['clump_err']
        
    wave_arr = np.zeros(len(filt_list))
    for i, filt in enumerate(filt_list):
        wave_arr[i] = get_pivot_wave(filt['name'])
        
    for i in range(len(clump_num)):
        flux_arr = np.zeros(len(filt_list))
        err_arr = np.zeros(len(filt_list))
        for j, filt in enumerate(filt_list):
            flux_arr[j] = filt['clump_flux'][i]
            err_arr[j] = filt['clump_err'][i]
        
        
        plt.figure(figsize = (10,5))
        plt.clf()
        plt.scatter(wave_arr, flux_arr*1e6, alpha = 0.5)
        plt.errorbar(wave_arr, flux_arr*1e6, yerr = err_arr*1e6, alpha = 0.5, fmt = 'none', capsize=5, capthick=2)
        plt.title('Clump {:d}'.format(clump_num[i]))
        plt.xlabel('Wavelength ($\mu$m)')
        plt.ylabel('Flux ($\mu$Jy)')
        
        savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/clump_SEDs/'
        savetext = 'clump{:d}_flux_SED_corrected.png'.format(clump_num[i])
        plt.savefig(savepath + savetext, dpi = 350)


    save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
    save_name = 'filter_flux_dendro_clump_mask_corrected.txt'
    
    ascii.write(table, save_path + save_name, overwrite = True)

