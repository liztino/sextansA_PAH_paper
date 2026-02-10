#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculates PAH flux in each filter using a mask from the dendrograms of another

TODO: add the filter fluxes (or maybe do that in a separate script?)

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.table import Table

from get_pivot_wave import get_pivot_wave

# load pah flux files

def load_PAH(filt_mid, k):
    datapath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
    dataname = 'SextansA_{:s}_k_method_pah_k_{:3.2f}.fits'.format(filt_mid, k)
    pah_hdr = fits.open(datapath + dataname)[0]
    pah_data = pah_hdr.data
    pah_head = pah_hdr.header
    
    errname = 'SextansA_{:s}_k_method_err_k_{:3.2f}.fits'.format(filt_mid, k)
    err_hdr = fits.open(datapath + errname)[0]
    err_data = err_hdr.data
    err_head = err_hdr.header
    
    conname = 'SextansA_{:s}_k_method_con_k_{:3.2f}.fits'.format(filt_mid, k)
    con_hdr = fits.open(datapath + conname)[0]
    con_data = con_hdr.data
    con_head = con_hdr.header    
    
    return pah_data, pah_head, err_data, con_data

# load master clump file
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_file = 'master_clump_cat.csv'

clumps = ascii.read(clump_path + clump_file)
clump_num = clumps['Clump']

# constants 
pixel_area = (3.0555555555555E-05)**2

# load filters
def load_filter(filt):
    
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
    
    # load the middle filter we will be continuum subtracting from
    filename = filepath + f'{filt}_reproject_to_F1500W_rot'
    hdu = fits.open(filename + '.fits')
    header = hdu[0].header
    data = hdu[0].data
    pivot = get_pivot_wave(filt)
    
    # add a field for clump flux
    clump_flux = np.zeros(len(clump_num))

    filt_dict = {'name': filt, 'data': data, 'header': header, 'wave': pivot, 'clump_flux': clump_flux}
        
    return filt_dict

F300M = load_filter('F300M')
F335M = load_filter('F335M')
F360M = load_filter('F360M')
F560W = load_filter('F560W')
F770W = load_filter('F770W')
F1000W = load_filter('F1000W')
F1130W = load_filter('F1130W')
F1500W = load_filter('F1500W')

filt_list = [F300M, F335M, F360M, F560W, F770W, F1000W, F1130W, F1500W]

# # quickly get a value for the filters at each pixel
# for filt in filt_list:
#     px = 833
#     py = 688
    
#     print('Filter {:s} = {:3.2f}'.format(filt['name'], filt['data'][py, px]))


# # load F1500W data
# path_F1500W = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject/'
# file_F1500W = 'F1500W_reproject_to_F1500W.fits'

# hdu_F1500W = fits.open(path_F1500W + file_F1500W)
# data_F1500W = hdu_F1500W[0].data
# clump_flux_F1500W = np.zeros(len(clump_num))

# uncertainties
pah3_sigma = 0.0036
pah7_sigma =  0.0066
pah11_sigma = 0.0138

#######################################
####### first set of k values #######
######################################
pah_3_k1, head_3_k1, err_3_k1, con_3_k1 = load_PAH('F335M', 2.07)
pah_7_k1, head_7_k1, err_7_k1, con_7_k1 = load_PAH('F770W', 4.33)
pah_11_k1, head_11_k1, err_11_k1, con_11_k1 = load_PAH('F1130W', 7.21)

clump_flux_3_k1 = np.zeros(len(clump_num))
clump_flux_7_k1 = np.zeros(len(clump_num))
clump_flux_11_k1 = np.zeros(len(clump_num))


clump_con_3_k1 = np.zeros(len(clump_num))
clump_con_7_k1 = np.zeros(len(clump_num))
clump_con_11_k1 = np.zeros(len(clump_num))

# error calculated thrrough just the sigma from the data
clump_flux_3_k1_err_sig = np.zeros(len(clump_num))
clump_flux_7_k1_err_sig = np.zeros(len(clump_num))
clump_flux_11_k1_err_sig = np.zeros(len(clump_num))

# error calculated propagating the uncertainty from k 
clump_flux_3_k1_err_k = np.zeros(len(clump_num))
clump_flux_7_k1_err_k = np.zeros(len(clump_num))
clump_flux_11_k1_err_k = np.zeros(len(clump_num))

#######################################
####### second set of k values #######
######################################
pah_3_k2, head_3_k2, err_3_k2, con_3_k2 = load_PAH('F335M', 4.45)
pah_7_k2, head_7_k2, err_7_k2, con_7_k2 = load_PAH('F770W', 5.84)
pah_11_k2, head_11_k2, err_11_k2, con_11_k2 = load_PAH('F1130W', 10.17)

clump_flux_3_k2 = np.zeros(len(clump_num))
clump_flux_7_k2 = np.zeros(len(clump_num))
clump_flux_11_k2 = np.zeros(len(clump_num))

# error calculated thrrough just the sigma from the data
clump_flux_3_k2_err_sig = np.zeros(len(clump_num))
clump_flux_7_k2_err_sig = np.zeros(len(clump_num))
clump_flux_11_k2_err_sig = np.zeros(len(clump_num))

# error calculated propagating the uncertainty from k 
clump_flux_3_k2_err_k = np.zeros(len(clump_num))
clump_flux_7_k2_err_k = np.zeros(len(clump_num))
clump_flux_11_k2_err_k = np.zeros(len(clump_num))


# loop through clump masks
for i in range(len(clump_num)):
    save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/clump_masks/'
    save_name = 'F770W_k_4.33_clump{:d}_mask.txt'.format(clump_num[i])
    
    mask = ascii.read(save_path + save_name)
    
    # if clump_num[i] == 16:
    #     print(mask)
    #######################################
    ####### first set of k values #######
    ######################################    
    clump_flux_3_k1[i] = np.nansum(pah_3_k1[mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_7_k1[i] = np.nansum(pah_7_k1[mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2
    clump_flux_11_k1[i] = np.nansum(pah_11_k1[mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2
    
    ######################################    
    clump_con_3_k1[i] = np.nansum(con_3_k1[mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2   
    clump_con_7_k1[i] = np.nansum(con_7_k1[mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2
    clump_con_11_k1[i] = np.nansum(con_11_k1[mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2
    
    
    # error prop will be flux_unc * sqrt(Npix) * constants
    # error calculated thrrough just the sigma from the data
    clump_flux_3_k1_err_sig[i] = np.sqrt(len(pah_3_k1[mask['col0'], mask['col1']])) * pah3_sigma * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_7_k1_err_sig[i] = np.sqrt(len(pah_7_k1[mask['col0'], mask['col1']])) * pah7_sigma * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_11_k1_err_sig[i] = np.sqrt(len(pah_11_k1[mask['col0'], mask['col1']])) * pah11_sigma * pixel_area * 1e6 * (np.pi/180)**2   

    # error calculated propagating the uncertainty from k 
    clump_flux_3_k1_err_k[i] = np.sqrt(np.nansum((err_3_k1[mask['col0'], mask['col1']])**2)) * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_7_k1_err_k[i] = np.sqrt(np.nansum((err_7_k1[mask['col0'], mask['col1']])**2)) * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_11_k1_err_k[i] = np.sqrt(np.nansum((err_11_k1[mask['col0'], mask['col1']])**2)) * pixel_area * 1e6 * (np.pi/180)**2   


    #######################################
    ####### second set of k values #######
    ######################################
    clump_flux_3_k2[i] = np.nansum(pah_3_k2[mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_7_k2[i] = np.nansum(pah_7_k2[mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2
    clump_flux_11_k2[i] = np.nansum(pah_11_k2[mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2
    
    # error prop will be flux_unc * sqrJt(Npix) * constants
    clump_flux_3_k2_err_sig[i] = np.sqrt(len(pah_3_k2[mask['col0'], mask['col1']])) * pah3_sigma * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_7_k2_err_sig[i] = np.sqrt(len(pah_7_k2[mask['col0'], mask['col1']])) * pah7_sigma * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_11_k2_err_sig[i] = np.sqrt(len(pah_11_k2[mask['col0'], mask['col1']])) * pah11_sigma * pixel_area * 1e6 * (np.pi/180)**2   

    # error calculated propagating the uncertainty from k 
    clump_flux_3_k2_err_k[i] = np.sqrt(np.nansum((err_3_k2[mask['col0'], mask['col1']])**2)) * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_7_k2_err_k[i] = np.sqrt(np.nansum((err_7_k2[mask['col0'], mask['col1']])**2)) * pixel_area * 1e6 * (np.pi/180)**2   
    clump_flux_11_k2_err_k[i] = np.sqrt(np.nansum((err_11_k1[mask['col0'], mask['col1']])**2)) * pixel_area * 1e6 * (np.pi/180)**2   


    # get flux for filters
    for filt in filt_list:
        filt['clump_flux'][i] = np.nansum(filt['data'][mask['col0'], mask['col1']]) * pixel_area * 1e6 * (np.pi/180)**2

# quick test to ensure math works
print('PAH', clump_flux_7_k1[0])
print('CON', clump_con_7_k1[0])
print('PAH + CON', clump_flux_7_k1[0] + clump_con_7_k1[0])
print('F770W', F770W['clump_flux'][0])



plot = False   

if plot: 
    plt.figure(2, figsize = (10,10))
    plt.clf()
    fig, axs = plt.subplots(nrows=3, ncols=1, num = 2,  sharex=True)
    N = np.arange(1, len(clump_flux_3_k1) + 1)
    
    axs[0].scatter(N, clump_flux_3_k1_err_sig, marker = 'o', alpha = 0.5, label = 'Sigma')
    axs[0].scatter(N, clump_flux_3_k1_err_k, marker = 'o', alpha = 0.5, label = 'Sigma+k')
    axs[0].legend(loc = 'best')
    axs[0].set_yscale('log')
    axs[0].set_title('3.3 Feature k = 2.07')
    
    axs[1].scatter(N, clump_flux_7_k1_err_sig, marker = 'o', alpha = 0.5, label = 'Sigma')
    axs[1].scatter(N, clump_flux_7_k1_err_k, marker = 'o', alpha = 0.5, label = 'Sigma+k')
    axs[1].set_yscale('log')
    axs[1].set_title('7.7 Feature k = 4.33')
    
    axs[2].scatter(N, clump_flux_11_k1_err_sig, marker = 'o', alpha = 0.5, label = 'Sigma')
    axs[2].scatter(N, clump_flux_11_k1_err_k, marker = 'o', alpha = 0.5, label = 'Sigma+k')
    axs[2].set_yscale('log')
    axs[2].set_title('11.3 Feature k = 7.21')
    
    plt.xlabel('Clump Number')
    plt.ylabel('Flux Error')
    plt.suptitle('Comparison of the Flux Error (First Set of k values)')
    
    # plt.subplots_adjust(wspace=0.5, hspace=0.01)
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    figname = 'SextansA_clump_flux_error_comparison_k1_new.png'
    plt.savefig(figpath + figname, dpi = 300)
    
    
    # plt.figure(3)
    # plt.clf()
    # N = np.arange(1, len(clump_flux_3_k1) + 1)
    # plt.errorbar(N, clump_flux_3_k1, yerr = clump_flux_3_k1_err_k, fmt = '.', alpha = 0.5)
    # plt.errorbar(N, clump_flux_3_k2, yerr = clump_flux_3_k2_err_k, fmt = '.', alpha = 0.5)
    
    plt.figure(3, figsize = (10,10))
    plt.clf()
    fig, axs = plt.subplots(nrows=3, ncols=1, num = 3,  sharex=True)
    N = np.arange(1, len(clump_flux_3_k1) + 1)
    
    axs[0].errorbar(N, clump_flux_3_k1, yerr = clump_flux_3_k1_err_k, fmt = '.', alpha = 0.5, label = 'k = 2.07')
    axs[0].errorbar(N, clump_flux_3_k2, yerr = clump_flux_3_k2_err_k, fmt = '.', alpha = 0.5, label = 'k = 4.45')
    axs[0].legend(loc = 'best')
    axs[0].set_yscale('log')
    axs[0].set_title('3.3 Feature')
    
    axs[1].errorbar(N, clump_flux_7_k1, yerr = clump_flux_7_k1_err_k, fmt = '.', alpha = 0.5, label = 'k = 4.33')
    axs[1].errorbar(N, clump_flux_7_k2, yerr = clump_flux_7_k2_err_k, fmt = '.', alpha = 0.5, label = 'k = 5.84')
    axs[1].legend(loc = 'best')
    axs[1].set_yscale('log')
    axs[1].set_title('7.7 Feature')
    
    axs[2].errorbar(N, clump_flux_11_k1, yerr = clump_flux_11_k1_err_k, fmt = '.', alpha = 0.5, label = 'k = 7.21')
    axs[2].errorbar(N, clump_flux_11_k2, yerr = clump_flux_11_k2_err_k, fmt = '.', alpha = 0.5, label = 'k = 10.17')
    axs[2].legend(loc = 'best')
    axs[2].set_yscale('log')
    axs[2].set_title('11.3 Feature')
    
    plt.xlabel('Clump Number')
    plt.ylabel('Flux')
    plt.suptitle('Comparison of the Flux')
    
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    figname = 'SextansA_clump_flux_errorbar_new.png'
    plt.savefig(figpath + figname, dpi = 300)
    
    # calculate values compared to the total PAH
    tot_PAH = clump_flux_3_k1 + clump_flux_7_k1 + clump_flux_11_k1
    sigma_PAH_3 = np.sum(clump_flux_3_k1)
    sigma_PAH_7 = np.sum(clump_flux_7_k1)
    sigma_PAH_11 = np.sum(clump_flux_11_k1)
    
    
    print(np.nansum(tot_PAH)/np.nansum(F1500W['clump_flux']))
    
    plt.figure(4)
    plt.clf()
    plt.scatter(N, tot_PAH/F1500W['clump_flux'])
    plt.xlabel('Clump Number')
    plt.ylabel('$\Sigma$PAH/F1500W', size = 'x-large')
    plt.axhline(0.36, c = 'k')
    
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    figname = 'SextansA_clump_Sigma_PAH_F1500W.png'
    plt.savefig(figpath + figname, dpi = 300)
    
    # comparing each feature to F1500W 
    
    plt.figure(5, figsize = (10,10))
    plt.clf()
    fig, axs = plt.subplots(nrows=3, ncols=1, num = 5,  sharex=True)
    N = np.arange(1, len(clump_flux_3_k1) + 1)
    
    axs[0].scatter(N, clump_flux_3_k1/F1500W['clump_flux'], marker = 'o', alpha = 0.5, label = 'k = 2.07')
    axs[0].set_title('3.3 Feature')
    axs[0].set_ylabel('3.3/F1500W')
    axs[0].axhline(1, c = 'k', lw = 0.5)
    
    axs[1].scatter(N, clump_flux_7_k1/F1500W['clump_flux'],  marker = 'o', alpha = 0.5, label = 'k = 4.33')
    axs[1].set_title('7.7 Feature')
    axs[1].set_ylabel('7.7/F1500W')
    axs[1].axhline(1, c = 'k', lw = 0.5)
    
    axs[2].scatter(N, clump_flux_11_k1/F1500W['clump_flux'],  marker = 'o', alpha = 0.5, label = 'k = 7.21')
    axs[2].set_title('11.3 Feature')
    axs[2].set_ylabel('11.3/F1500W')
    axs[2].axhline(1, c = 'k', lw = 0.5)
    
    plt.xlabel('Clump Number')
    plt.suptitle('PAH Flux Compared to F1500W')
    
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    figname = 'SextansA_clump_PAH_to_F1500W_new.png'
    plt.savefig(figpath + figname, dpi = 300)
    
    
    plt.figure(6, figsize = (10,10))
    plt.clf()
    fig, axs = plt.subplots(nrows=3, ncols=1, num = 6,  sharex=True)
    N = np.arange(1, len(clump_flux_3_k1) + 1)
    
    axs[0].scatter(N, clump_flux_3_k1/clump_con_3_k1, marker = 'o', alpha = 0.5, label = 'k = 2.07')
    axs[0].set_title('3.3 Feature')
    axs[0].set_ylabel('3.3 PAH/3.3 Con')
    axs[0].axhline(1, c = 'k', lw = 0.5)
    
    axs[1].scatter(N, clump_flux_7_k1/clump_con_7_k1,  marker = 'o', alpha = 0.5, label = 'k = 4.33')
    axs[1].set_title('7.7 Feature')
    axs[1].set_ylabel('7.7PAH/7.7 Con')
    axs[1].axhline(1, c = 'k', lw = 0.5)
    
    axs[2].scatter(N, clump_flux_11_k1/clump_con_11_k1,  marker = 'o', alpha = 0.5, label = 'k = 7.21')
    axs[2].set_title('11.3 Feature')
    axs[2].set_ylabel('11.3 PAH/11.3 Con')
    axs[2].axhline(1, c = 'k', lw = 0.5)
    
    plt.xlabel('Clump Number')
    plt.suptitle('PAH Flux Compared to F1500W')
    
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    figname = 'SextansA_clump_PAH_to_continuum_new.png'
    plt.savefig(figpath + figname, dpi = 300)

    # create an SED of each clump
    
    # axs[1].scatter(N, clump_flux_7_k1_err_sig, marker = 'o', alpha = 0.5, label = 'Sigma')
    # axs[1].scatter(N, clump_flux_7_k1_err_k, marker = 'o', alpha = 0.5, label = 'Sigma+k')
    # axs[1].set_yscale('log')
    # axs[1].set_title('7.7 Feature k = 4.33')
    
    # axs[2].scatter(N, clump_flux_11_k1_err_sig, marker = 'o', alpha = 0.5, label = 'Sigma')
    # axs[2].scatter(N, clump_flux_11_k1_err_k, marker = 'o', alpha = 0.5, label = 'Sigma+k')
    # axs[2].set_yscale('log')
    # axs[2].set_title('11.3 Feature k = 7.21')
    
    # plt.xlabel('Clump Number')
    # plt.ylabel('Flux Error')
    # plt.suptitle('Comparison of the Flux Error (First Set of k values)')
    
    # # plt.subplots_adjust(wspace=0.5, hspace=0.01)
    # figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/compare_clump_flux/'
    # figname = 'SextansA_clump_flux_error_comparison_k1.png'
    # plt.savefig(figpath + figname, dpi = 300)
    
    # plt.scatter(N, (clump_flux_3_k1 - clump_flux_3_k2)/clump_flux_3_k2, marker = 'o', alpha = 0.5)

    
# # make final table that will be wrapped into clump properties table
# make_table = True

# if make_table:
        
#     # save EVERYTHING then only use parts
#     arrays = [clump_num, clump_flux_3_k1, clump_flux_3_k1_err_k,  clump_flux_3_k2, clump_flux_3_k2_err_k, clump_con_3_k1,
#               clump_flux_7_k1, clump_flux_7_k1_err_k,  clump_flux_7_k2, clump_flux_7_k2_err_k, clump_con_7_k1,
#               clump_flux_11_k1, clump_flux_11_k1_err_k,  clump_flux_11_k2, clump_flux_11_k2_err_k, clump_con_11_k1]
        
#     names = ['clump_num', 'clump_flux_3_k1', 'clump_flux_3_k1_err_k', 'clump_flux_3_k2', 'clump_flux_3_k2_err_k', 'clump_con_3_k1', 
#              'clump_flux_7_k1', 'clump_flux_7_k1_err_k', 'clump_flux_7_k2', 'clump_flux_7_k2_err_k', 'clump_con_7_k1', 
#              'clump_flux_11_k1', 'clump_flux_11_k1_err_k', 'clump_flux_11_k2', 'clump_flux_11_k2_err_k', 'clump_con_11_k1']
    
#     for filt in filt_list:
#         arrays.append(filt['clump_flux'])
#         names.append(filt['name'] + '_flux')
    
#     table = Table(arrays, names = names)
#     save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
#     save_name = 'dendro_mask_clump_flux_with_filt.txt'
    
#     ascii.write(table, save_path + save_name, overwrite = True)

