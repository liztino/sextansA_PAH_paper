#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Performs continuum subtraction with the new method that accounts for PAH contamination in the continuum filters

@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
# import reproject
# from astropy.wcs import WCS
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.ion()

# custom functions 
from get_pivot_wave import get_pivot_wave
import k_eq
import k_eq_with_err

def k_consub(filt_mid, filt_low, filt_up, k, k_err, save = True):

    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
    
    # load the middle filter we will be continuum subtracting from
    mid_name = filepath + f'{filt_mid}_reproject_to_F1500W_rot'
    hdu_mid = fits.open(mid_name + '.fits')
    header_mid = hdu_mid[0].header
    data_mid = hdu_mid[0].data
    pivot_mid = get_pivot_wave(filt_mid)

    # load the upper filter     that we will use to match to the other filters
    up_name = filepath + f'{filt_up}_reproject_to_F1500W_rot'
    hdu_up = fits.open(up_name  + '.fits')
    header_up = hdu_up[0].header
    data_up = hdu_up[0].data
    pivot_up = get_pivot_wave(filt_up)

    # load the lowest wavelength filter
    low_name = filepath + f'{filt_low}_reproject_to_F1500W_rot'
    hdu_low = fits.open(low_name + '.fits')
    header_low = hdu_low[0].header
    data_low = hdu_low[0].data
    pivot_low = get_pivot_wave(filt_low)
    
    # # really dumb math check to ensure k method is working 
    # # test low first
    # data_low = 5
    # data_mid = 10
    # data_up = 8
    # f1_err = 0.1
    # f2_err = 0.1
    # f3_err = 0.1
    # k_err = 0.4
    # consub = k_eq_with_err.get_pah_low(data_low, data_mid, data_up, pivot_low, pivot_mid, pivot_up, k, f1_err, f2_err, f3_err, k_err)
    # print('PAH', consub['pah'])
    # print('CON', consub['con'])
    # print('PAH + CON', consub['pah'] + consub['con'])
    
    # # test up 
    # data_low = 5
    # data_mid = 10
    # data_up = 8
    # f1_err = 0.1
    # f2_err = 0.1
    # f3_err = 0.1
    # k_err = 0.4
    # consub = k_eq_with_err.get_pah_up(data_low, data_mid, data_up, pivot_low, pivot_mid, pivot_up, k, f1_err, f2_err, f3_err, k_err)
    # print('PAH', consub['pah'])
    # print('CON', consub['con'])
    # print('PAH + CON', consub['pah'] + consub['con'])
    
    
    # print(pivot_mid, pivot_up, pivot_low)
    if filt_mid == 'F770W':
        f1_err = 0.004516
        f2_err = 0.004589
        f3_err = 0.006651
        # k_err = 0.73
        pah_unc = 0.0066
        
        consub = k_eq_with_err.get_pah_low(data_low, data_mid, data_up, pivot_low, pivot_mid, pivot_up, k, f1_err, f2_err, f3_err, k_err)
    else: 
        # 3.3
        if filt_mid == 'F335M':
            f1_err = 0.003054
            f2_err = 0.002090
            f3_err = 0.003054
            pah_unc = 0.0035
        
        # 11.3
        elif filt_mid == 'F1130W':
            f1_err = 0.006651
            f2_err = 0.011505
            f3_err = 0.020477
            pah_unc = 0.0138
        
            # k_err = 1.24
        
        consub = k_eq_with_err.get_pah_up(data_low, data_mid, data_up, pivot_low, pivot_mid, pivot_up, k, f1_err, f2_err, f3_err, k_err)
    
    if save:
        save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
        save_name = 'SextansA_{:s}_k_method_pah_k_{:3.2f}.fits'.format(filt_mid, k)
        fits.writeto(save_path + save_name, consub['pah'], header = header_mid, overwrite = True)
        
        save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
        save_name = 'SextansA_{:s}_k_method_err_k_{:3.2f}.fits'.format(filt_mid, k)
        fits.writeto(save_path + save_name, consub['pah_err'], header = header_mid, overwrite = True)

        save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
        save_name = 'SextansA_{:s}_k_method_con_k_{:3.2f}.fits'.format(filt_mid, k)
        fits.writeto(save_path + save_name, consub['con'], header = header_mid, overwrite = True)
        
    save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
    save_name = 'SextansA_{:s}_k_method_slope_k_{:3.2f}.fits'.format(filt_mid, k)
    fits.writeto(save_path + save_name, consub['slope'], header = header_mid, overwrite = True)
        
    print(consub, np.nanmax(consub['pah']))
    
    plt.figure(3,  figsize = (30, 5))
    plt.clf()
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows = 1, ncols = 4, num = 3, figsize = (20, 5))
    
    im = ax1.imshow(consub['pah'], origin = 'lower', vmin = 0, vmax = 0.05, cmap = 'magma')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar(im, cax=cax)
    ax1.set_xlim(820, 850)
    ax1.set_ylim(680, 710)
    ax1.set_title('PAH Flux')
    ax1.tick_params(axis='both', labelleft=False, labelbottom = False)

    avg_unc = np.nanmean(consub['pah_err'][680:710, 820:850])
    im = ax2.imshow(consub['pah_err'], origin = 'lower', vmin = 0, vmax = 0.05, cmap = 'magma')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar(im, cax=cax)
    ax2.set_xlim(820, 850)
    ax2.set_ylim(680, 710)
    ax2.set_title('PAH Uncertainty, avg = {:5.4f}'.format(avg_unc))
    ax2.tick_params(axis='both', labelleft=False, labelbottom = False)

    
    im = ax3.imshow(consub['pah']/consub['pah_err'], origin = 'lower', vmin = 0, vmax = 10, cmap = 'rainbow')
    ax3.contour(consub['pah']/consub['pah_err'], levels = [3,5, 10], colors= 'gray', linewidths = 2)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar(im, cax=cax)    
    ax3.set_xlim(820, 850)
    ax3.set_ylim(680, 710)
    ax3.set_title('SNR Map')
    ax3.tick_params(axis='both', labelleft=False, labelbottom = False)
    
    im = ax4.imshow(consub['pah']/pah_unc, origin = 'lower', vmin = 0, vmax = 10, cmap = 'rainbow')
    ax4.contour(consub['pah']/pah_unc, levels =[3, 5, 10], colors = 'gray', linewidths = 2)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = fig.colorbar(im, cax=cax)    
    ax4.set_xlim(820, 850)
    ax4.set_ylim(680, 710)
    ax4.set_title('SNR Map from single SNR')
    ax4.tick_params(axis='both', labelleft=False, labelbottom = False)
    
    fig.suptitle('{:s} PAH: k = {:3.2f} Unc =  {:5.4f}'.format(filt_mid, k, pah_unc))
    
    savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/'
    savename = '{:s}_k_method_k{:3.2f}_uncertainty.pdf'.format(filt_mid, k)
    # plt.savefig(savepath + savename)




# k_consub('F335M', 'F300M', 'F360M', k=4.45, k_err = 0.39, save = True)
# k_consub('F770W', 'F560W', 'F1000W', k=5.84, k_err = 0.73, save = True)
k_consub('F1130W', 'F1000W', 'F1500W', k=7.21, k_err = 0.92, save = True)




# # for IC 1613
# def k_consub(filt_mid, filt_low, filt_up, k, save = True):

#     filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613/aug23_reduction/reproject_north/'
    
#     # load the middle filter we will be continuum subtracting from
#     mid_name = filepath + f'IC1613_jw2391_{filt_mid}_reproject_mcf_north'
#     hdu_mid = fits.open(mid_name + '.fits')
#     header_mid = hdu_mid[0].header
#     data_mid = hdu_mid[0].data
#     pivot_mid = get_pivot_wave(filt_mid)

#     # load the upper filter     that we will use to match to the other filters
#     up_name = filepath + f'IC1613_jw2391_{filt_up}_reproject_mcf_north'
#     hdu_up = fits.open(up_name  + '.fits')
#     header_up = hdu_up[0].header
#     data_up = hdu_up[0].data
#     pivot_up = get_pivot_wave(filt_up)

#     # load the lowest wavelength filter
#     low_name = filepath + f'IC1613_jw2391_{filt_low}_reproject_mcf_north'
#     hdu_low = fits.open(low_name + '.fits')
#     header_low = hdu_low[0].header
#     data_low = hdu_low[0].data
#     pivot_low = get_pivot_wave(filt_low)
    
#     # print(pivot_mid, pivot_up, pivot_low)
#     if filt_mid == 'F770W':
#         consub = k_eq.get_pah_low(data_low, data_mid, data_up, pivot_low, pivot_mid, pivot_up, k)
#     else: 
#         consub = k_eq.get_pah_up(data_low, data_mid, data_up, pivot_low, pivot_mid, pivot_up, k)
    
#     if save:
#         save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613/aug23_reduction/quick_k_consub/'
#         save_name = f'IC1613_{filt_mid}_k_method_pah.fits'
#         fits.writeto(save_path + save_name, consub['pah'], header = header_mid, overwrite = True)

#         save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/IC1613/aug23_reduction/quick_k_consub/'
#         save_name = f'IC1613_{filt_mid}_k_method_con.fits'
#         fits.writeto(save_path + save_name, consub['con'], header = header_mid, overwrite = True)


# k_consub('F335M', 'F300M', 'F360M', k=1.6, save = True)
# k_consub('F770W', 'F560W', 'F1000W', k=5.87, save = True)
# k_consub('F1130W', 'F1000W', 'F1500W', k=10.17, save = True)


# # pah = (F1000W_pivot - F560W_pivot)/ (F1000W_pivot - F560W_pivot - k*F1000W_pivot - k*F770W_pivot) * (
# #     (F770W_data - (((F1000W_data * (F770W_pivot - F560W_pivot)))/(F1000W_pivot - F560W_pivot)) + (((F560W_data * (F770W_pivot - F560W_pivot)))/(F1000W_pivot - F560W_pivot))))


# def fp2_v2(f1, f2, f3, lam1, lam2, lam3, k):
#     # equation calculated through Julia's math
#     fp1 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(1 - k - ((lam2 - lam1)/(lam3 - lam1)))
#     fp2 = fp1*k
    
#     print('fp1', fp1)
#     print('fp2', fp2)
    
#     return fp2

# def fc2_v2(f1, f2, f3, lam1, lam2, lam3, k):
#     # equation calculated through Julia's math
#     fp1 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(1 - k - ((lam2 - lam1)/(lam3 - lam1)))
#     fc1 = f1 - fp1    
#     fc2 = (((f3 - fc1)/(lam3 - lam1)) * (lam2 - lam1)) + fc1
    
#     print('fc1', fc1)
#     print('fc2', fc2)
#     print('slope', (((f3 - fc1)/(lam3 - lam1))))
#     # print('delta_lam', lam2 - lam1)
    
#     return fc2

# k = 5.505

# pah = fp2_v2(F560W_data, F770W_data, F1000W_data, F560W_pivot, F770W_pivot, F1000W_pivot, k)
# con = fc2_v2(F560W_data, F770W_data, F1000W_data, F560W_pivot, F770W_pivot, F1000W_pivot, k)


# plt.figure(2)
# plt.clf()
# plt.imshow(pah, origin = 'lower', vmin = 0, vmax = 0.01)


# save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images_rot/'
# save_name = 'SextansA_F770W_pah_newmethod.fits'
# fits.writeto(save_path + save_name, pah, header = F770W_head, overwrite = True)

# save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images_rot/'
# save_name = 'SextansA_F770W_con_newmethod.fits'
# fits.writeto(save_path + save_name, con, header = F770W_head, overwrite = True)

# # open up the original subtraction, the interpolation method
# inpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images_rot/'
# inname = 'SextansA_F770W_sub_match_rot_pah.fits'
# data = fits.open(inpath + inname)[0].data

# ratio = pah/data

# save_name = 'SextansA_F770W_pah_ratio_methods.fits'
# fits.writeto(save_path + save_name, ratio, header = F770W_head, overwrite = True)


