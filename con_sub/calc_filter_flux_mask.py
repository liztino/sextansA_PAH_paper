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
from astropy.wcs import WCS
from scipy.ndimage import binary_dilation, binary_fill_holes, generate_binary_structure, binary_erosion

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
clump_file = 'master_clump_cat.csv'

clumps = ascii.read(clump_path + clump_file)
clump_num = clumps['Clump']

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
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/data/aug24_final_reduction/'
    
    if filt in miri_filts:
        filepath = filepath + 'miri/'
        
        filename = glob.glob(filepath + f'*{filt}*i2d.fits')[0]
        
    elif filt in nircam_filts:
        filepath = filepath + 'nircam/'
        
        if filt in lw_filts:
            filename = glob.glob(filepath + '*{:s}*jhat_i2d.fits'.format(filt.lower()))[0]
        else:
            filename = glob.glob(filepath + '*{:s}*i2d.fits'.format(filt.lower()))[0]
    
    # load the filter and info
    hdu = fits.open(filename)
    header = hdu['SCI'].header
    data = hdu['SCI'].data
    err = hdu['ERR'].data
    
    filt_dict = {'name': filt, 'data': data, 'header': header, 'err': err}
        
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

filt_list = [F200W]
filt_list = [F115W, F150W, F200W, F300M, F335M, F360M, F560W, F770W, F1000W, F1130W, F1500W]


# load F1500W data to convert the mask ra and dec
path_F1500W = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
file_F1500W = 'F1500W_reproject_to_F1500W.fits'

hdu_F1500W = fits.open(path_F1500W + file_F1500W)
data_F1500W = hdu_F1500W[0].data
head_F1500W = hdu_F1500W[0].header

F1500W_pixel_scale = head_F1500W['CDELT1']**2

w_F1500W = WCS(head_F1500W)

# first loop through the filters
for filt in filt_list:
    
    # prep the wcs
    w_filt = WCS(filt['header'])
    
    filt['clump_flux'] = np.zeros(len(clump_num))
    filt['clump_err'] = np.zeros(len(clump_num))
    
    # then loop through clumps
    for i in range(len(clump_num)):
        
        # if clump_num[i] == 1:
            # load clump file
        clumppath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/clump_masks/'
        clumpname = 'F770W_k_4.33_clump{:d}_mask.txt'.format(clump_num[i])
        
        mask = ascii.read(clumppath + clumpname)
                    
        # convert the mask pixel coordinates to world
        ra, dec = w_F1500W.pixel_to_world_values(mask['col1'], mask['col0'])
        
        # convert the ra_dec to pixel coordinates for the given filter
        x_filt, y_filt = w_filt.world_to_pixel_values(ra, dec)
        
        # convert to integer
        x_filt = x_filt.astype(int)
        y_filt = y_filt.astype(int)
        
        # create a mask
        filt_mask = np.zeros_like(filt['data'], dtype = bool)
        filt_mask[y_filt, x_filt]= True
        
        # now fill the holes
        # insane code that chatgpt gave me
        struct_element = generate_binary_structure(2, 2)
        expanded_mask = binary_dilation(filt_mask, structure=struct_element, iterations=5)
        filled_mask = binary_fill_holes(expanded_mask)
        filled_mask = binary_erosion(filled_mask,  structure=struct_element, iterations=5)
        
        
        # sum the filter 
        pixel_area = filt['header']['CDELT1']**2
        flux_sum = np.nansum(filt['data'][filled_mask]) * pixel_area * 1e6 * (np.pi/180)**2   
        
        # add to the dictionary 
        filt['clump_flux'][i] = flux_sum
        
        # now calculate the uncertainty 
        err = np.sqrt(np.nansum(np.square(filt['err'][filled_mask]))) * pixel_area * 1e6 * (np.pi/180)**2   
        filt['clump_err'][i] = err
        
        print(filt['name'], flux_sum, err)
        
        # calculate sizes to compare to the excess
        len_orig_mask = len(mask)
        len_new_mask = len((filt['data'][filled_mask]))
        
        # in degrees
        area_orig_mask =  F1500W_pixel_scale * len_orig_mask
        area_new_mask =  pixel_area * len_new_mask
        
        filt['clump_area'][i] = area_new_mask * (180/np.pi)**2   
        
        print('Area orig', area_orig_mask)
        print('Area new', area_new_mask)
        
        diff = area_new_mask - area_orig_mask
        
        fiducial_flux = 7.2061648e-06/(diff * 1e6 * (np.pi/180)**2)
        
        # print('diff', diff)
            
        #     # fiducial_flux = 0.06
            
        #     # flux_add = fiducial_flux * diff * 1e6 * (np.pi/180)**2 
            
        #     print('additional flux', fiducial_flux)
    
            
        #     # code to plot and test
        #     plot = False
        #     if plot:
        #         plt.figure(1)
        #         plt.clf()
        #         plt.imshow(filt_mask, vmin = 0, vmax = 0.1, origin = 'lower')
        #         # plt.xlim(694, 725)
        #         # plt.ylim(495, 530)
                
        #         plt.figure(2)
        #         plt.clf()
        #         plt.imshow(filt['data'], vmin = 0, vmax = 0.1, origin = 'lower')
        #         # plt.xlim(694, 725)
        #         # plt.ylim(495, 530)
                
        #         plt.figure(3)
        #         plt.clf()
        #         plt.imshow(filled_mask, vmin = 0, vmax = 0.1, origin = 'lower')
        #         # plt.xlim(694, 725)
        #         # plt.ylim(495, 530)
                
        #         masked_data = np.copy(filt['data'])
        #         masked_data[~filled_mask] = np.nan
        #         masked_data[filled_mask] = masked_data[filled_mask]
                
        #         plt.figure(4)
        #         plt.clf()
        #         plt.imshow(masked_data, vmin = 0, vmax = 0.1, origin = 'lower')
        #         # plt.xlim(694, 725)
        #         # plt.ylim(495, 530)
    


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
        savetext = 'clump{:d}_flux_SED.png'.format(clump_num[i])
        # plt.savefig(savepath + savetext, dpi = 350)


    save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
    save_name = 'filter_flux_dendro_clump_mask.txt'
    
    # ascii.write(table, save_path + save_name, overwrite = True)

