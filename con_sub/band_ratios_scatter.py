#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates band ratio plots and images from the data


@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import reproject
from astropy.wcs import WCS
import numpy.ma as ma
import astropy.visualization.wcsaxes as vis
import astropy.units as u 
import matplotlib
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm
import seaborn as sns
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.ion()

# custom functions 
from get_pivot_wave import get_pivot_wave
import k_eq

def load_filter(filt):
    
    filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject_rot/'
    
    # load the middle filter we will be continuum subtracting from
    filename = filepath + f'{filt}_reproject_to_F1500W_rot'
    hdu = fits.open(filename + '.fits')
    header = hdu[0].header
    data = hdu[0].data
    pivot = get_pivot_wave(filt)
    
    filt_dict = {'name': filt, 'data': data, 'header': header, 'wave': pivot}
    
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


# load the continuum subtracted PAH map to make a mask
# using the F1130W to include the most flux
filt = 'F770W'
cutoff = 0.003
consub_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/'
savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/maps/'

k_method_file = consub_path + f'k_method/SextansA_{filt}_k_method_pah.fits'
k_method_hdu = fits.open(k_method_file)
k_method_data = k_method_hdu[0].data
k_method_head = k_method_hdu[0].header

# define the pixel region to use 
# see the box defined in the region file 
x1 = 793
x2 = 948
y1 = 651
y2 = 798

mask1 = np.zeros_like(k_method_data, dtype = bool)
mask2 = np.zeros_like(k_method_data, dtype = bool)

mask1[y1:y2, x1:x2] = True

ind = np.where(k_method_data > cutoff)
mask2[ind] = True

mask = np.logical_and(mask1, mask2)


for filt_dict in filt_list:  

    filt_dict['mask_arr'] = ma.masked_array(filt_dict['data'], mask = ~mask, fill_value = np.nan)
    filt_dict['mask'] = filt_dict['mask_arr'].flatten()

# also grab the PAH continuum subtracted data
pah_3 = {'name': 'F335M'}
pah_7 = {'name': 'F770W'}
pah_11 = {'name': 'F1130W'}

pah_filt_dict = [pah_3 , pah_7, pah_11]

for pah_dict in pah_filt_dict:
    filt = pah_dict['name']
    
    k_method_file = consub_path + f'k_method/SextansA_{filt}_k_method_pah.fits'
    k_method_hdu = fits.open(k_method_file)
    k_method_data_filt = k_method_hdu[0].data
    
    pah_dict['data'] = k_method_data_filt
    pah_dict['mask_arr'] = ma.masked_array(pah_dict['data'], mask = ~mask, fill_value = np.nan)
    pah_dict['mask_arr'][pah_dict['mask_arr'] < 0] = np.nan
    pah_dict['mask'] = pah_dict['mask_arr'].flatten()


plt.figure(4, figsize = (9,6))
plt.clf()


def plot_D21_mods(U_cut=1, color = 'k'):
    
    # load the band ratios from the D21 k analysis for PAHs
    
    D21_sp = ''
    inpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/D21_models/'
    D21_F335M_table = D21_sp + 'F335M_D21_models_with_k_consub.txt'
    D21_F335M = ascii.read(inpath + D21_F335M_table )  

    D21_F770W_table = D21_sp + 'F770W_D21_models_with_k_consub.txt'
    D21_F770W = ascii.read(inpath + D21_F770W_table )  

    D21_F1130W_table = D21_sp + 'F1500W_D21_models_with_k_consub.txt'
    D21_F1130W = ascii.read(inpath + D21_F1130W_table )  

    # cut each table to have U<3
    D21_F335M = D21_F335M[D21_F335M['U'] == U_cut]
    D21_F770W = D21_F770W[D21_F770W['U'] == U_cut]
    D21_F1130W = D21_F1130W[D21_F1130W['U'] == U_cut]
    
    D21_arr = [D21_F335M, D21_F770W, D21_F1130W]
    
    size_arr = ['sma','std', 'lrg']
    ion_arr = ['lo', 'st', 'hi']
    
    prop_F335M = np.zeros((3,3))
    prop_F770W = np.zeros((3,3))
    prop_F1130W = np.zeros((3,3))
    
    prop_arr = [prop_F335M, prop_F770W, prop_F1130W]
    
    prop_label = np.zeros((3,3), dtype = '<U11')
    
    
    i = 0
    for size in size_arr:
        j = 0
        for ion in ion_arr:
    
            ind_F335M = np.where((D21_F335M['size'] == size_arr[i]) & (D21_F335M['ion'] == ion_arr[j]) )
            prop_F335M[i][j] = np.nanmean(D21_F335M[ind_F335M]['pah_k_method'])
            
            ind_F770W = np.where((D21_F770W['size'] == size_arr[i]) & (D21_F770W['ion'] == ion_arr[j]) )
            prop_F770W[i][j] = np.nanmean(D21_F770W[ind_F770W]['pah_k_method'])
                            
            ind_F1130W = np.where((D21_F1130W['size'] == size_arr[i]) & (D21_F1130W['ion'] == ion_arr[j]) )
            prop_F1130W[i][j] = np.nanmean(D21_F1130W[ind_F1130W]['pah_k_method'])
            
            prop_label[i][j] = size_arr[i] + '-' + ion_arr[j]
            
            j = j+1
            
        i = i + 1

    
    # plot the models on top
    plt.scatter(prop_F335M/prop_F770W, prop_F335M/prop_F1130W, color = color, s = 30, alpha = 0.7)
    plt.plot(prop_F335M/prop_F770W, prop_F335M/prop_F1130W, color = color, alpha = 0.7)
    for i in range(3):
        plt.plot(prop_F335M[i][:]/prop_F770W[i][:], prop_F335M[i][:]/prop_F1130W[i][:], color = color, alpha = 0.7)

    # label = False
    
    # if label:
    #     i = 0
    #     for size in size_arr:
    #         j = 0
    #         for ion in ion_arr:
    #             plt.text(prop_F335M[i][j]/prop_F770W[i][j], prop_F335M[i][j]/prop_F1130W[i][j], s = prop_label[i][j])
                
    #             j = j+1
                
    #         i = i + 1



U_list = np.arange(0, 7, 1)

c_list = sns.color_palette("viridis", len(U_list))

for i in range(len(U_list)):

    plot_D21_mods(U_cut=U_list[i], color = c_list[i]) 
 
# fig = plt.gcf()
# fig.subplots_adjust(bottom=0.5)
cmap = mpl.colors.ListedColormap(c_list)
norm = mpl.colors.BoundaryNorm(U_list, cmap.N)
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb2 = mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                ticks=U_list,
                                boundaries=[0] + U_list + [13],
                                spacing='proportional',
                                orientation='vertical')
cb2.set_label('log U, Radiation Field')

# tick_locs = (np.arange(U_list) + 0.5)*(U_list-1)/len(U_list)
# cb2.set_ticks(tick_locs)
   
# norm = simple_norm(pah_3['mask_arr'], 'sqrt', min_cut = 0.002, max_cut = np.nanmax(pah_3['mask_arr']))
# im = plt.scatter(pah_3['mask']/pah_7['mask'], pah_3['mask']/pah_11['mask'], c = pah_3['mask'], cmap = 'magma_r', norm = norm, s = 3, alpha = 0.5)
# cax = plt.colorbar(im)
# cax.set_label('3.3 PAH Flux', size = 'x-large')



ax.set_xlabel('3.3/7.7 PAH', size = 'x-large')
ax.set_ylabel('3.3/11.3 PAH', size = 'x-large')
ax.minorticks_on()

ax.semilogx()
ax.semilogy()

ax.set_xlim(0.03, 4)
ax.set_ylim(0.03, 1.2)

ax.scatter(0.4245467412922441, 0.15613038367480647, marker = '*', color = 'purple', s = 200)

#load the regions with sums
sumpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/'
# sumname = 'SextansA_reg_sum_diffcut.txt'

# sum_table = ascii.read(sumpath + sumname)

# ax.scatter(sum_table['pah_3']/sum_table['pah_7'], sum_table['pah_3']/sum_table['pah_11'], marker = 's', s = 100, color = 'blue')

# for i in range(len(sum_table)):
#     ax.text(sum_table['pah_3'][i]/sum_table['pah_7'][i], sum_table['pah_3'][i]/sum_table['pah_11'][i], s = sum_table['reg'][i])



# sumname = 'SextansA_reg_sum.txt'

# sum_table = ascii.read(sumpath + sumname)

# ax.scatter(sum_table['pah_3']/sum_table['pah_7'], sum_table['pah_3']/sum_table['pah_11'], marker = 's', s = 100, color = 'red')
# for i in range(len(sum_table)):
#     ax.text(sum_table['pah_3'][i]/sum_table['pah_7'][i], sum_table['pah_3'][i]/sum_table['pah_11'][i], s = sum_table['reg'][i])


sumname = 'SextansA_reg_sum_real_cutoff.txt'
sum_table = ascii.read(sumpath + sumname)

ax.scatter(sum_table['pah_3']/sum_table['pah_7'], sum_table['pah_3']/sum_table['pah_11'], marker = 's', s = 100, color = 'red')
# for i in range(len(sum_table)):
    # ax.text(sum_table['pah_3'][i]/sum_table['pah_7'][i], sum_table['pah_3'][i]/sum_table['pah_11'][i], s = sum_table['reg'][i])



savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/band_ratios/'

label= False
if label:
    savename = 'Orig_SextansA_PAH_band_ratio_scatter_D21_model_label.pdf'
else:
    savename = 'Orig_SextansA_PAH_band_ratio_scatter_D21_model.pdf'
    
# plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1)


