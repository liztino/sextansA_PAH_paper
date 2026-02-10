#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Reprojects images onto each other 

@author: ejtino
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from reproject import reproject_interp, reproject_adaptive, reproject_exact
from astropy.table import Table
from scipy.interpolate import interp1d

def regrid_cube(order, file):

    sm_path = '/Users/Eliz/Documents/UMD/Research/stinytim/cube/sm_cube/'
    reg = 'n76'    
    v = ''
        
    cubename = reg + '_' + order + '_sm_to_irs_35.fits'
    uncname = reg + '_' + order + '_sm_to_irs_35_unc.fits'
    
    cubehdu = fits.open(sm_path + cubename)
    cube = cubehdu[0].data
    cube_head = cubehdu[0].header
    
    unchdu = fits.open(sm_path + uncname)
    unc = unchdu[0].data
    unc_head = unchdu[0].header

    cubepath = '/Users/Eliz/Documents/UMD/Research/cubism/cubes/' + reg + '/'
    cubename = reg + '_' + order + '_cube_2020' + v + '.fits'
    wavehdu = fits.open(cubepath + cubename)
    wave_data = wavehdu[1].data
    wave = wave_data[0][0]
    wave = wave.flatten()
    
    print(wave[0], wave[-1])
    
    # file_path = '/Users/Eliz/Documents/UMD/Research/pahfit/linemaps/fits/n76/' + order + '/' + file
    # file_hdu = fits.open(file_path)
    # file_head = file_hdu[0].header
    
    # LL1_path = '/Users/Eliz/Documents/UMD/Research/pahfit/linemaps/fits/n76/LL1/n76_LL1_SiII.fits'
    # LL1_hdu = fits.open(LL1_path)
    # LL1_head = LL1_hdu[0].header
    
    # shape = np.shape(LL1_hdu[0].data)
    
    # # print(np.shape(file_hdu[0].data))
    # # print(np.shape(cube))
    # # print(np.shape(LL1_hdu[0].data))
    
    
    # l = np.shape(cube)[0]
    # new_cube = np.zeros((l, shape[0], shape[1]))
    # new_unc = np.zeros((l, shape[0], shape[1]))
    
    # for i in range(l):
    
    #     new_cube[i,:,:] = reproject_interp((cube[i,:,:], file_head), LL1_head, return_footprint=False)
    #     new_unc[i,:,:] = reproject_interp((unc[i,:,:], file_head), LL1_head, return_footprint=False)

    # pri = fits.PrimaryHDU(new_cube, header = LL1_head)
    # pri.add_checksum()
    
    # # tab = fits.BinTableHDU(wavehdu[1].data)
    # tab = Table(np.transpose(wave))
    # bintab = fits.table_to_hdu(tab)
    
    # hdu_list = fits.HDUList([pri, bintab])
    
    # savepath = '/Users/Eliz/Documents/UMD/Research/stinytim/cube/regrid_cube/'
    # savename = reg + '_' + order + '_regrid.fits'

    # hdu_list.writeto(savepath + savename, overwrite = True)    

    # return wave, new_cube, new_unc, LL1_head

# LL1_wave, LL1_cube, LL1_unc, LL1_head = regrid_cube('LL1', 'n76_LL1_SiII.fits')
# LL2_wave, LL2_cube, LL2_unc, LL2_head = regrid_cube('LL2', 'n76_LL2_NeIII.fits')
# SL1_wave, SL1_cube, SL1_unc, SL1_head = regrid_cube('SL1', 'n76_SL1_SIV.fits')
# SL2_wave, SL2_cube, SL2_unc, SL2_head = regrid_cube('SL2', 'n76_SL2_H2_S5.fits')

regrid_cube('LL1', 'n76_LL1_SiII.fits')
regrid_cube('LL2', 'n76_LL2_NeIII.fits')
regrid_cube('SL1', 'n76_SL1_SIV.fits')
regrid_cube('SL2', 'n76_SL2_H2_S5.fits')

# x = 75
# y = 25

# plt.figure(1)
# plt.clf()
# plt.errorbar(LL1_wave, LL1_cube[:, y, x], yerr = LL1_unc[:,y,x], label = 'LL1')
# plt.errorbar(LL2_wave, LL2_cube[:, y, x], yerr = LL2_unc[:,y,x], label = 'LL2')
# plt.errorbar(SL1_wave, SL1_cube[:, y, x], yerr = SL1_unc[:,y,x], label = 'SL1')
# plt.errorbar(SL2_wave, SL2_cube[:, y, x], yerr = SL2_unc[:,y,x], label = 'SL2')
# plt.legend(loc = 'best')

# # finding overlap between wave arrays and creating a final wave array

# ############## SL1 & SL2 #############
# def get_overlap(wave1, cube1, wave2, cube2):
    
#     # wave1 = SL2_wave
#     # cube1 = SL2_cube
#     # wave2 = SL1_wave
#     # cube2 = SL2_cube
    
#     ovr1 = np.where(wave1 > wave2[0])[0]
#     ind1 = ovr1[0] - 2
#     ind2 = ovr1[-1] + 1
    
#     ovr2 = np.where(wave2 < wave1[-1])[0]
    
#     # print(ovr2)
    
#     wave_ovr = np.zeros_like(cube2[ovr2,:,:])
#     wave_vals = wave2[ovr2]
    
#     for x in range(np.shape(LL1_cube)[2]):
#         for y in range(np.shape(LL1_cube)[1]):
#             s = interp1d(wave1[ind1:ind2], cube1[ind1:ind2,y,x], kind = 'cubic') 
            
#             wave_ovr[ovr2, y, x] = np.nanmean([s(wave_vals),cube2[ovr2, y, x]])
            
    
#     len1 = np.shape(cube1)[0]
#     len2 = np.shape(cube2)[0]
    
#     length = len1 + len2 - len(ovr2) - len(ovr1) + len(wave_vals)
    
#     ovr_cube = np.zeros((length,np.shape(cube1)[1], np.shape(cube1)[2]))
#     ovr_wave = np.zeros(length)
    
#     ovr_cube[:ovr1[0], :, :] = cube1[:ovr1[0], :,:]
#     ovr_cube[ovr1[0]:ovr1[0]+len(ovr2), :, :] = wave_ovr[:,:,:]
#     ovr_cube[ovr1[0]+len(ovr2):, :, :] = cube2[len(ovr2):, :, :]

#     ovr_wave[:ovr1[0]] = wave1[:ovr1[0]]
#     ovr_wave[ovr1[0]:ovr1[0]+len(ovr2)] = wave_vals
#     ovr_wave[ovr1[0]+len(ovr2):] = wave2[len(ovr2):]

#     return ovr_wave, ovr_cube


# SL_wave, SL_unc = get_overlap(SL2_wave, SL2_unc, SL1_wave, SL1_unc)
# SL_LL_wave, SL_LL_unc = get_overlap(SL_wave, SL_unc, LL2_wave, LL2_unc)

# final_wave, final_unc = get_overlap(SL_LL_wave, SL_LL_unc, LL1_wave, LL1_unc)

# del LL1_head['LINE']
# del LL1_head['LINEWAV']

# pri = fits.PrimaryHDU(final_unc, header = LL1_head)
# pri.add_checksum()

# # tab = fits.BinTableHDU(wavehdu[1].data)
# tab = Table(np.transpose(final_wave))
# bintab = fits.table_to_hdu(tab)

# hdu_list = fits.HDUList([pri, bintab])

# savepath = '/Users/Eliz/Documents/UMD/Research/stinytim/cube/regrid_cube/'
# savename = 'n76_all_orders_regrid_unc.fits'

# hdu_list.writeto(savepath + savename, overwrite = True)    

        