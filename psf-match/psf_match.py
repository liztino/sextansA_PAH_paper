#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

PSF matches from given PSF files

@author: Eliz
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
import scipy.fftpack as fp
import scipy.ndimage as nd
import os
from astropy.convolution import convolve_fft
from astropy.io import ascii
from circularize_laura import circularize
from center import center
from fwhm import fwhm
from regrid import regrid
from radial_data import radial_data
from shift_psf import shift_psf


def psf_match(psf1_file, psf2_file, name):

    # # psf1_file = 'SIII-19_psf'  
    # psf1_file = 'SIV-11_psf'
    # # psf2_file = 'SIII-33_psf'
    # psf2_file = 'SiII-35_psf'
        
    print('p1')
    # # path for the psf files
    psf_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/psfs/'

    # open psf files
    psf1_hdu = fits.open(psf_path + psf1_file + '.fits')
    psf1 = psf1_hdu[0].data
    psf1_head = psf1_hdu[0].header
    
    psf2_hdu = fits.open(psf_path + psf2_file + '.fits')
    psf2 = psf2_hdu[0].data
    psf2_head = psf2_hdu[0].header
    
    # input psf
    # psf1_name = 'SIII_19_v3.fits'
    # # psf1_name = 'SIII_19_10arcmin_samp2.fits'
    # # psf1_name = 'SIV_10_15arcmin.fits'
    
    
    # # desired psf
    # psf2_name = 'SIII_33_v1.fits'
    # # psf2_name = 'SIII_33_10arcmin_samp2.fits'
    # # psf2_name = 'SIII_33_15arcmin.fits'
    
    
    # filepath = '/Users/Eliz/Documents/UMD/Research/stinytim/test/'
    
    # psf1_hdu = fits.open(filepath + psf1_name)
    # psf1 = psf1_hdu[0].data
    # psf1_head = psf1_hdu[0].header
    
    # psf2_hdu = fits.open(filepath + psf2_name)
    # psf2 = psf2_hdu[0].data
    # psf2_head = psf2_hdu[0].header
    
    # identifies names of lines
    psf2_name = psf2_file.split('.')[0]
    
    
    # create figure path
    psf_name = name + '_to_' + psf2_name
    figpath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/kernels/' + psf_name + '/'
    if not os.path.exists(figpath):
        os.makedirs(figpath)
        
    print(np.sum(psf1))
    print(np.sum(psf2))
        
    # pix_size = psf2_head['PIXSCALX']
    

    # regrid each psf 
    pix_size = psf1_head['PIXELSCL']
    # psf1 = regrid(psf1, psf1_head['PIXSCALX'], pix_size)
    psf2 = regrid(psf2, psf2_head['PIXELSCL'], pix_size)
    
    print(len(psf1))
    print(len(psf2))
        
    print(center(psf1))
    print(center(psf2))
    
    print(psf1_head['PIXELSCL'])
    print(psf2_head['PIXELSCL'])

    
    # grid_size = 2855
    # resize to grid size of 3645 pixels or 729 arcsec
    def resize(img, grid_size):    
        # trim if greater than grid size
        if np.shape(img)[0] > grid_size:
            xcen, ycen = center(img)
            
            ind = int((grid_size)/2)
            img = img[(xcen-ind):(xcen+ind+1), (ycen-ind):(ycen+ind+1)]
            
            return img
         
        # pad with zeros if less than grid dize
        elif np.shape(img)[0] < grid_size:
            shape = np.shape(img)
            
            pad_val = grid_size - shape[0]
            print(pad_val)
            
            new_array = np.zeros((grid_size, grid_size), dtype = float)
            new_shape = np.shape(new_array)
                
            xstart_ker = int((new_shape[0] - shape[0])/2.) + 1 
            xend_ker = xstart_ker + shape[0]
            ystart_ker = int((new_shape[1] - shape[1])/2.) + 1
            yend_ker = ystart_ker + shape[1]
            
            print(xstart_ker, xend_ker, ystart_ker, yend_ker)
            
            new_array[xstart_ker:xend_ker, ystart_ker:yend_ker] = img
            
            return new_array
        
        # do nothing if shape is grid size
        else:
            return img
        
    # grid_size = len(psf2)
    # grid_size = 3645
    # grid_size = 1859
    
    # psf1 = resize(psf1, grid_size)
    # psf2 = resize(psf2, grid_size)

    
    if (len(psf1) >= len(psf2)):
        grid_size = len(psf1)
        psf2 = resize(psf2, grid_size)
    else:
        grid_size = len(psf2)
        psf1 = resize(psf1, grid_size)
    
    print(np.sum(psf1))
    print(np.sum(psf2))
    
    print(len(psf1))
    print(len(psf2))
    
    print(center(psf1))
    print(center(psf2))
    
    # normalize psfs to one
    psf1 = psf1 / np.sum(psf1)
    psf2 = psf2 / np.sum(psf2)
    
    # plot initial PSF
    vmin = 0
    vmax = 0.01 * np.nanmax(psf1)
    norm = simple_norm(psf1, 'log', min_cut = vmin, max_cut = vmax)
    
    fig1 = plt.figure(1)
    plt.clf() 
    ax1 = fig1.add_subplot(1,2,1)
    ax1.imshow(psf1, norm = norm, origin = 'lower')
    ax1.title.set_text('Input PSF1')
    
    ax2 = fig1.add_subplot(1,2,2)
    ax2.imshow(psf2, norm = norm, origin = 'lower')
    ax2.title.set_text('Desired PSF2')
        
    figname = 'orig_psfs'
    plt.savefig(figpath + figname + '.pdf')
    
    # circularize and plot circularized psfs
    # psf1 = circularize(psf1)
    psf1 = psf1 / np.sum(psf1)
    
    # psf2 = circularize(psf2)
    psf2 = psf2 / np.sum(psf2)
    
    fig1 = plt.figure(13)
    plt.clf() 
    ax1 = fig1.add_subplot(1,2,1)
    ax1.imshow(psf1, norm = norm, origin = 'lower')
    ax1.title.set_text('Input PSF1')
    
    ax2 = fig1.add_subplot(1,2,2)
    ax2.imshow(psf2, norm = norm, origin = 'lower')
    ax2.title.set_text('Desired PSF2')
        
    figname = 'circ_psfs'
    plt.savefig(figpath + figname + '.pdf')
    
    # psf1 = shift_psf(psf1)
    # psf2 = shift_psf(psf2)
    
    # psf1 = nd.shift(psf1, 1)
    # psf2 = nd.shift(psf2, 1)
    
    # vmin = 0
    # vmax = 0.1 * np.nanmax(psf1)
    # norm = simple_norm(psf1, 'log', min_cut = vmin, max_cut = vmax)
    
    # fig2 = plt.figure(2)
    # plt.clf() 
    # ax1 = fig2.add_subplot(1,2,1)
    # ax1.imshow(psf1, norm = norm, origin = 'lower')
    # ax1.title.set_text('Input PSF1')
    
    # ax2 = fig2.add_subplot(1,2,2)
    # ax2.imshow(psf2, norm = norm, origin = 'lower')
    # ax2.title.set_text('Desired PSF2')
    # figname = 'circ_psfs'
    # plt.savefig(figpath + figname + '.pdf')
    
    # perform the FFT
    fft_psf1 = fp.fft2(psf1)
    fft_psf2 = fp.fft2(psf2)
    
    # plotting FFT results
    mag_psf1 = np.sqrt(fft_psf1.real**2 + fft_psf1.imag**2)
    mag_psf2 = np.sqrt(fft_psf2.real**2 + fft_psf2.imag**2)
    
    vmin = 0
    vmax = 0.1 * np.nanmax(mag_psf1)
    norm = simple_norm(mag_psf1, 'log', min_cut = vmin, max_cut = vmax)
    
    fig3 = plt.figure(3)
    plt.clf()
    ax1 = fig3.add_subplot(1,2,1)
    ax1.imshow(fp.fftshift(mag_psf1), origin = 'lower', norm = norm)
    ax1.title.set_text('Magnitude of PSF1 FFT ')
    
    ax2 = fig3.add_subplot(1,2,2)
    ax2.imshow(fp.fftshift(mag_psf2), origin = 'lower', norm = norm)
    ax2.title.set_text('Magnitude of PSF2 FFT ')
    
    figname = 'ft_mag'
    plt.savefig(figpath + figname + '.pdf')
    
    sh_fft_psf1 = fp.fftshift(fft_psf1)
    sh_fft_psf2 = fp.fftshift(fft_psf2)
    
    phi_psf1 = np.arctan(sh_fft_psf1.imag/sh_fft_psf1.real)
    phi_psf2 = np.arctan(sh_fft_psf2.imag/sh_fft_psf2.real)
    
    norm = simple_norm(phi_psf1, 'log')
    
    fig4 = plt.figure(4)
    plt.clf()
    ax1 = fig4.add_subplot(1,2,1)
    ax1.imshow(phi_psf1, origin = 'lower')
    ax1.title.set_text('Phase of PSF1 FFT ')
    
    ax2 = fig4.add_subplot(1,2,2)
    im = ax2.imshow(phi_psf2, origin = 'lower')
    ax2.title.set_text('Phase of PSF2 FFT ')
    
    cbaxes = fig4.add_axes([0.15, 0.08, 0.7, 0.05]) 
    plt.colorbar(im, cax = cbaxes, orientation='horizontal')
    
    figname = 'ft_phase'
    plt.savefig(figpath + figname + '.pdf')
    
    # First filter for both psfs 
    xfreq = fp.fftfreq(np.shape(fft_psf1)[0], pix_size) * 2 * np.pi
    yfreq = fp.fftfreq(np.shape(fft_psf1)[1], pix_size) * 2 * np.pi

    filt = np.zeros_like(fft_psf1)
    freq = np.zeros_like(fft_psf1)
    
    width1 = fwhm(psf1, psf1_head)
    fwhm1 = width1['xarcsec']
    
    kappa_list = np.arange(0.8, 1.25, 0.05)
    kappa_list =[1]
    # kappa_list = [0.95]
    W_list = np.zeros(len(kappa_list))
    D_list = np.zeros(len(kappa_list))
    
    
    for i in range(len(kappa_list)):
        kappa = kappa_list[i]
    
        kappapath = figpath + 'k_{:4.3f}/'.format(kappa)
        if not os.path.exists(kappapath):
            os.makedirs(kappapath)
    
        k_h = kappa * 2*np.pi /fwhm1
        k_l = 0.7*k_h
        
        for x in range(len(xfreq)):
            for y in range(len(yfreq)):
                freq[x,y] = np.sqrt(xfreq[x]**2 + yfreq[y]**2)
                
                    
                # main filter for psf1
                if freq[x,y] <= k_l:
                    filt[x,y] = 1.
                    
                if ((freq[x,y] <= k_h) and (freq[x,y]>=k_l)):
                    filt[x,y] = 0.5 * (1 + np.cos(np.pi * ((freq[x,y] - k_l)/(k_h - k_l))))
                    
                if freq[x,y] >= k_h:
                    filt[x,y] = 0.
        
        div = abs((1./fft_psf1)*filt)
        
        vmin = -0.1* np.nanmax(div)
        vmax = 0.8 * np.nanmax(div)
        norm = simple_norm(div, 'log', min_cut = 0, max_cut = vmax)
        
        fig6 = plt.figure(6)
        plt.clf()
        ax1 = fig6.add_subplot(1,2,1)
        ax1.imshow(fp.fftshift(abs(filt)), origin = 'lower')
        ax1.title.set_text('Filt PSF1 k_h = {:4.3f}'.format(k_h))
        
        ax2 = fig6.add_subplot(1,2,2)
        im = ax2.imshow(fp.fftshift(abs((1/fft_psf1)*filt)) , origin = 'lower', norm = norm)
        ax2.title.set_text('FFT(PSF1) * phi')
        
        figname = 'psf1_filter'
        plt.savefig(kappapath + figname + '.pdf')
        
        # create kernel
        ker = fp.ifft2((fft_psf2/fft_psf1) * filt)
        ker = fp.fftshift(ker.real)
        # ker = circularize(ker)
        
        # print(np.sum(ker))
        
        ker = ker/np.sum(ker)
        
        neg_sum = np.sum(ker[ker < 0])
        
        print(ker)
        
        plt.figure(8)
        plt.clf()
        plt.imshow(ker)
                
        vmin = -0.1* np.nanmax(ker)
        vmax = 0.1 * np.nanmax(ker)
        norm = simple_norm(ker, 'linear', min_cut = vmin, max_cut = vmax)
        
        xcen_k, ycen_k = center(ker)
        
        plt.figure(8)
        plt.clf()
        plt.imshow(ker, norm = norm, origin = 'lower')
        plt.title('Kernel')
        # plt.xlim(xcen_k - 50, xcen_k + 50)
        # plt.ylim(ycen_k - 50, ycen_k + 50)
        
        figname = 'kernel'
        plt.savefig(kappapath + figname + '.pdf')
        
        conv = convolve_fft(psf1, ker, boundary='wrap', allow_huge=True)
        
        psf2_rad = radial_data(psf2)
        conv_rad = radial_data(conv)
        
        psf2_r = psf2_rad.r * pix_size
        conv_r = conv_rad.r * pix_size
        
        
        plt.figure(9)
        plt.clf()
        plt.plot(psf2_r, psf2_rad.median, label = 'Orig PSF2')
        plt.plot(conv_r, conv_rad.median, ls = '--', label = 'Test PSF2')
        plt.legend(loc = 'best')
        plt.ylabel(r'PSF ($\Psi(\Theta)$)')
        plt.xlabel(r'$\Theta$ (arcsec)')
        
        plt.xlim(0,2)
        
        figname = 'rad_prof'
        plt.savefig(kappapath + figname + '.pdf')
        
        rad1 = psf2_r  * psf2_rad.median / (np.max(psf2_r  * psf2_rad.median))
        rad2 = conv_r * conv_rad.median / (np.max(conv_r * conv_rad.median))
        
        D = np.sum(np.abs(conv - psf2))
        
        print('kappa = {:4.3f}'.format(kappa))
        print('D = {:10.8f}'.format(D))
        print('W- = {:10.8f}'.format(neg_sum))
        
        plt.figure(10)
        plt.clf()
        plt.plot(psf2_r , rad1, label = 'Orig PSF2')
        plt.plot(conv_r, rad2, ls = '--', label = 'Test PSF2')
        plt.legend(loc = 'best')
        plt.ylabel(r'$\Theta \ \Psi(\Theta) / max(\Theta \ \Psi)$')
        plt.xlabel(r'$\Theta$ (arcsec)')
        
        plt.xlim(0,2)
        
        ax = plt.gca()
        text = 'kappa= {:4.3f} \nD = {:8.7f} \nW- = {:8.7f}'.format(kappa, D, neg_sum)
        plt.text(0.70, 0.70,text,transform = ax.transAxes)
        
        figname = 'aniano_prof'
        plt.savefig(kappapath + figname + '.pdf')
        
        diff_im = conv-psf2
        
        vmin = -0.1* np.nanmax(diff_im)
        vmax = 0.1 * np.nanmax(diff_im)
        norm = simple_norm(diff_im, 'linear', min_cut = vmin, max_cut = vmax)
        
        plt.figure(11)
        plt.clf()
        im = plt.imshow(diff_im, norm = norm, origin = 'lower')
        plt.colorbar(im)
        plt.title('Convolved - PSF2')
        
        figname = 'diff_image'
        plt.savefig(kappapath + figname + '.pdf')
        
        
        psf1_head['NAXIS1'] = len(ker)
        psf1_head['NAXIS2'] = len(ker)
        
        kername = psf_name + '_ker_{:4.2f}'.format(kappa)
        
        h = fits.PrimaryHDU(ker, header = psf1_head)
        
        h.writeto(kappapath + kername + '.fits', overwrite = True)
        
        D_list[i] = D
        W_list[i] = neg_sum
        
    print(kappa_list)
    print(D_list)
    print(W_list)
        
    W_list = abs(W_list)
    
    plt.figure(12)
    plt.clf()
    ax = plt.gca()
    ax.plot(kappa_list, D_list, '.', c = 'blue')
    lns1 = ax.plot(kappa_list, D_list, '-', label = 'D', c = 'blue')
    ax.set_xlabel('kappa')
    ax.set_ylabel('D')
    ax.set_ylim(0,0.1)
    # ax.semilogy()
    
    ax2 = ax.twinx()
    ax2.plot(kappa_list, W_list, '.', c = 'orange')
    lns2 = ax2.plot(kappa_list, W_list, '-', label = 'W-', c = 'orange')
    # ax2.semilogy()
    ax2.set_ylabel('W-')
    ax2.set_ylim(0,1.5)
    
    
    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='best')
    
    
    figname = 'kappa'
    plt.savefig(figpath + figname + '.pdf')
    
    savedata = np.transpose(np.array([kappa_list, D_list, W_list]))
    savename = 'kappa_results.txt'
    
    header = 'kappa \t D \t W-'
    np.savetxt(figpath + savename, savedata, header = header, delimiter = '\t')

# psf1_file = 'SIII_33-33_psf'
# psf2_file = 'SiII-35_psf'
# name = 'name'
# psf_match(psf1_file, psf2_file, name)

# line_path = '/Users/Eliz/Documents/UMD/Research/pahfit/'
# line_list = 'SL1_lines.txt'
# lines = ascii.read(line_path + line_list, names = ['names', 'wave'])
# lines = lines[::-1]

# for i in range(len(lines)):
# psf1_file = lines['names'][i] + '-' + str(round(lines['wave'][i])) + '_psf'

# psf1_file = 'F335M_webbpsfv130_rot'
# psf2_file = 'F360M_webbpsfv130_rot'
# psf_match(psf1_file, psf2_file, psf1_file)
    
# name = lines['names'][i] + '-' + str(round(lines['wave'][i]))
# print('Working on ' + psf1_file)
# psf_match(psf1_file, psf2_file, psf1_file)

# filts = ['F560W',  'F770W', 'F1000W', 'F1130W']
# filts = ['F1130W']
# filts = ['F300M', 'F335M', 'F360M', 'F560W',  'F770W', 'F1000W', 'F1130W']
filts = ['F150W', 'F200W']
# filts = [ 'F560W',  'F770W', 'F1000W', 'F1130W']



# # # val = np.arange(5,14,1)
for filt in filts:
    
    psf1_file = f'{filt}_webbpsfv130_rot'
    psf2_file = 'F1500W_webbpsfv130_rot'
        
    # name = lines['names'][i] + '-' + str(round(lines['wave'][i]))
    print('Working on ' + psf1_file)
    psf_match(psf1_file, psf2_file, psf1_file)