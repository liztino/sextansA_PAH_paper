#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates dendrograms of the continuum subtracted PAH emission

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from astrodendro import Dendrogram, pp_catalog
from astrodendro.analysis import PPStatistic
from astropy import units as u
from astropy.wcs import WCS
import matplotlib.patheffects as path_effects
import seaborn as sns
import astropy.visualization.wcsaxes as vis
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle



from astropy.visualization import make_lupton_rgb
from astropy.visualization import make_rgb, ManualInterval, AsinhStretch, ImageNormalize


# try the 7.7 first as an intermediate test
filt_mid = 'F770W'
k = 4.33
datapath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/consub_images/k_method_new/'
dataname = 'SextansA_{:s}_k_method_pah_k_{:3.2f}.fits'.format(filt_mid, k)

pah_hdr = fits.open(datapath + dataname)[0]
pah_data = pah_hdr.data
pah_head = pah_hdr.header

if filt_mid == 'F335M':
	pah_sigma =  0.0036
elif filt_mid == 'F770W':
	pah_sigma =  0.0066
elif filt_mid == 'F1130W':
	pah_sigma = 0.0138

d = Dendrogram.compute(pah_data, min_value = pah_sigma, min_delta = 3*pah_sigma, min_npix = 4)

# load all the other PAH data so we can plot it 
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

pah_3, head_3, err_3, con_3 = load_PAH('F335M', 2.07)
pah_7, head_7, err_7, con_7 = load_PAH('F770W', 4.33)
pah_11, head_11, err_11, con_11 = load_PAH('F1130W', 7.21)


# this is the lupton approach
# it doesn't actually look that good, because 11.3 is the clear brightest so it doesn't show areas in which all 3 PAH detections are in clearly
w = WCS(head_3)
min_vals = [-0.02,-0.02,-0.02]
max_vals = [1.0, 0.5, 0.15]
# image = make_lupton_rgb(pah_11, pah_7, pah_3, stretch=0.1, minimum = min_vals, Q=1)
# intervals = [ManualInterval(vmin=-0.02, vmax = 1.0), 
#              ManualInterval(vmin=-0.02, vmax = 0.5), 
#              ManualInterval(vmin=-0.02, vmax = 0.15)]
intervals = [ManualInterval(vmin=-0.02, vmax = 0.1), 
             ManualInterval(vmin=-0.02, vmax = 0.1), 
             ManualInterval(vmin=-0.02, vmax = 0.1)]
pah_rgb = make_rgb(pah_11, pah_7, pah_3, interval=intervals, stretch=AsinhStretch(a=0.3))


intervals = [ManualInterval(vmin=-0.02, vmax = 0.2), 
             ManualInterval(vmin=-0.02, vmax = 0.2), 
             ManualInterval(vmin=-0.02, vmax = 0.2)]
con_rgb = make_rgb(con_11, con_7, con_3, interval=intervals, stretch=AsinhStretch(a=0.3))
# image = make_rgb(pah_11, pah_7, pah_3, interval=intervals)


fig1 = plt.figure(1, figsize = (8.5,5))
plt.clf()

# plot the continuum filters 
ax1 = fig1.add_subplot(1, 3, 1, projection = w)
im = ax1.imshow(con_rgb, origin = 'lower')

# plot the PAH
ax2 = fig1.add_subplot(1, 3, 2, projection = w)
im = ax2.imshow(pah_rgb, origin = 'lower')

# plot the PAH
ax3 = fig1.add_subplot(1, 3, 3, projection = w)
im = ax3.imshow(pah_rgb, origin = 'lower')

# # load and plot the Halpha
# filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
# # filename = filepath + 'MUSE_SexA_pipeline_Halpha_proj_F1500W'
# filename = filepath + 'Sextans_A_I_Ha_d2009_proj_F1500W'
# Halpha = fits.open(filename + '.fits')

# w_ha = WCS(Halpha[0].header)

# ax4 = fig1.add_subplot(2, 2, 4, projection = w_ha)
# norm = ImageNormalize(Halpha[0].data, stretch=AsinhStretch(), vmin=0, vmax=400)
# cmap = sns.color_palette('magma', as_cmap=True)
# cmap.set_bad(color='black')
# im = ax4.imshow(Halpha[0].data, origin = 'lower', cmap = cmap, norm=norm)


# # load and plot the Halpha
# filepath =  '/Users/etarantino/Documents/JWST_DATA_PAHS/SEXA_ANCILLARY/'
# filename = filepath + 'd75hmrms'
# lt_halpha = fits.open(filename + '.fits')

# w_lt_ha = WCS(lt_halpha[0].header)

# fig4 = plt.figure(4)
# plt.clf()
# ax_lt = fig4.add_subplot(1,1,1, projection = w_lt_ha)
# ax_lt.imshow(lt_halpha[0].data, origin = 'lower')

# finding the mask for each clump 
# 	# get the clumps that are only in the previous catalog with detections in all PAH features
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_file = 'template_clump_F770W.csv'

clumps = ascii.read(clump_path + clump_file, format = 'csv')
filt_idx = '{:s}_k_{:3.2f}'.format(filt_mid, k)
clump_list = clumps[filt_idx]

# load the final clump table, since it has the center
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_file = 'clump_props_pah_filt_flux_corrected.txt'
clumps = ascii.read(clump_path + clump_file)


# now try a plot the clump 
wcs = WCS(pah_head)

# plt.figure(5)
# plt.clf()
# fig, ax = plt.subplots(subplot_kw=dict(projection=wcs), num = 5)
# ax.imshow(pah_data, origin = 'lower', vmin = -0.05, vmax = 0.1)


axes = [ax1, ax2, ax3]

# loop through the clumps
for i in range(len(clump_list)):  
    
    # get structure from dendrogram stuff
    structure = d[clump_list[i]]
    mask = structure.get_mask()
    mask_hdu = fits.PrimaryHDU(mask.astype('short'), pah_head)


    for ax in axes:

        # plot number    
        if clumps['clump_num'][i] == 15:
            # plot contour 
            cs = ax.contour(mask_hdu.data, colors = 'gainsboro', linewidths = 0.5) 
            
            # ax.text(clumps['ra_cen'][i]-0.0001, clumps['dec_cen'][i]-0.0008, s = '{:d}'.format(clumps['clump_num'][i]), transform=ax.get_transform('fk5'), c = 'k', fontsize = 'medium', 
            #         fontweight='bold', path_effects=[path_effects.withStroke(linewidth=1, foreground='white')])
            
            
        elif clumps['clump_num'][i] == 16:
            # plot contour 
            cs = ax.contour(mask_hdu.data, colors = 'gainsboro', linewidths = 0.5) 
            # ax.text(clumps['ra_cen'][i]-0.0002, clumps['dec_cen'][i]-0.0005, s = '{:d}'.format(clumps['clump_num'][i]), transform=ax.get_transform('fk5'), c = 'k', fontsize = 'medium', 
            #         fontweight='bold', path_effects=[path_effects.withStroke(linewidth=1, foreground='white')])
            
        elif clumps['clump_num'][i] == 1:
            # plot contour 
            cs = ax.contour(mask_hdu.data, colors = 'gainsboro', linewidths = 0.5) 
            # ax.text(clumps['ra_cen'][i]+0.0006, clumps['dec_cen'][i]+0.0001, s = '{:d}'.format(clumps['clump_num'][i]), transform=ax.get_transform('fk5'), c = 'k', fontsize = 'medium', 
            #         fontweight='bold', path_effects=[path_effects.withStroke(linewidth=1, foreground='white')])
            
        elif clumps['clump_num'][i] == 2:
            # plot contour 
            cs = ax.contour(mask_hdu.data, colors = 'gainsboro', linewidths = 0.5) 
            # ax.text(clumps['ra_cen'][i]+0.00025, clumps['dec_cen'][i]+0.0002, s = '{:d}'.format(clumps['clump_num'][i]), transform=ax.get_transform('fk5'), c = 'k', fontsize = 'medium', 
            #         fontweight='bold', path_effects=[path_effects.withStroke(linewidth=1, foreground='white')])
            
        else:
            # plot contour 
            cs = ax.contour(mask_hdu.data, colors = 'gainsboro', linewidths = 0.5) 
            # ax.text(clumps['ra_cen'][i]+0.0006, clumps['dec_cen'][i]+0.0003, s = '{:d}'.format(clumps['clump_num'][i]), transform=ax.get_transform('fk5'), c = 'k', fontsize = 'medium', 
            #         fontweight='bold', path_effects=[path_effects.withStroke(linewidth=1, foreground='white')])
    
    
for ax in axes:
    
    xlim1 = 700
    xlim2 = 1150
    ylim1 = 400
    ylim2 = 912
    
    # PSF/beam info 
    # in arcseconds
    # if ax == ax4:
    #     fwhm = (1.1 * u.arcsec).to(u.deg)
    # else:
    #     fwhm = (0.488 * u.arcsec).to(u.deg)
    # beam_angle = 0 * u.deg
    
    ax.set_xlim(xlim1, xlim2)
    ax.set_ylim(ylim1, ylim2)
    
    # # add scalebar and such 
    # vis.add_beam(ax = ax, header = None, major = fwhm,  minor = fwhm, angle = beam_angle, color = 'k', frame = True)
    # vis.add_scalebar(ax = ax, length = 5 * u.arcsec, corner = 'top right', label = '5\" = 35 pc', color = 'white')
    
    # # add the Shi region
    # shi = SkyCoord('10:11:06.2', '-04:42:23', unit=(u.hourangle, u.deg))
    # s = SphericalCircle(shi, (32/2) * u.arcsec, edgecolor = 'white', facecolor = 'none',  transform=ax.get_transform('fk5'), alpha = 0.5)
    # ax.add_patch(s)



# load the YSOs
anc_path = '/Users/etarantino/Documents/JWST/01619/YSOs/'
cat_file = 'ysos_with_colorcuts.txt'
cat = ascii.read(anc_path + cat_file)
c = SkyCoord(cat['RA'], cat['DEC'], frame = 'fk5', unit = (u.deg, u.deg))

ax3.scatter(c.ra.deg, c.dec.deg, marker = 'x', c = 'mediumturquoise', transform = ax.get_transform('world'), s = 20, zorder = 0.0001)


# # set limits for for the Halpha map
# # first get the RA and DEC for the regular
# RA1, DEC1 = w.wcs_pix2world((xlim1, ylim1))
# RA2, DEC2 = w.wcs_pix2world((xlim2, ylim2))

# # now translate to pixels
# (xlim1_ha, ylim1_ha) = w_ha.wcs_world2pix((RA1, DEC1))
# (xlim2_ha, ylim2_ha) = w_ha.wcs_world2pix((RA2, DEC2))

# ax4.set_xlim(xlim1_ha, xlim2_ha)
# ax4.set_ylim(ylim1_ha, ylim2_ha)


ax1.set_ylabel('DEC')
# ax1.yaxis.set_visible(True)
# ax1.tick_params(axis='y', labelleft=True)

ax3.set_ylabel('DEC')
# ax1.yaxis.set_visible(True)
# ax1.tick_params(axis='y', labelleft=True)

ax1.set_xlabel('RA')

ax2.set_xlabel('RA')
ax3.set_xlabel('RA')

# ax4.set_xlabel('RA')


ax2.yaxis.set_visible(False)
ax2.tick_params(axis='y', labelleft=False)

ax3.yaxis.set_visible(False)
ax3.tick_params(axis='y', labelleft=False)

# ax1.xaxis.set_visible(False)
# ax1.tick_params(axis='x', labelbottom=False)

# ax2.xaxis.set_visible(False)
# ax2.tick_params(axis='x', labelbottom=False)


plt.subplots_adjust(wspace = 0.1, hspace = 0.05)


ax1.set_title('a.) Continuum Images')
ax2.set_title('b.) Continuum-Subtracted PAH Images')
ax3.set_title('c.) YSO candidates')
# ax4.set_title('d.) Ionized Gas H$\\alpha$')

savepath = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
savename = 'clumps_on_RGB_PAH_with_YSOs.png'
plt.savefig(savepath + savename, bbox_inches='tight',pad_inches = 0.1, dpi = 400)



