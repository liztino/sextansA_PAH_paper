#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate radial profiles of the more extended clumps to determine how the F1500W and PAH fluxes vary as a function of distance

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

from astropy.visualization import simple_norm
from astropy.modeling.models import Gaussian2D
from photutils.centroids import centroid_quadratic
from photutils.datasets import make_noise_image
from photutils.profiles import RadialProfile

from get_pivot_wave import get_pivot_wave

plt.ion()

# load clump table
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_file = 'clump_props_pah_filt_flux_corrected.txt'
clumps = ascii.read(clump_path + clump_file)

clump = 14
pix_rad = 30


# load all the filters and PAH
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

# load the PAH now
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

pah_3_k1, head_3_k1, err_3_k1, con_3_k1 = load_PAH('F335M', 2.07)
pah_7_k1, head_7_k1, err_7_k1, con_7_k1 = load_PAH('F770W', 4.33)
pah_11_k1, head_11_k1, err_11_k1, con_11_k1 = load_PAH('F1130W', 7.21)
pah3_sigma = 0.0036
pah7_sigma =  0.0066
pah11_sigma = 0.0138
F1500W_sigma = 0.020477

# create a 2D array for the F1500W error
F1500W_err = np.ones_like(F1500W['data'])*F1500W_sigma

# set up the profile arrays
xycen = (round(clumps[clump]['x_cen']), round(clumps[clump]['y_cen']))
edge_radii = np.arange(pix_rad)
beam_fwhm = 0.488

# get pixel scale to convert to arcseconds
pixscale = (head_11_k1['CDELT2'] * u.deg).to(u.arcsec).value

# try using the photutils radial profile generator
rp = RadialProfile(pah_11_k1, xycen, edge_radii, error=err_11_k1)
rp_F1500W = RadialProfile(F1500W['data'], xycen, edge_radii, error = F1500W_err)
rp_con = RadialProfile(con_11_k1, xycen, edge_radii, error = err_11_k1)

rp_max_val = np.nanmax(rp.profile)
rp_max_val_F1500W = np.nanmax(rp_F1500W.profile)

# get normalization scaling factors
rp_max = np.nanmax(rp.profile)
rp_con_max = np.nanmax(rp_con.profile)
rp_F1500W_max = np.nanmax(rp_F1500W.profile)

# normalize
rp.normalize(method = 'sum')
rp_F1500W.normalize(method = 'sum')
rp_con.normalize(method = 'sum')

# get normalization scaling factors
rp_max_norm = np.nanmax(rp.profile)
rp_con_max_norm = np.nanmax(rp_con.profile)
rp_F1500W_max_norm = np.nanmax(rp_F1500W.profile)

plt.figure(4, figsize=(15, 5))
plt.clf()
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5), num = 4)
ax1.plot(rp.radius * pixscale, rp.profile, zorder = 5, c = '#0072B2', label = '11.3 PAH')
ax1.axhline(pah11_sigma*(rp_max_norm/rp_max),  c = '#0072B2', lw = 0.5, ls = '--')
ax1.fill_between(rp.radius * pixscale, rp.profile + rp.profile_error, rp.profile - rp.profile_error, zorder = 0.01, color = '#0072B2', alpha = 0.1)

ax1.plot(rp_con.radius * pixscale, rp_con.profile, zorder = 5, c = '#009E73', label = '11.3 CON')
ax1.fill_between(rp_con.radius * pixscale, rp_con.profile + rp_con.profile_error,  rp_con.profile - rp_con.profile_error, zorder = 0.01, color = '#33B58F', alpha = 0.1)


ax1.plot(rp_F1500W.radius * pixscale, rp_F1500W.profile, zorder = 5, c = '#E69F00', label = 'F1500W')
ax1.axhline(F1500W_sigma*(rp_F1500W_max_norm/rp_F1500W_max),  c = '#E69F00', lw = 0.5, ls = '--')
ax1.fill_between(rp_F1500W.radius * pixscale, rp_F1500W.profile + rp_F1500W.profile_error, rp_F1500W.profile - rp_F1500W.profile_error, zorder = 0.01, color = '#F0B84A', alpha = 0.1)

ax1.set_xlim(0, edge_radii[-1] * pixscale)
ax1.axvline(beam_fwhm, c = 'k', lw = 0.5, ls = '--', zorder = 0.000001)

ax1.legend(loc = 'best')

def pixel(arcsec):
    return arcsec / pixscale

def arcsec(pixel):
    return pixel * pixscale

axsec = ax1.secondary_xaxis('top', functions=(pixel, arcsec))
axsec.set_xlabel('r (pixels)', size = 'x-large')
ax1.set_xlabel('r (arcseconds)', size = 'x-large')

ax1.set_ylabel('Normalized Flux', size = 'x-large')

ax1.minorticks_on()
axsec.minorticks_on()

ax1.set_title('Clump {:d}'.format(clumps['clump_num'][clump]), y = 0.9)


ind_ap = np.linspace(5, pix_rad-2, 3)
norm = simple_norm(pah_11_k1, vmin = -0.05, vmax = rp_max_val, stretch = 'asinh')
ax2.imshow(pah_11_k1, norm=norm, origin='lower')
rp.apertures[int(ind_ap[0])].plot(ax=ax2, color='C0', lw=2)
rp.apertures[int(ind_ap[1])].plot(ax=ax2, color='C1', lw=2)
rp.apertures[-1].plot(ax=ax2, color='C3', lw=2)
# rp.apertures[25].plot(ax=ax, color='C4', lw=2)
ax2.set_xlim(xycen[0] - pix_rad - 20, xycen[0] + pix_rad + 20)
ax2.set_ylim(xycen[1] - pix_rad - 20, xycen[1] + pix_rad + 20)
ax2.scatter(xycen[0], xycen[1], c = 'r', marker = 'x')
ax2.set_title('PAH 11.3')


norm = simple_norm(F1500W['data'], vmin = -0.05, vmax = rp_max_val_F1500W, stretch = 'asinh')
ax3.imshow(F1500W['data'], norm=norm, origin='lower')
rp.apertures[int(ind_ap[0])].plot(ax=ax3, color='C0', lw=2)
rp.apertures[int(ind_ap[1])].plot(ax=ax3, color='C1', lw=2)
rp.apertures[-1].plot(ax=ax3, color='C3', lw=2)
# rp.apertures[25].plot(ax=ax, color='C4', lw=2)
ax3.set_xlim(xycen[0] - pix_rad - 20, xycen[0] + pix_rad + 20)
ax3.set_ylim(xycen[1] - pix_rad - 20, xycen[1] + pix_rad + 20)
ax3.scatter(xycen[0], xycen[1], c = 'r', marker = 'x')
ax3.set_title('F1500W')

save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/radial_profiles/'
save_name = 'clump_{:d}_radprof_11.3_PAH_unnorm.png'.format(clumps['clump_num'][clump])
# plt.savefig(save_path + save_name, bbox_inches='tight',pad_inches = 0.1, dpi = 300)


# plot just the profile alone
plt.figure(5, figsize=(7, 6))
plt.clf()
ax1 = plt.gca()
ax1.plot(rp.radius * pixscale, rp.profile, zorder = 5, c = '#0072B2', label = '11.3 $\mathrm{\mu}$m PAH')
ax1.axhline(pah11_sigma*(rp_max_norm/rp_max),  c = '#0072B2', lw = 0.5, ls = '--')
ax1.fill_between(rp.radius * pixscale, rp.profile + rp.profile_error, rp.profile - rp.profile_error, zorder = 0.01, color = '#0072B2', alpha = 0.1)

ax1.plot(rp_con.radius * pixscale, rp_con.profile, zorder = 5, c = '#009E73', label = '11.3 $\mathrm{\mu}$m CON')
ax1.fill_between(rp_con.radius * pixscale, rp_con.profile + rp_con.profile_error,  rp_con.profile - rp_con.profile_error, zorder = 0.01, color = '#33B58F', alpha = 0.1)
ax1.axhline(pah11_sigma*(rp_con_max_norm/rp_con_max),  c = '#E69F00', lw = 0.5, ls = '--')

ax1.plot(rp_F1500W.radius * pixscale, rp_F1500W.profile, zorder = 5, c = '#E69F00', label = 'F1500W')
ax1.axhline(F1500W_sigma*(rp_F1500W_max_norm/rp_F1500W_max),  c = '#E69F00', lw = 0.5, ls = '--')
ax1.fill_between(rp_F1500W.radius * pixscale, rp_F1500W.profile + rp_F1500W.profile_error, rp_F1500W.profile - rp_F1500W.profile_error, zorder = 0.01, color = '#F0B84A', alpha = 0.1)

ax1.set_xlim(0.05, edge_radii[-2] * pixscale)
ax1.axvline(beam_fwhm, c = 'k', lw = 0.5, ls = '--', zorder = 0.000001)

ax1.legend(loc = 'best')

def pixel(arcsec):
    return arcsec / pixscale

def arcsec(pixel):
    return pixel * pixscale

axsec = ax1.secondary_xaxis('top', functions=(pixel, arcsec))
axsec.set_xlabel('r (pixels)', size = 'x-large')
ax1.set_xlabel('r (arcseconds)', size = 'x-large')

ax1.set_ylabel('Normalized Flux', size = 'x-large')

ax1.minorticks_on()
axsec.minorticks_on()

# ax1.set_title('Clump {:d}'.format(clumps['clump_num'][clump]), y = 0.9)

save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/radial_profiles/'
save_name = 'clump_{:d}_radprof_11.3_PAH_paper.png'.format(clumps['clump_num'][clump])
plt.savefig(save_path + save_name, bbox_inches='tight',pad_inches = 0.1, dpi = 300)


