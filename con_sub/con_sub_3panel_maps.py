#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates plots of the continuum subtracted PAH features to compare to one another

@author: etarantino
"""
import aplpy
import matplotlib.pyplot as plt
import astropy.units as u

plt.ion()

con_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/con_sub_images/'
reproj_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/psf-match/reproject/'

pah_3_file = 'SextansA_F335M_sub_match_s23_pah.fits'
pah_7_file = 'SextansA_F770W_sub_match_pah.fits'
pah_11_file = 'SextansA_F1130W_sub_match_pah.fits'

# pah_3_file = 'SextansA_F335M_sub_match_s23_pah.fits'
# pah_7_file = 'SextansA_F335M_sub_match_pah.fits'
# pah_11_file = 'SextansA_F335M_method_diff.fits'

# filt = 'F335M'
# file1 = f'{filt}_reproject_to_F1500W.fits'
# file2 = f'SextansA_{filt}_sub_match_con.fits'
# file3 = f'SextansA_{filt}_sub_match_pah.fits'

plt.figure(1, figsize=(15, 7))
plt.clf()
fig = plt.figure(figsize=(15, 7), num = 1)

f1 = aplpy.FITSFigure(con_path + pah_3_file, figure=fig, subplot=[0.1,0.1,0.25,0.8],  north = True)
f2 = aplpy.FITSFigure(con_path + pah_7_file, figure=fig, subplot=[0.38,0.1,0.25,0.8],  north = True)
f3 = aplpy.FITSFigure(con_path + pah_11_file, figure=fig, subplot=[0.66,0.1,0.25,0.8],  north = True)

# f1 = aplpy.FITSFigure(reproj_path + file1, figure=fig, subplot=[0.1,0.1,0.25,0.8],  north = True)
# f2 = aplpy.FITSFigure(con_path + file2, figure=fig, subplot=[0.38,0.1,0.25,0.8],  north = True)
# f3 = aplpy.FITSFigure(con_path + file3, figure=fig, subplot=[0.66,0.1,0.25,0.8],  north = True)

vmin = -0.01
vmax = 0.1

f2.axis_labels.hide_y()
f2.tick_labels.hide_y()

f3.axis_labels.hide_y()
f3.tick_labels.hide_y()

cen_x = 152.7780119
cen_y = -4.7041120

f1.recenter(cen_x, cen_y, width=0.005,height=0.005)
f2.recenter(cen_x, cen_y, width=0.005,height=0.005)
f3.recenter(cen_x, cen_y, width=0.005,height=0.005)


fig.canvas.draw()

f1.add_scalebar(1 * u.arcsec, corner = 'bottom left', color = 'white')
f1.scalebar.set_label("1'' = 7 pc")

f2.add_scalebar(1 * u.arcsec, corner = 'bottom left', color = 'white')
f2.scalebar.set_label("1'' = 7 pc")

f3.add_scalebar(1 * u.arcsec, corner = 'bottom left', color = 'white')
f3.scalebar.set_label("1'' = 7 pc")

# f1.add_beam()
# f1.beam.show(major = 0.5*u.arcsec, minor = 0.5*u.arcsec, corner  = 'bottom right')

f1.show_colorscale(cmap = 'magma', vmin = vmin, vmax = vmax)
f2.show_colorscale(cmap = 'magma', vmin = vmin, vmax = vmax)
f3.show_colorscale(cmap = 'magma', vmin = vmin, vmax = vmax)

f1.add_colorbar()
f1.colorbar.hide()

f2.add_colorbar()
f2.colorbar.hide()

f3.add_colorbar(location = 'right')
f3.colorbar.set_axis_label_text('MJy/sr')
f3.colorbar.set_pad(0.1)

save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/plots/'
# save_name = f'{filt}_PAH_con_sub_compare.pdf'
save_name = 'SextansA_consub_all_filt_compare.pdf'

plt.savefig(save_path + save_name, overwrite = True)
