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

path = '/Users/etarantino/Documents/JWST_DATA_PAHS/data/'


plt.figure(1, figsize=(15, 7))
plt.clf()
fig = plt.figure(figsize=(15, 7), num = 1)

file1 = 'sextansa_f33.fits'
file2 = 'sextansa_f33_cont.fits'

# f1 = aplpy.FITSFigure(con_path + pah_3_file, figure=fig, subplot=[0.1,0.1,0.25,0.8],  north = True)
# f2 = aplpy.FITSFigure(con_path + pah_7_file, figure=fig, subplot=[0.38,0.1,0.25,0.8],  north = True)
# f3 = aplpy.FITSFigure(con_path + pah_11_file, figure=fig, subplot=[0.66,0.1,0.25,0.8],  north = True)

f1 = aplpy.FITSFigure(path + file1, figure=fig, subplot=[0.1,0.1,0.4,0.8],  north = True)
f2 = aplpy.FITSFigure(path + file2, figure=fig, subplot=[0.5,0.1,0.4,0.8],  north = True)

vmin = -0.01
vmax = 0.1

f2.axis_labels.hide_y()
f2.tick_labels.hide_y()


f1.recenter(152.7770516, -4.7067905, width = 0.02, height = 0.02)
f2.recenter(152.7770516, -4.7067905, width = 0.02, height = 0.02)



fig.canvas.draw()

# f1.add_scalebar(1 * u.arcsec, corner = 'bottom left', color = 'white')
# f1.scalebar.set_label("1'' = 7 pc")

# f2.add_scalebar(1 * u.arcsec, corner = 'bottom left', color = 'white')
# f2.scalebar.set_label("1'' = 7 pc")

# f3.add_scalebar(1 * u.arcsec, corner = 'bottom left', color = 'white')
# f3.scalebar.set_label("1'' = 7 pc")

# f1.add_beam()
# f1.beam.show(major = 0.5*u.arcsec, minor = 0.5*u.arcsec, corner  = 'bottom right')

f1.show_colorscale(cmap = 'magma_r', vmin = vmin, vmax = vmax)
f2.show_colorscale(cmap = 'magma_r', vmin = vmin, vmax = vmax)

f1.add_colorbar()
f1.colorbar.hide()

f2.add_colorbar()
f2.colorbar.set_pad(0.1)

# f2.colorbar.set_box([0.9, 0.1, 0.03, 0.8],box_orientation='vertical')

f2.colorbar.set_axis_label_text('MJy sr$^{-1}')
# f3.colorbar.set_pad(0.04)
# f3.colorbar.set_fraction(0.046)

save_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/plots/'
save_name = 'F335M_S23_consub_orig_cosmicdust.pdf'
plt.savefig(save_path + save_name, overwrite = True)
