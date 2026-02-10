from tweakwcs import JWSTgWCS
from jwst import datamodels
import numpy as np
import glob

import matplotlib.pyplot as plt
from astropy.stats import sigma_clip

RAD2ARCSEC = 3600.0 * np.rad2deg(1.0)

filter = "F560W"
work_dir = '../ic1613/miri/'

twfiles = np.sort(glob.glob(work_dir + f"./{filter}/stage3/*i2d.fits"))

shifts = np.zeros((2, len(twfiles)))
for k, cfile in enumerate(twfiles):
    # extract x,y shifts from the aligned image:
    aligned_model = datamodels.open(cfile)
    matrix = aligned_model.meta.wcs.forward_transform['tp_affine'].matrix.value
    cshift = RAD2ARCSEC * aligned_model.meta.wcs.forward_transform['tp_affine'].translation.value
    shifts[:, k] = cshift
    print(cfile, cshift)

plt.plot(range(len(twfiles)), shifts[0, :], "b-")
plt.plot(range(len(twfiles)), shifts[1, :], "g-")

# determine the aveage values for each tile
for k in range(1):
    k1 = k*10
    k2 = k1 + 10
    avex = np.average(sigma_clip(shifts[0, k1:k2]))
    avey = np.average(sigma_clip(shifts[1, k1:k2]))
    print(k, avex, avey)
    plt.axhline(avex, linestyle="dotted", color="b")
    plt.axhline(avey, linestyle="dotted", color="g")