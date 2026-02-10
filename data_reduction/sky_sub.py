
import copy
import warnings
import numpy as np


from astropy.modeling import models, fitting

from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import fits
from astropy.wcs import WCS

# from astropy.wcs.utils import proj_plane_pixel_scales
from regions import Regions
from astropy.coordinates import SkyCoord

from jwst import datamodels
from jwst.datamodels import dqflags
from tweakwcs import JWSTgWCS

import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, ManualInterval, SqrtStretch

import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, ManualInterval, SqrtStretch


# modified from Karl Gorden's SMACS Field MIRI processing

def plot_stack(istack, vmin, vmax):
    data_2d = istack[:,:,0]

    length = istack.shape[2]


    norm = ImageNormalize(
        data_2d, interval=ManualInterval(vmin=vmin, vmax=vmax), stretch=SqrtStretch()
    )

    fig, axs = plt.subplots(2,2, figsize=(12, 10))

    for i, ax in enumerate(fig.axes):
        im = ax.imshow(istack[:,:,i], origin="lower", norm=norm, cmap=plt.get_cmap("magma_r"))

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, label="Flux", cax = cbar_ax)
    plt.xlabel("Pixel column")
    plt.ylabel("Pixel row")




def make_sky(
    files,
    subfiles=None,
    scalebkg=False,
    exclude_above=None,
    exclude_delta=None,
    ds9regions=None, 
    save_bkgr=False
):
    """
    Make sky background by sigma clipping in image coordinates and subtract it
    from all the input files.

    Parameters
    ----------
    files : strs
       Array of cal files to use to create the sky image
    scalebkg : boolean
       Scale each image by its median to the average value [default=False]
    exclude_above : float
       Exclude data above this value from the sky creation
    exclude_delta : float
       Exclude data above the median bkg + this value from sky creation
    ds9regions : ds9 region file
       Exclude pixels inside ds9 regions from sky creation
    """
    if ds9regions is not None:
        ereg = Regions.read(ds9regions, format="ds9")
        # for creg in ereg:
        #     creg.radius *= 0.5

    istack = None
    for k, cfile in enumerate(files):
        print(f"processing {cfile}")
        cdata = datamodels.open(cfile)
        if istack is None:
            isize = cdata.data.shape
            istack = np.empty((isize[0], isize[1], len(files)))
            istackmed = np.empty((len(files)))
        tdata = cdata.data

        # remove all the non imager data
        # bdata = cdata.dq & dqflags.pixel["DO_NOT_USE"] > 0
        # tdata[bdata] = np.NaN

        if exclude_above is not None:
            tdata[tdata > exclude_above] = np.NaN

        if ds9regions is not None:
            # radeg, decdeg = cdata.meta.wcs([500., 600.], [500., 400.])
            # skycoord = SkyCoord(radeg, decdeg, unit='deg')
            # print(skycoord)
            # get standard WCS info from FITS header
            # with warnings.catch_warnings():
            #     warnings.simplefilter("ignore")
            #     t = fits.open(cfile)
            #     cwcs = WCS(t[1].header)
            # cwcs = cdata.meta.wcs.to_fits()

            fits_header, fits_hdulist = cdata.meta.wcs.to_fits()
            cwcs = WCS(fits_header)  # <-- "astropy" wcs

            pixx = np.arange(isize[1])
            pixy = np.arange(isize[0])
            imagex, imagey = np.meshgrid(pixx, pixy)
            imagera, imagedec = cwcs.wcs_pix2world(imagex, imagey, 0)
            # imagera, imagedec = cwcs.pixel_to_world(imagex, imagey, 0)
            skycoord = SkyCoord(imagera, imagedec, unit="deg")
            for creg in ereg:
                inoutimage = creg.contains(skycoord, cwcs)
                tdata[inoutimage] = np.NaN
            cdata.data = tdata
            cdata.write(cfile.replace("cal.fits", "cal_mask.fits"))
            # fits.writeto("test.fits", inoutimage * 1., overwrite=True)

        istackmed[k] = np.nanmedian(tdata)
        print(f"median sky = {istackmed[k]}")

        if exclude_delta is not None:
            tdata[tdata > istackmed[k] + exclude_delta] = np.NaN

        istack[:, :, k] = tdata

    # adjust the levels to the median
    # allows for data taken at different times with different backgrounds
    medsky = np.mean(istackmed)
    if scalebkg:
        for k in range(len(files)):
            istack[:, :, k] += medsky - istackmed[k]
            print(k, np.nanmedian(istack[:, :, k]))
    else:
        print("Not scaling individual images to median bkg")

    skyflat_mean, skyflat_median, skyflat_std = sigma_clipped_stats(
        istack, sigma_lower=3, sigma_upper=1, axis=2
    )



    # subtract the sky properly adjusted from the data
    if subfiles is None:
        subfiles = files
    for k, cfile in enumerate(subfiles):
        cdata = datamodels.open(cfile)
        cdata.data -= skyflat_mean
        if scalebkg:
            print(cfile, medsky - istackmed[k])
            cdata.data += medsky - istackmed[k]
        else:
            print(cfile)
        ndata = np.isnan(cdata.data)
        cdata.data[ndata] = 0.0
        cdata.dq[ndata] = cdata.dq[ndata] & dqflags.pixel["DO_NOT_USE"]
        cdata.write(cfile.replace("_cal.fits", "_skysub_cal.fits"))

    if save_bkgr:
        fits.writeto(cfile.replace("_cal.fits", "_sky_median.fits"), skyflat_median, overwrite = True)
        fits.writeto(cfile.replace("_cal.fits", "_sky_mean.fits"), skyflat_mean, overwrite = True)
        fits.writeto(cfile.replace("_cal.fits", "_sky_std.fits"), skyflat_std, overwrite = True)


    return skyflat_mean, skyflat_median, skyflat_std

