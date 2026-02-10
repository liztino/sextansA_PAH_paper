from image1overf_gen import sub1fimaging #See note.
from jwst import datamodels
from astropy.io import fits
import glob
import numpy as np 
from jwst.datamodels import dqflags

# nircam_cal_files=glob.glob('*_cal.fits')
# print(nircam_cal_files)

def run_1fcor(nircam_cal_files):
    for i in range(0, len(nircam_cal_files)):
        cal2file = nircam_cal_files[i]
        cal21overffile = cal2file.replace('_cal.fits', '_1fcor_cal.fits')
        print ('Running 1/f correction on {} to produce {}'.format(cal2file,cal21overffile))
        with fits.open(cal2file) as cal2hdulist:
            sigma_bgmask=3.0
            sigma_1fmask=2.0
            splitamps=False
            correcteddata = sub1fimaging(cal2hdulist,sigma_bgmask,sigma_1fmask,splitamps)

            # remove the reference pixel part of the array 
            ref_pix = cal2hdulist['DQ'].data & dqflags.pixel['REFERENCE_PIXEL'] == dqflags.pixel['REFERENCE_PIXEL']
            ref_xx, ref_yy = np.where(ref_pix==False)

            xind1 = ref_xx[0]
            xind2 = ref_xx[-1]+1
            yind1 = ref_yy[0]
            yind2 = ref_yy[-1]+1

            cal2hdulist['SCI'].data[xind1:xind2, yind1:yind2] = correcteddata

            cal2hdulist.writeto(cal21overffile, overwrite=True)

    print('Done!')
