'''

Script to query MAST for the uncal files from a desired program

'''

from astroquery.mast import Observations
from download_step import DownloadStep
import glob
import os
import numpy as np
import shutil
from astropy.io import ascii, fits

NIRCam = False
MIRI = True

### MAST query NIRCam ###
if NIRCam: 
	prop_id = '2391'
	target = 'Sextans A'
	download_dir = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/uncal_aug24'
	instrument_name = 'NIRCAM/IMAGE'
	# command to get only uncal files
	calib_level = [0, 1]
	extension = 'fits'
	login = True
	api_key = ' bce4c11ca6cd49e59d6a88bba79a77b5'
	product_type = None

	test_down = DownloadStep(target = target, prop_id=prop_id, download_dir = download_dir, 
		instrument_name = instrument_name, calib_level = calib_level, extension=extension, 
		login = True, api_key = api_key, product_type=product_type)

	test_down.do_step()

### MAST query MIRI ###
if MIRI: 
	prop_id = '2391'
	target = 'Sextans A'
	download_dir = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/uncal_miri_aug24'
	instrument_name = 'MIRI/IMAGE'
	# command to get only uncal files
	calib_level = [0, 1]
	extension = 'fits'
	login = True
	api_key = ' bce4c11ca6cd49e59d6a88bba79a77b5'
	product_type = None

	test_down = DownloadStep(target = target, prop_id=prop_id, download_dir = download_dir, 
		instrument_name = instrument_name, calib_level = calib_level, extension=extension, 
		login = True, api_key = api_key, product_type=product_type)

	test_down.do_step()

	main_dir = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/miri_aug24/'

	files = glob.glob(download_dir + '/mastDownload/JWST/*/*uncal.fits')

	for uncal in files:

		hdu = fits.open(uncal)
		head = hdu[0].header
		filt = head['FILTER']

		workdir = main_dir + filt 

		print(workdir)
		if not os.path.isdir(workdir):
			os.mkdir(workdir)

		workdir = workdir + '/stage0/'

		if not os.path.isdir(workdir):
			os.mkdir(workdir)

		shutil.move(uncal,workdir)

