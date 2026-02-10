

#Modify the path to a directory on your machine
import os
os.environ["CRDS_PATH"] = "/user/etarantino/crds_cache/"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# Packages that allow us to get information about objects:
import asdf
import copy
import shutil
import glob


# Numpy library:
import numpy as np

# To read association file
import json
from jwst.associations import load_asn

# Astropy tools:
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.visualization import ImageNormalize, ManualInterval, LogStretch

import matplotlib.pyplot as plt
import matplotlib as mpl

# List of possible data quality flags
from jwst.datamodels import dqflags

# The entire pipeline
from jwst.pipeline import calwebb_detector1
from jwst.pipeline import calwebb_image2
from jwst.pipeline import calwebb_image3

from NIRCam_scripts import sky_sub

# To read association file
import json
from jwst.associations import load_asn
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base


# set up the observation info and where the uncal files are
inpath = '/astro/dust_kg/etarantino/JWST_PAHS_2391/ic1613/uncal/' 
prog = 'jw02391-ic1613'

# where the pipeline files will go (.cal, .rate, mosaics, etc)
outpath = '/astro/dust_kg/etarantino/JWST_PAHS_2391/ic1613/nircam_aug/'

# HST reference catalog
ref_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/ic1613/'
ref_name = 'HST16513_IC-1613_catalog_cut.fits'

# flags for code 
run_detector1 = True
run_image2 = True
run_image3 = True

# stage 1
keep_group = False
snowballs = True

# stage 2
skip_resample = True

# stage 3
abs_cat = True
sky_match = True

# background fitting
bkgr_sub = False
scalebkg = False
savebkg = False

filt_list = ['F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M']
short_filts = {'F115W', 'F150W', 'F200W',}
long_filts = {'F300M', 'F335M', 'F360M'}

# grab uncal files 
uncal_files = glob.glob(inpath + f'*uncal.fits')

for uncal_file in uncal_files:

    uncal_head = fits.open(uncal_file)[0].header
    filt = uncal_head['FILTER']
    pupil = uncal_head['PUPIL']

    if pupil == 'CLEAR':
        filt = filt
    else:
        filt = filt + '-' + pupil

    workdir = outpath  + filt + '/'

    # make directory for filter if it doesn't exist
    if not os.path.isdir(workdir): 
        os.mkdir(workdir)

    print('WORKING ON ' + uncal_file + ' FILTER ' + filt)

    ###############################
    ########### step 1 ############
    ###############################

    if run_detector1: 
        # setup dictionary for parameters for stage 1
        det1_dict = {}

        # allows computations for a saturated ramp with one good sample
        if keep_group:
            det1_dict['ramp_fit'] = {'suppress_one_group': False}

        # can help remove snowballs, expands the number of pixels that will be flagged around large cosmic ray events
        if snowballs:
            det1_dict['jump'] = {'expand_large_events': True}

        det1_dict['jump'] = {'save_results': True}

        # runs the pipeline with all the settings set above
        calwebb_detector1.Detector1Pipeline.call(uncal_file, steps=det1_dict, output_dir=workdir, save_results=True) 

    ################################
    ############ step 2 ############
    ################################

if run_image2:

    for filt in filt_list:

        workdir = outpath + filt + '/'

        rate_files = glob.glob(workdir + f'*rate.fits')

        for rate_file in rate_files:

            # dictionary for parameters for stage 2
            im2_dict = {}

            # skip the resample step because it will be done in stage 3
            if skip_resample:
                im2_dict["resample"] = {"skip":True}

            # runs the pipeline with all the settings set above
            calwebb_image2.Image2Pipeline.call(rate_file, steps=im2_dict, output_dir=workdir, save_results=True)



################################
############ step 3 ############
################################

if run_image3:
    for filt in filt_list:

        # the folders use uppercase letters for filters
        up_filt = filt.upper()

        # pipeline working directory 
        workdir = outpath + up_filt + '/'

        # resampling factors given from Martha
        if filt in short_filts:
            pix_scale = 0.021
            # det_array = short_det
            short = True
        elif filt in long_filts:
            pix_scale = 0.042
            # det_array = long_det
            short = False
        else:
            print('Wrong Filters!!!')

        if bkgr_sub:
            # for det in det_array:
            #     print('Working on ', det)
            #     cal_files = glob.glob(workdir + f"*{det}_cal.fits")
            #     sky_bkg = sky_sub(cal_files, scalebkg = scalebkg, savebkg = savebkg)

            nircam_files = glob.glob(workdir +  "*_modsub_cal.fits")
            # nircam_asn_name = workdir + f'{prog}_nircam_{filt}_asn_modsub'   

        else:
            if short:
                nircam_files = (glob.glob(workdir + f"*n?c??_cal.fits"))

            else:
                nircam_files = (glob.glob(workdir + f"*long_cal.fits"))
                
            nircam_asn_name = workdir + f'{prog}_nircam_{filt}'   

        nircam_asn = asn_from_list.asn_from_list(nircam_files, rule=DMS_Level3_Base, product_name=nircam_asn_name)

        nircam_asn_file = f'{nircam_asn_name}.json'
        with open(nircam_asn_file, 'w') as outfile:
            name, serialized = nircam_asn.dump(format='json')
            outfile.write(serialized)

        # # get list of all files in directory
        # files = os.listdir(workdir)

        # # find file that has "asn" keyword
        # dir_text = [k for k in files if 'asn.json' in k]
        # if len(dir_text) == 1:
        #     asn_file = dir_text[0]
        # else:
        #     print('TOO MANY ASSOCIATION FILES, STOPPING')
        #     break

        # setup the directory of step paramters
        im3_dict = {}

        # set pixel scale to recommended values
        im3_dict['resample'] = {'pixel_scale': pix_scale}

        # debugging the outlier detection step
        # im3_dict['outlier_detection'] = {'save_intermediate_results': True, "save_results": True, "snr": "3.0 2.0", "grow": 3}

        # use the supplied catalog as the astrometry reference catalog
        if abs_cat:
            im3_dict["tweakreg"] = {"save_results": True, "abs_refcat": ref_path + ref_name, "save_catalogs": True, 'save_abs_catalog':True, 'output_dir': workdir}

        # matching the background in pipeline step, don't use if background is turned on
        if not sky_match:
            im3_dict["skymatch"] = {"skip": True}

        calwebb_image3.Image3Pipeline.call(nircam_asn_file, steps=im3_dict, output_dir=workdir, save_results=True)

        # # copy the mosaic into a separate folder
        # mosaic_path = outpath + '/mosaic/'

        # if not os.path.isdir(mosaic_path): 
        #     os.mkdir(mosaic_path)

        # mosaic_file = asn_text +  '-' + filt + '_i2d.fits'
        # shutil.copyfile(workdir + mosaic_file, mosaic_path + mosaic_file)










