

#Modify the path to a directory on your machine
import os
os.environ["CRDS_PATH"] = "/user/etarantino/crds_cache/"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
os.environ["CRDS_CONTEXT"] = "jwst_1256.pmap"

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

# custom written code
from correct1f_gen import run_1fcor

from jhat import jwst_photclass,st_wcs_align, st_wcs_align_batch
wcs_align = st_wcs_align()


# set up the observation info and where the uncal files are
inpath = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/uncal_aug24/mastDownload/JWST/' 
prog = 'jw02391007001'
asn_text = 'jw02391-o007_t002_nircam_v2'

# where the pipeline files will go (.cal, .rate, mosaics, etc)
outpath = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/nircam_jan25_cleanrun/'

# reference catalog for step 3 of the pipeline
ref_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/'
ref_name = 'HST16104_SEXTANS-A_catalog_5k.ecsv'

from jhat import jwst_photclass,st_wcs_align, st_wcs_align_batch
wcs_align = st_wcs_align()

# use both or only one module
modb = False

# flags for code 
run_detector1 = False
run_image2 = False
run_image3 = True

# stage 1
keep_group = True
snowballs = True
cpufraction = '1'

# stage 2
skip_resample = True
cor_1f = True

# stage 3
abs_cat = True
sky_match = False
run_jhat = False

# background fitting
bkgr_sub = False
scalebkg = True
savebkg = True


# generate observation dictionary containing the dithers, detectors, filters, etc

if modb:
    nircam = {'obs': ['03101', '05101', '07101'],
          'filt': [ 'f115w','f150w', 'f200w', 'f300m', 'f335m','f360m'],
          'det': ['nrcblong', 'nrcb1', 'nrcb2', 'nrcb3', 'nrcb4'],
          'exp': ['00001', '00002', '00003', '00004']
          }

else: 

    nircam = {'obs': ['03101', '05101', '07101'],
              'filt': ['f360m', 'f115w','f150w', 'f200w', 'f300m', 'f335m'],
              'det': ['nrcalong', 'nrcblong', 'nrca1', 'nrca2', 'nrca3', 'nrca4', 'nrcb1', 'nrcb2', 'nrcb3', 'nrcb4'],
              'exp': ['00001', '00002', '00003', '00004']
              }



short_filts = {'f115w', 'f150w', 'f200w'}
long_filts = {'f300m', 'f335m', 'f360m'}
short_det = ['nrca1', 'nrca2', 'nrca3', 'nrca4', 'nrcb1', 'nrcb2', 'nrcb3', 'nrcb4']
long_det = ['nrcalong', 'nrcblong']
short_det_b = ['nrcb1', 'nrcb2', 'nrcb3', 'nrcb4']
long_det_b = ['nrcblong']


# loop through each observation, then detector, then dither 
for obs in nircam['obs']:
    for det in nircam['det']:
        for exp in nircam['exp']:

            # set up files and paths
            name = prog + '_' + obs + '_' + exp + '_' + det
            uncal_path = inpath + name + '/'
            uncal_file = name + '_uncal.fits'

            # read the file and find the corresponding filter
            uncal_head = fits.open(uncal_path + uncal_file)[0].header
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

            print('WORKING ON ' + name + ' FILTER ' + filt)

            ###############################
            ########### step 1 ############
            ###############################

            if run_detector1: 
                # setup dictionary for parameters for stage 1
                det1_dict = {}

                # allows computations for a saturated ramp with one good sample
                det1_dict['ramp_fit'] = {'suppress_one_group': not keep_group, 'maximum_cores': cpufraction}


                # can help remove snowballs, expands the number of pixels that will be flagged around large cosmic ray events
                det1_dict['jump'] = {'expand_large_events': snowballs, 'save_results': True, 'maximum_cores' : cpufraction}

                # runs the pipeline with all the settings set above
                calwebb_detector1.Detector1Pipeline.call(uncal_path + uncal_file, steps=det1_dict, output_dir=workdir, save_results=True) 

            ################################
            ############ step 2 ############
            ################################

            if run_image2:
                rate_file = name + '_rate.fits'

                # dictionary for parameters for stage 2
                im2_dict = {}

                # skip the resample step because it will be done in stage 3
                if skip_resample:
                    im2_dict["resample"] = {"skip":True}

                # runs the pipeline with all the settings set above
                calwebb_image2.Image2Pipeline.call(workdir + rate_file, steps=im2_dict, output_dir=workdir, save_results=True)



################################
############ step 3 ############
################################

if run_image3:
    for filt in nircam['filt']:

        # the folders use uppercase letters for filters
        up_filt = filt.upper()

        # pipeline working directory 
        workdir = outpath + up_filt + '/'


        # resampling factors given from Martha
        if filt in short_filts:
            pix_scale = 0.021

            if modb:
                det_array = short_det_b
            else:
                det_array = short_det

            short = True
        elif filt in long_filts:
            pix_scale = 0.042
            if modb:
                det_array = long_det_b
            else:
                det_array = long_det

            short = False
        else:
            print('Wrong Filters!!!')

        # run 1f correction if needed
        # identify all of the *cal files
        if cor_1f:
            # needed or else it will also grab the 1fcorrected files 
            if short:
                nircam_files = (glob.glob(workdir + f"*nrc??_cal.fits"))
            else:
                nircam_files = (glob.glob(workdir + f"*long_cal.fits"))

            run_1fcor(nircam_files)

            nircam_files = glob.glob(workdir + f"*1fcor_cal.fits")
            nircam_asn_name = workdir + asn_text + f'{filt}_1fcor_orig'   

        else:
            # needed or else it will also grab the 1fcorrected files 
            if short:
                nircam_files = (glob.glob(workdir + f"*nrc??_cal.fits"))
            else:
                nircam_files = (glob.glob(workdir + f"*long_cal.fits"))

            nircam_asn_name = workdir + asn_text + f'{filt}'  



        if bkgr_sub:
            for det in det_array:
                print('Working on ', det)
                cal_files = glob.glob(workdir + f"*{det}_1fcor_cal.fits")
                sky_bkg = sky_sub(cal_files, scalebkg = scalebkg, savebkg = savebkg)

            # get background files and make asn files
            if modb:
                if short:
                    nircam_files = (glob.glob(workdir + f"*nrcb?_1fcor_bkgrsub_cal.fits"))
                else:
                    nircam_files = (glob.glob(workdir + f"*blong_1fcor_bkgrsub_cal.fits"))

                nircam_asn_name = workdir + f'{prog}_nircam_{filt}_asn_modb_1fcor_bkgrsub'   
            
            else: 
                nircam_files = (glob.glob(workdir + f"*1fcor_bkgrsub_cal.fits"))

                nircam_asn_name = workdir + f'{prog}_nircam_{filt}_asn_1fcor_bkgrsub'   

        else:
            if modb:
                if short:
                    nircam_files = (glob.glob(workdir + f"*nrcb?_1fcor_cal.fits"))
                else:
                    nircam_files = (glob.glob(workdir + f"*blong_1fcor_cal.fits"))

                nircam_asn_name = workdir + f'{prog}_nircam_{filt}_asn_modb_1fcor'   
            
            else: 
                if short:
                    nircam_files = (glob.glob(workdir + f"*nrc??_1fcor_cal.fits"))
                else:
                    nircam_files = (glob.glob(workdir + f"*long_1fcor_cal.fits"))

                nircam_asn_name = workdir + f'{prog}_nircam_{filt}_asn_1fcor'   

        # code for jhat alignment
        if run_jhat:
            if bkgr_sub:
                if modb:
                    if short:
                        nircam_files = (glob.glob(workdir + f"*nrcb?_1fcor_bkgrsub_cal.fits"))
                    else:
                        nircam_files = (glob.glob(workdir + f"*blong_1fcor_bkgrsub_cal.fits"))

                    nircam_asn_name = workdir + f'{prog}_nircam_{filt}_asn_modb_1fcor_bkgrsub_jhat'   
                
                else: 
                    nircam_files = (glob.glob(workdir + f"*1fcor_bkgrsub_cal.fits"))

                    nircam_asn_name = workdir + f'{prog}_nircam_{filt}_asn_1fcor_bkgrsub_jhat'   

            savedir = f'../sexa/nircam_jan25_cleanrun/{up_filt}/'


            for align_im in nircam_files:

                d2d_max = 0.3
                # xshift = 6
                # yshift = -1
                xshift = 0
                yshift = 0

                wcs_align.run_all(align_im,
                      telescope='jwst',
                      outsubdir=savedir,
                      overwrite=True,
                      d2d_max=d2d_max,
                      showplots=0,
                      saveplots=2,
                      refcatname=ref_path + ref_name,
                      refcat_racol='RA',
                      refcat_deccol='DEC',
                      refcat_magcol='mag',
                      refcat_magerr_col = 'err',
                      histocut_order='dxdy',
                      sharpness_lim=(0.3,0.9),
                      roundness1_lim=(-0.7, 0.7),
                      SNR_min= 3,
                      iterate_with_xyshifts = True, 
                      xshift = xshift,
                      yshift = yshift)

            jhat_files = glob.glob(workdir + f"*_jhat.fits")
            nircam_files = jhat_files            


        print(nircam_files)
        nircam_asn = asn_from_list.asn_from_list(nircam_files, rule=DMS_Level3_Base, product_name=nircam_asn_name )

        nircam_asn_file = f'{nircam_asn_name}_v2.json'
        with open(nircam_asn_file, 'w') as outfile:
            name, serialized = nircam_asn.dump(format='json')
            outfile.write(serialized)

        # setup the directory of step paramters
        im3_dict = {}

        # set pixel scale to recommended values
        im3_dict['resample'] = {'pixel_scale': pix_scale}

        # debugging the outlier detection step
        # im3_dict['outlier_detection'] = {'save_intermediate_results': True, "save_results": True, "snr": "3.0 2.0", "grow": 3}

        # use the supplied catalog as the astrometry reference catalog
        if abs_cat:

            im3_dict["tweakreg"] = {"snr_threshold": 5, "expand_refcat": True, "tolerance": 0.1, "save_results": True, "abs_refcat": ref_path + ref_name, "save_catalogs": True, 'save_abs_catalog':True, 'output_dir': workdir}

        elif run_jhat:
            im3_dict["tweakreg"] = {'skip': True}


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










