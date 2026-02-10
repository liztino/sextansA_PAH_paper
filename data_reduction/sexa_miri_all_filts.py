import os, sys
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import glob
import copy
os.environ["CRDS_PATH"] = "/user/etarantino/crds_cache/"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

# Packages that allow us to get information about objects:
import asdf
import copy
import shutil

# Numpy library:
import numpy as np

# To read association file
import json
from jwst.associations import load_asn
from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

# Astropy tools:
from astropy.io import fits
# from astropy.utils.data import download_file
# from astropy.visualization import ImageNormalize, ManualInterval, LogStretch

import matplotlib.pyplot as plt
import matplotlib as mpl

# List of possible data quality flags
from jwst.datamodels import dqflags
from jwst import datamodels

# custom helper routines from Karl
from jwst_helpers import show_image, overlay_catalog
from miri_helpers import miri_detector1, miri_image2, miri_image3
from miri_clean import shift_cal_wcs, make_sky

# The entire pipeline
from jwst.pipeline import calwebb_detector1
from jwst.pipeline import calwebb_image2
from jwst.pipeline import calwebb_image3


# variables controlling filter info
# filt_list = ['F560W']
filt_list = ['F560W', 'F770W', 'F1000W', 'F1130W', 'F1500W']

work_dir = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/miri_aug24'

# variables for stage 3
scalebkg = True
savebkg = True
pix_scale = 0.11
rotation = None

# shifts determine from previous runs with treakreg or with other filters
#   analysis using Analyze_teakreg_shifts.ipynb
#   run this notebook after running shortest wavelength dataset with tweakreg=True and align_to_gaia=True
#     and setting all the tile?_shifts to [0.0, 0.0]
#   use the resulting shifts here and set tweakreg=False and align_to_gaia=False
shifts = [-0.23363470808898834, 0.011315765654765218]

# HST reference catalog
ref_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/'
ref_name = 'HST16104_SEXTANS-A_catalog_cut.fits'

# reference catalog for F1130W
ref_path_F1130W = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/'
ref_name_F1130W = 'JWST2391_F1000W_catalog.ecsv'

run_detector1 = False
run_image2 = False
run_image3 = True

cpufraction = '1'
make_bkg = True

def main():

    for filter in filt_list:

        if filter == 'F1130W':
            align_short = False

            ref_path = ref_path_F1130W
            ref_name = ref_name_F1130W 

        else: 
            align_short = True


        # Since the F1130W filter was taken months after the others, we can't use the GS trick to align
        # instead, we use the F1000W catalog 
        if filter == 'F1130W':
            align_short = False
            # reference catalog for F1000W filter
            ref_path = ref_path_F1130W
            ref_name = ref_name_F1130W
        elif filter == 'F560W':
            align_short = False
            # HST reference catalog
            ref_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/'
            ref_name = 'HST16104_SEXTANS-A_catalog_cut.fits'

        else:
            align_short = True

        # create stage directories if they do not exist
        for k in range(4):
            cpath = f"/{filter}/stage{k}"
            if not os.path.exists(work_dir + cpath):
                os.makedirs(work_dir + cpath)

        # prep stage 1
        miri_uncal_files = glob.glob(work_dir + f"/{filter}/stage0/*uncal.fits")
        output_dir = work_dir + f'/{filter}/stage1'
        print(miri_uncal_files)

        # run stage 1
        if run_detector1:
            print('running detector1..."')      
            miri_detector1(miri_uncal_files, output_dir, cpufraction=cpufraction, save_jump_info=True)

        # prep stage 2
        miri_rate_files = glob.glob(work_dir + f"/{filter}/stage1/*rate.fits")
        output_dir = work_dir + f'/{filter}/stage2'
        print(miri_rate_files)

        # run stage 2
        if run_image2:
            miri_image2(miri_rate_files, output_dir)

        if align_short:
            miri_cal_files = glob.glob(work_dir + f"/{filter}/stage2/*_mirimage_cal.fits")

            for cfile in miri_cal_files:
                shift_cal_wcs(cfile, shifts)

            # make sky and subtract it
            if make_bkg:
                miri_cal_files = glob.glob(work_dir + f"/{filter}/stage2/*_mirimage_fixed_wcs_cal.fits")
                sky_bkg = make_sky(miri_cal_files, scalebkg = scalebkg, savebkg = savebkg)

                miri_files = glob.glob(work_dir + f"/{filter}/stage2/*_fixed_wcs_skysub_cal.fits")
                output_dir =  work_dir + f"/{filter}/stage3/"
                print(miri_files)

                miri_asn_name = output_dir + f'/miri_{filter}_stage3_asn_skysub'   

                miri_asn = asn_from_list.asn_from_list(miri_files, rule=DMS_Level3_Base, product_name=miri_asn_name)

                miri_asn_file = f'{miri_asn_name}.json'
                with open(miri_asn_file, 'w') as outfile:
                    name, serialized = miri_asn.dump(format='json')
                    outfile.write(serialized)

            if run_image3:
                output_dir =  work_dir + f"/{filter}/stage3/"
                miri_asn_name = output_dir + f'/miri_{filter}_stage3_asn_skysub' 
                miri_asn_file = f'{miri_asn_name}.json'
                miri_image3(miri_asn_file, output_dir, tweakreg=False, 
                        pixel_scale=pix_scale)
        else:

            # make sky and subtract it
            if make_bkg:
                miri_cal_files = glob.glob(work_dir + f"/{filter}/stage2/*_mirimage_cal.fits")
                sky_bkg = make_sky(miri_cal_files, scalebkg = scalebkg, savebkg = savebkg)

                miri_files = glob.glob(work_dir + f"/{filter}/stage2/*_skysub_cal.fits")
                output_dir =  work_dir + f"/{filter}/stage3/"
                print(miri_files)

                if rotation is not None:
                    miri_asn_name = output_dir + f'miri_{filter}_stage3_nirproj_asn_skysub'
                else: 
                    miri_asn_name = output_dir + f'/miri_{filter}_stage3_asn_skysub'   

                miri_asn = asn_from_list.asn_from_list(miri_files, rule=DMS_Level3_Base, product_name=miri_asn_name)

                miri_asn_file = f'{miri_asn_name}.json'
                with open(miri_asn_file, 'w') as outfile:
                    name, serialized = miri_asn.dump(format='json')
                    outfile.write(serialized)

            if run_image3:
                output_dir =  work_dir + f"/{filter}/stage3/"
                miri_asn_name = output_dir + f'/miri_{filter}_stage3_asn_skysub' 
                miri_asn_file = f'{miri_asn_name}.json'
                miri_image3(miri_asn_file, output_dir, tweakreg=True, 
                            pixel_scale=pix_scale, absolute_cat = ref_path + ref_name)


if __name__ == '__main__':
    sys.exit(main())



