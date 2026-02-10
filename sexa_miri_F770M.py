import os
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
from astropy.utils.data import download_file
from astropy.visualization import ImageNormalize, ManualInterval, LogStretch

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
filter = 'F560W'
basename = 'jw02391005001'
work_dir = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/miri'

# variables for stage 3
scalebkg = True
savebkg = True
pix_scale = 0.11
rotation = None

# HST reference catalog
ref_path = '/astro/dust_kg/etarantino/JWST_PAHS_2391/sexa/'
ref_name = 'HST16104_SEXTANS-A_catalog_cut.fits'

# # shifts determine from previous runs with treakreg or with other filters
# #   analysis using Analyze_teakreg_shifts.ipynb
# #   run this notebook after running shortest wavelength dataset with tweakreg=True and align_to_gaia=True
# #     and setting all the tile?_shifts to [0.0, 0.0]
# #   use the resulting shifts here and set tweakreg=False and align_to_gaia=False
# #tile1_shifts = [0.0, 0.0]
# #tweakreg = True
# #align_to_gaia = True
# tile1_shifts = [0.006003914950907428, -1.6397513463822861]
# tweakreg = False
# align_to_gaia = False

# display range for uncalibrated images
drange_raw = [17., 50.]
# display range for calibrated images
drange_cal = [5., 40.]
dmap = "afmhot"  # same as ds9 bb

run_detector1 = False
run_image2 = False
run_image3 = True
make_bkg = True

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
    miri_detector1(miri_uncal_files, output_dir)

# prep stage 2
miri_rate_files = glob.glob(work_dir + f"/{filter}/stage1/*rate.fits")
output_dir = work_dir + f'/{filter}/stage2'
print(miri_rate_files)

# run stage 2
if run_image2:
    miri_image2(miri_rate_files, output_dir)

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
    miri_image3(miri_asn_file, output_dir, tweakreg=True, 
                pixel_scale=pix_scale, absolute_cat = ref_path + ref_name)



