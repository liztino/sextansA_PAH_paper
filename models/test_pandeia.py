#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Trying to see if the Pandeia engine makes it easy to do synthetic photometry

@author: etarantino
"""
import numpy as np 
import matplotlib.pyplot as plt
import os
from pandeia.engine.instrument_factory import InstrumentFactory\

# # The following section is only needed if the PYSYN_CDBS environment variable is not already set.
# # The PYSYN_CDBS environment variable should point to the path of the CDBS data files
# import os
# location_of_cdbs = "/path/to/cdbs/files"
# os.environ['PYSYN_CDBS'] = location_of_cdbs
# # End section
os.environ['PYSYN_CDBS']
 
# # The following section is only needed if the pandeia_refdata environment variable is not already set
# # The pandeia_refdata environment variable should point to the path of the pandeia reference data
# import os
# location_of_pandeia_refdata = "/path/to/pandeia/refdata"
# os.environ['pandeia_refdata'] = location_of_pandeia_refdata
# # End section
 
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation
 
configuration = build_default_calc('jwst', 'miri', 'imaging')
scene = {}
scene['position'] = {'x_offset': 0., 'y_offset': 0., 'orientation': 0., 'position_parameters': ['x_offset', 'y_offset', 'orientation']}
scene['shape'] = {'geometry': 'point'}
scene['spectrum'] = {'name': 'Phoenix Spectrum', 'spectrum_parameters': ['sed', 'normalization']}
scene['spectrum']['sed'] = {'sed_type': 'phoenix', 'key': 'g2v'}
scene['spectrum']['normalization'] = {'type': 'jwst', 'bandpass': 'miri,imaging,f560w', 'norm_flux': 2., 'norm_fluxunit': 'mjy'}
scene['spectrum']['lines'] = []
scene['spectrum']['extinction'] = {'bandpass': 'j', 'law': 'mw_rv_31', 'unit': 'mag', 'value': 0}
configuration['scene'][0] = scene
configuration['configuration']['instrument']['filter'] = 'f1000w'
 
report = perform_calculation(configuration)
print(report['scalar']['sn'])