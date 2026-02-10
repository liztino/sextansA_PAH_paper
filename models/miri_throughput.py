#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates the transmission functions for the MIRI filters through the pandeia engine
Need to do this because ST makes everything complicated...

@author: etarantino
"""
import numpy as np 
import matplotlib.pyplot as plt
from pandeia.engine.instrument_factory import InstrumentFactory


# tested to make sure that the detector parameters do not matter for the F770W filter 
# something is super broken with the F2550W filter, leaving it out for now...

# filters = ['F560W', 'F770W', 'F1000W', 'F1130W', 'F1280W', 'F1500W']
filters = ['F1280W']
lower = np.array([3, 5, 8, 9, 10, 12])
upper = np.array([8, 10, 12, 13, 15, 18])

lower[0] = 8
upper[0] = 18

plt.figure(1)
plt.clf()

for i in range(len(filters)):
    filt = filters[i]
    filt_low = filt.lower()
    print('Working on ', filt)

    # MIRI Imaging
    conf= {
        "detector": {
            "nexp": 1,
            "ngroup": 25,
            "nint": 15,
            "readout_pattern": "fastr1",
            "subarray": "full"
        },
        "instrument": {
            "aperture": "imager",
            "filter": f"{filt_low}",
            "instrument": "miri",
            "mode": "imaging"
        },
    }
     
    # set up your wavelengths
    wave = np.arange(lower[i], upper[i], 0.001) 
     
    # create a configured instrument
    instrument_factory = InstrumentFactory(config=conf)
     
    # where conf is a configuration dictionary for a calculation:
    #     see below for examples
     
    # get the throughput of the instrument over the desired wavelength range
    eff = instrument_factory.get_total_eff(wave)
    print(eff)

    # plt.plot(wave, eff, label = f'{filt}', ls = '--')

    ind = np.where(eff != 0)[0]
    print(ind)
    delta = 10
    i1 = ind[0] - delta
    i2 = ind[-1] + delta

    final_eff = eff[i1:i2]
    final_wave = wave[i1:i2] 

    plt.plot(final_wave, final_eff, label = f'{filt}', ls = '--')

    savedir = '/Users/etarantino/Documents/JWST/filts/miri/'
    savename = f'{filt}_mean_system_throughput.txt'

    data = np.array([final_wave, final_eff]).T

    np.savetxt(savedir + savename, data, header = 'Microns  Throughput')
    

plt.legend(loc = 'best')
plt.show()

