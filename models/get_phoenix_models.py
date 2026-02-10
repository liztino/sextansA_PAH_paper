#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Find the best match Phoenix Model to the NIRCAM SW filters: F115W, F150W, and F200W

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u

os.environ['PYSYN_CDBS']
import stsynphot as stsyn
from synphot.models import Empirical1D
from synphot import units
import synphot

plt.ion()

# constants
c = 3e14        # in microns/sec
c_A = 3e18      # in Angstrom/sec
h = 6.62e-27    # in erg s 

# parameters for phoenix models
temp = 3000
logg = 4.0 
Z = -0.2

save_spec = False
calc_through = False


regname = 'SexA_pah_small_box'
regdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/obs/'
regfile = regdir + f'sexa_SED_region_{regname}.txt'
jwst_data = ascii.read(regfile)


savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/phoenix_more_filts/'
modeldir = savedir + 'spec/'
throughdir = savedir + 'throughput/'

if not os.path.exists(throughdir):
    os.makedirs(throughdir)    

if not os.path.exists(modeldir):
    os.makedirs(modeldir)    
  

def get_phoenix(temp, logg, Z, const):
    phoenix = {}
    # phoenix models are in angstroms for wavelenth
    # and erg cm^-2 s^-1 Angstrom^-1 for flux 
    ph = stsyn.grid_to_spec('phoenix', temp, Z, logg)  
    ph_wave = ph.waveset.value
    ph_flux = ph(ph.waveset).value/(4*np.pi)

    ph_flux = ph_flux * const
    
    phoenix['microns'] = ph_wave/1e4
    phoenix['angstroms'] = ph_wave
    phoenix['Flambda'] = ph_flux
    
    return phoenix

# setting up the filters
# just using the SW nircam filters
filt_list = ['F115W', 'F150W', 'F200W']
filt_dict = {}
jwst_dict = {}

for filt_name in filt_list:
    # load filter info
    filts_dir = '/Users/etarantino/Documents/JWST/filts/nircam/'
    filt_file = f'{filts_dir}/{filt_name}_mean_system_throughput.txt'
    filt = ascii.read(filt_file)

    # create filter object in synphot
    filt_feature = synphot.SpectralElement(Empirical1D, points = filt['Microns']*1e4, lookup_table = filt['Throughput'], keep_neg = True)

    # save to the filter dictionary
    filt_dict[filt_name] = filt_feature
    
    # find the JWST for the filter
    ind = np.where(jwst_data['filt'] == filt_name)[0]
    MJy_flux = jwst_data['avg'][ind]
    
    # convert the flux to Flam units then save to a dictionary to use later
    Flam_flux = (MJy_flux * 1e6 * u.Jy).to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(filt_feature.pivot()))
    jwst_dict[filt_name] = Flam_flux
    

# normalizing values so the models are in the right ballpark
# in MJy/sr
phoenix_norm = 0.05

# convert to Flam units to stay consistent
phoenix_norm_Flam = (phoenix_norm * 1e6 * u.Jy).to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(1 * u.micron))

temp_list = np.arange(5000, 6000, 500)
logg_list = np.arange(0, 5.5, 0.5)

chisq_list = []
temp_master = []
logg_master = []
const_master = []
for i in range(len(temp_list)):
    # get initial model to find good normalization list 
    ph_model = get_phoenix(temp_list[i], logg_list[0], Z, 1)
    max_val = np.nanmax(ph_model['Flambda'])
    norm = (phoenix_norm_Flam/max_val).value
    
    # generate range of normalization values 
    # const_list = np.logspace(np.log10(0.1 * norm), np.log10(10*norm), 10)
    const_list = np.logspace(-24.5, -23.7, 20)
    # plt.figure()
    # plt.clf()
    # plt.title('T={:4.0f}'.format(temp_list[i]))
    
    for j in range(len(logg_list)):
        for k in range(len(const_list)):
            # create model
            # ph_model = get_phoenix(temp_list[i], logg_list[j], Z, const_list[k])
            
            if save_spec:
                # save model
                ph_name = 'phoenix_{:4.0f}_{:2.1f}_{:3.2f}.txt'.format(temp_list[i], logg_list[j], np.log10(const_list[k]))
                ascii.write(ph_model, modeldir + ph_name, overwrite = True )
                
            if calc_through:
                # turn the total spectrum into a source spectrum object in synphot 
                sp = synphot.SourceSpectrum(Empirical1D, points = ph_model['angstroms'], lookup_table = ph_model['Flambda'] * units.FLAM, keep_neg = True)
                
                # loop through all filters
                filt_through = np.zeros(len(filt_list))
                jwst_vals = np.zeros(len(filt_list))
                for l, filt_name in enumerate(filt_list):
                    # create an observation and calculate the effective stimulation
                    obs = synphot.Observation(sp, filt_dict[filt_name])
                    filt_through[l] = obs.effstim(flux_unit='flam').value
                    
                    # get jwst value 
                    jwst_vals[l] = jwst_dict[filt_name].value
                    
                # save the filter information
                through_array = np.array([filt_list, filt_through]).T
                through_name = 'phoenix_throughput_{:4.0f}_{:2.1f}_{:3.2f}.txt'.format(temp_list[i], logg_list[j], np.log10(const_list[k]))
                ascii.write(through_array, throughdir + through_name, names = ['filter', 'throughput'], overwrite = True)
                
                # calculate the chisquared 
                chisq_list.append(np.nansum((jwst_vals - filt_through)**2/(filt_through)))
                # print(jwst_vals)
                # print(filt_through)
                # print(np.nansum((jwst_vals - filt_through)**2/(filt_through)))
                
                # create array with all the parameters 
                temp_master.append(temp_list[i])
                logg_master.append(logg_list[j])
                const_master.append(const_list[k])
            else:
                through_name = 'phoenix_throughput_{:4.0f}_{:2.1f}_{:3.2f}.txt'.format(temp_list[i], logg_list[j], np.log10(const_list[k]))
                through_data = ascii.read(throughdir + through_name)
                
                filt_through = np.zeros(len(filt_list))
                jwst_vals = np.zeros(len(filt_list))
                for l, filt_name in enumerate(filt_list):
                    ind = np.where(through_data['filter'] == filt_name)[0]
                    filt_through[l] = through_data['throughput'][ind]
                    
                    # get jwst value 
                    jwst_vals[l] = jwst_dict[filt_name].value
                    
                # calculate the chisquared 
                chisq_list.append(np.nansum((jwst_vals - filt_through)**2/(filt_through)))
                
                # print(jwst_vals)
                # print(filt_through)
                # print(np.nansum((jwst_vals - filt_through)**2/(filt_through)))
                
                # create array with all the parameters 
                temp_master.append(temp_list[i])
                logg_master.append(logg_list[j])
                const_master.append(const_list[k])
                

final_table = np.array([temp_master, logg_master, const_master, chisq_list]).T
final_table_name = f'phoenix_models_{regname}_master_table_Jy.txt'
ascii.write(final_table, savedir + final_table_name, names = ['T', 'logg', 'Const', 'Chisq'], overwrite = True)            
                
                

            

