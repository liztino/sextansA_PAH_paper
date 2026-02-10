#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Takes the D21 and BH22 astrodust + PAH models for the dust and PAH signal
Create grid of Phoenix models for the stellar continuum 
Uses the JWST throughput curves to find prediction for JWST data 

@author: etarantino
"""
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import scipy.interpolate as inter

os.environ['PYSYN_CDBS']
import pysynphot as S
import stsynphot as stsyn
from synphot import SourceSpectrum, SpectralElement, units
from synphot.models import Empirical1D
import synphot

plt.ion()

# constants
c = 3e14        # in microns/sec
c_A = 3e18      # in Angstrom/sec
h = 6.62e-27    # in erg s 

# column density in Sextans A to normalize the Draine models
# Max integrated HI flux near reg in Jy/beam * m/s 
HI_flux = 250

# beam in arcseconds
bmaj = 7.6
bmin = 6.5

# calculate the column density 
N_HI = (1.105e21 * HI_flux)/(bmaj * bmin)

# returns the corresponding D21 model in microns and MJy/sr
def get_D21_model(qpah, U, ion, size, N_HI):
    D21_model = {}
    
    # load Draine model
    drainepath = '/Users/etarantino/Documents/PAHs/Draine2021_models/'
    model = 'BC03_Z0.0004_10Myr'
    # drainename = 'pahspec.out_bc03_z0.0004_1e7_0.50_st_sma'
    drainename = 'pahspec.out_bc03_z0.0004_1e7_{:01.2f}_{:s}_{:s}'.format(U, ion, size)
    header = ['wave', 'total',   'Astrodust',   'PAH^+',    'PAH^0']
    
    # D21 models are given in 4pi nu jnu -> erg s^-1 H^-1
    draine_data = ascii.read(drainepath + model + '/' + drainename, names = header, data_start = 7)
    
    orig_qpah = 3.51
    factor = qpah/orig_qpah
    
    # convert from nu * P_nu to P_nu
    Hz = c/(draine_data['wave'])
    
    # extract and convert from erg/s/H to erg/s/cm^2/Hz/sr
    astrodust = (draine_data['Astrodust'].value/(Hz * 4 * np.pi)) * N_HI
    pah_pos = (draine_data['PAH^+'].value/(Hz * 4 * np.pi)) * N_HI
    pah_neut = (draine_data['PAH^0'].value/(Hz * 4 * np.pi)) * N_HI
    total = (draine_data['total'].value/(Hz * 4 * np.pi)) * N_HI
    
    # convert to the given qpah
    new_pah = factor * pah_pos + factor * pah_neut
    new_tot = new_pah + astrodust
    
    draine_wave = draine_data['wave']
    
    # convert to MJy
    erg_MJy = 1e-17
    draine_flux = new_tot/erg_MJy
    
    # plt.plot(draine_wave, draine_flux, label = 'U={:01.2f}'.format(U))
    # plt.xlim(0, 30)
    # plt.ylim(0,10)
    # plt.title('D21 Model MJy')
    
    cutoff = np.where(draine_wave < 30)[0]
    
    D21_model['microns'] = draine_wave[cutoff]
    D21_model['MJy_sr'] = draine_flux[cutoff]
    
    # convert to Flambda 
    angstrom = draine_data['wave'].value*1e4
    Flambda = new_tot * c_A/(angstrom)**2
    
    D21_model['Flambda'] = Flambda[cutoff]
    D21_model['angstroms'] = angstrom[cutoff]
    
    return D21_model

# normalizing values so the models are in the right ballpark
# in MJy/sr
D21_norm = 0.2
phoenix_norm = 0.05

# convert to Flam units to stay consistent
D21_norm_Flam = (D21_norm * 1e6 * u.Jy).to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(10 * u.micron))

U_list = np.arange(0, 7.5, 0.5)
qpah_list = np.arange(0, 2.5, 0.5)

# for i in range(len(U_list)):

#     plt.figure()
#     plt.clf()
#     plt.title('U={:01.2f}'.format(U_list[i]))
    
#     for j in range(len(qpah_list)):
#         D21_model = get_D21_model(qpah_list[j], U_list[i], 'hi', 'lrg', N_HI)
        
#         plt.plot(D21_model['microns'], D21_model['Flambda'])
        
#     plt.xlim(0, 20)
    

D21_model = get_D21_model(0.5, 4, 'hi', 'lrg', 1e19) 

temp = 5500
logg = 5.0
Z = -0.2
# phoenix models are in angstroms for wavelenth
# and erg cm^-2 s^-1 Angstrom^-1 for flux 
ph = stsyn.grid_to_spec('phoenix', temp, Z, logg)  
ph_wave = ph.waveset.value
ph_flux = ph(ph.waveset).value/(4*np.pi)

mod_fac = 6.233570019315545e-25
ph_flux = ph_flux * mod_fac

plt.figure(2)
plt.clf()
plt.plot(ph_wave, ph_flux)

plt.figure(3)
plt.clf()
plt.plot(D21_model['angstroms'], D21_model['Flambda'])


# interpolate from D21 data
ph_interp = inter.CubicSpline(ph_wave, ph_flux)

# stop interpolation at 10 micron
stop = 10
ind = np.where(D21_model['angstroms'] > stop * 1e4)[0]
vals = D21_model['angstroms'] [ind]
ph_vals = np.zeros(np.shape(D21_model['angstroms']))
ph_vals[ind] = ph_interp(vals)

tot_mod = ph_interp(D21_model['angstroms']) +  D21_model['Flambda']


plt.figure(1)
plt.clf()
plt.plot(D21_model['angstroms'], D21_model['Flambda'], color = 'b', alpha = 0.7, label = 'D21 Model')
plt.plot(ph_wave, ph_flux, color = 'r', alpha = 0.7, label = 'Phoenix Model')
plt.plot(D21_model['angstroms'], ph_interp(D21_model['angstroms']), color = 'purple', alpha = 0.7, label = 'Phoenix Interpolated')
plt.plot(D21_model['angstroms'], tot_mod, c = 'k', alpha = 0.7, label = 'Total')
plt.legend(loc = 'best')
plt.xlim(0, 200000)
# plt.ylim(0,2e-5)




savedir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/grids/D21_new_phoenix/'
modeldir = savedir + 'spec/'
throughdir = savedir + 'throughput/'

if not os.path.exists(throughdir):
    os.makedirs(throughdir)    

if not os.path.exists(modeldir):
    os.makedirs(modeldir)    

save_spec = False
calc_through = False

regname = 'SexA_pah_small_box'
regdir = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/models/obs/'
regfile = regdir + f'sexa_SED_region_{regname}.txt'
jwst_data = ascii.read(regfile)

miri_filts = {'F1500W', 'F1000W', 'F560W', 'F770W', 'F1130W'}
nircam_filts = {'F115W', 'F150W', 'F200W', 'F300M', 'F335M', 'F360M'}

filt_list = ['F300M', 'F335M', 'F360M', 'F560W', 'F770W', 'F1000W', 'F1130W', 'F1500W']
MJy_flux = np.zeros(len(filt_list))
Flam_flux = np.zeros(len(filt_list))
micron_mid = np.zeros(len(filt_list))
filt_dict = {}
jwst_dict = {}
jwst_dict_Jy = {}
filt_wave = np.zeros(len(filt_list))

for i, filt_name in enumerate(filt_list):
    # load filter info
    if filt_name in nircam_filts:
        filts_dir = '/Users/etarantino/Documents/JWST/filts/nircam/'
    else:
        filts_dir = '/Users/etarantino/Documents/JWST/filts/miri/'
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
    jwst_dict_Jy[filt_name] = MJy_flux
    
    filt_wave[i] = filt_feature.pivot().value
    
# normalizing values so the models are in the right ballpark
# in MJy/sr
D21_norm = 0.2

# convert to Flam units to stay consistent
D21_norm_Flam = (D21_norm * 1e6 * u.Jy).to(u.erg / u.s / u.AA / u.cm**2, equivalencies=u.spectral_density(10 * u.micron))

# setup phoenix model so it can be added to the future D21 models
# best fit phoenix model was found through get_phoenix_models.py script
temp = 5500
logg = 5.0
Z = -0.2
# phoenix models are in angstroms for wavelenth
# and erg cm^-2 s^-1 Angstrom^-1 for flux 
ph = stsyn.grid_to_spec('phoenix', temp, Z, logg)  
ph_wave = ph.waveset.value
ph_flux = ph(ph.waveset).value/(4*np.pi)

mod_fac = 6.233570019315545e-25
ph_flux = ph_flux * mod_fac

# setup interpolation
ph_interp = inter.CubicSpline(ph_wave, ph_flux)

# grab a D21 model so we can use the x-axis
D21_model = get_D21_model(0.5, 4, 'hi', 'lrg', 1e19) 
xx = D21_model['angstroms']

# interpolate over the D21 model x-axis
ph_xx = ph_interp(xx)

# U_list = np.arange(0, 7.5, 0.5)
# qpah_list = np.arange(0, 2.5, 0.5)
U_list = np.arange(0, 5, 0.5)
qpah_list = np.arange(0, 2, 0.5)
ion_list = ['lo', 'st', 'hi']
size_list = ['sma', 'std', 'lrg']

chisq_list = []
chisq_Jy = []
U_master = []
qpah_master = []
const_master = []
ion_master = []
size_master = []
for i in range(len(U_list)):
    # get initial model to find good normalization list 
    D21_model = get_D21_model(qpah_list[0], U_list[i], ion_list[0], size_list[0], 1)
    max_val = np.nanmax(D21_model['Flambda'])
    norm = (D21_norm_Flam/max_val).value
    
    # generate range of normalization values 
    # const_list = np.logspace(np.log10(0.1 * norm), np.log10(10*norm), 10)
    const_list = np.logspace(19, 23, 40)
    
    # plt.figure()
    # plt.clf()
    # plt.title('T={:4.0f}'.format(temp_list[i]))
    
    for j in range(len(qpah_list)):
        for k in range(len(const_list)):
            for l in range(len(ion_list)):
                for m in range(len(size_list)):

                    # create model
                    D21_model = get_D21_model(qpah_list[j], U_list[i], ion_list[l], size_list[m], const_list[k]) 
                    full_model = {}
                    full_model['microns'] = D21_model['microns']
                    full_model['angstroms'] = D21_model['angstroms']
                    full_model['D21_MJy'] = D21_model['MJy_sr']
                    full_model['D21_Flam'] = D21_model['Flambda']
                    full_model['phoenix'] = ph_xx
                    full_model['total'] = ph_xx + D21_model['Flambda']
                    
                    if save_spec:
                        # save model
                        full_name = 'D21_combo_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.txt'.format(U_list[i], qpah_list[j], np.log10(const_list[k]), ion_list[l], size_list[m])
                        ascii.write(full_model, modeldir + full_name, overwrite = True )
                        
                    if calc_through:
                        # turn the total spectrum into a source spectrum object in synphot 
                        sp = synphot.SourceSpectrum(Empirical1D, points = full_model['angstroms'], lookup_table = full_model['total'] * units.FLAM, keep_neg = True)
                        
                        # loop through all filters
                        filt_through = np.zeros(len(filt_list))
                        filt_through_Jy = np.zeros(len(filt_list))
                        jwst_vals = np.zeros(len(filt_list))
                        jwst_vals_Jy = np.zeros(len(filt_list))
                        for f, filt_name in enumerate(filt_list):
                            # create an observation and calculate the effective stimulation
                            obs = synphot.Observation(sp, filt_dict[filt_name])
                            filt_through[f] = obs.effstim(flux_unit='flam').value
                            filt_through_Jy[f] = (filt_through[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
                            filt_through_Jy[f] = filt_through_Jy[f]/1e6
                            
                            # get jwst value 
                            jwst_vals[f] = jwst_dict[filt_name].value
                            jwst_vals_Jy[f] = jwst_dict_Jy[filt_name].value
                            
                        # save the filter information
                        through_array = np.array([filt_list, filt_through]).T
                        through_name = 'D21_combo_throughput_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.txt'.format(U_list[i], qpah_list[j], np.log10(const_list[k]), ion_list[l], size_list[m])
                        ascii.write(through_array, throughdir + through_name, names = ['filter', 'throughput'], overwrite = True)
                        
                        # calculate the chisquared 
                        chisq_list.append(np.nansum((jwst_vals - filt_through)**2/(filt_through)))
                        chisq_Jy.append(np.nansum((jwst_vals_Jy - filt_through_Jy)**2/(filt_through_Jy)))

                        # print(jwst_vals)
                        # print(filt_through)
                        # print(np.nansum((jwst_vals - filt_through)**2/(filt_through)))
                        
                        
                        U_master.append(U_list[i])
                        qpah_master.append(qpah_list[j])
                        const_master.append(const_list[k])
                        ion_master.append(ion_list[l])
                        size_master.append(size_list[m])
                        
                    else:
                        through_name =  'D21_combo_throughput_{:01.2f}_{:01.2f}_{:3.1f}_{:s}_{:s}.txt'.format(U_list[i], qpah_list[j], np.log10(const_list[k]), ion_list[l], size_list[m])
                        through_data = ascii.read(throughdir + through_name)
                        
                        # loop through all filters
                        filt_through = np.zeros(len(filt_list))
                        filt_through_Jy = np.zeros(len(filt_list))
                        jwst_vals = np.zeros(len(filt_list))
                        jwst_vals_Jy = np.zeros(len(filt_list))
                        for f, filt_name in enumerate(filt_list):
                            ind = np.where(through_data['filter'] == filt_name)[0]
                            filt_through[f] = through_data['throughput'][ind]
                            
                            filt_through_Jy[f] = (filt_through[f] * u.erg / u.s / u.AA / u.cm**2).to(u.Jy, equivalencies=u.spectral_density(filt_wave[f] * u.AA)).value
                            filt_through_Jy[f] = filt_through_Jy[f]/1e6
                            
                            # get jwst value 
                            jwst_vals[f] = jwst_dict[filt_name].value
                            jwst_vals_Jy[f] = jwst_dict_Jy[filt_name].value
                        
                        # calculate the chisquared 
                        # chisq_list.append(np.nansum((jwst_vals - filt_through)**2/(filt_through)))
                        # chisq_Jy.append(np.nansum((jwst_vals_Jy - filt_through_Jy)**2/(filt_through_Jy)))
                        
                        # calc reduced chisquared
                        param = 5
                        chisq_list.append(np.nansum((jwst_vals - filt_through)**2)/(len(filt_through) - param))
                        chisq_Jy.append(np.nansum((jwst_vals_Jy - filt_through_Jy)**2)/(len(filt_through_Jy) - param))
    

                        # print(jwst_vals)
                        # print(filt_through)
                        # print(np.nansum((jwst_vals - filt_through)**2/(filt_through)))
                        
                        
                        U_master.append(U_list[i])
                        qpah_master.append(qpah_list[j])
                        const_master.append(const_list[k])
                        ion_master.append(ion_list[l])
                        size_master.append(size_list[m])
                

final_table = np.array([U_master, qpah_master, const_master, ion_master, size_master, chisq_list, chisq_Jy]).T
final_table_name = f'D21_combo_{regname}_models_master_table_v2.txt'
ascii.write(final_table, savedir + final_table_name, names = ['U', 'qpah', 'const', 'ion', 'size', 'chisq', 'chisq_Jy'], overwrite = True)            
                

