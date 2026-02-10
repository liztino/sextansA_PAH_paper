#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Find the best match Phoenix Model to the NIRCAM SW filters: F115W, F150W, and F200W

@author: etarantino
"""

import numpy as np
import emcee
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import corner
import os
from astropy.io import ascii

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

save_spec = True
calc_through = True

# set up wavelengths
wave = [1.15409, 1.5009, 1.9875]

# load filter transmission curves
# use the SW NIRCam filters
filt_list = ['F115W', 'F150W', 'F200W']
filts_dir = '/Users/etarantino/Documents/JWST/filts/nircam'


# future Liz-- look at equation 5 in Gordon+ 2022 to do the bandpass math instead of using synphot
def model(theta):
    temp, logg, Z, const  = theta
    ph_model = {}
    # phoenix models are in angstroms for wavelenth
    # and erg cm^-2 s^-1 Angstrom^-1 for flux 
    ph = stsyn.grid_to_spec('phoenix', temp, Z, logg)  
    ph_wave = ph.waveset.value
    ph_flux = ph(ph.waveset).value/(4*np.pi)

    ph_flux = ph_flux * const
    
    ph_model['microns'] = ph_wave/1e4
    ph_model['angstroms'] = ph_wave
    ph_model['Flambda'] = ph_flux
    
    # create a spectrum so we can evaluate
    sp = synphot.SourceSpectrum(Empirical1D, points = ph_model['angstroms'], lookup_table = ph_model['Flambda'] * units.FLAM, keep_neg = True)
    
    # use the SW NIRCam filters
    filt_list = ['F115W', 'F150W', 'F200W']
    
    filt_through = np.zeros(len(filt_list))
    for i, filt_name in enumerate(filt_list):
        
        filts_dir = '/Users/etarantino/Documents/JWST/filts/nircam'
        filt_file = f'{filts_dir}/{filt_name}_mean_system_throughput.txt'
        filt = ascii.read(filt_file)
    
        # create filter object in synphot
        filt_feature = synphot.SpectralElement(Empirical1D, points = filt['Microns']*1e4, lookup_table = filt['Throughput'], keep_neg = True)
    
        # create an observation and calculate the effective stimulation
        obs = synphot.Observation(sp, filt_feature)
        filt_through[i] = obs.effstim(flux_unit='flam').value
    
    return filt_through

def lnlike(theta, x, y, yerr):
    LnLike = -0.5 * np.nansum(((y - model(theta))/yerr)**2)
    return LnLike

def lnprior(theta):
    temp, logg, Z, const = theta
    
    if (temp < 12000) or (temp > 2300):
        if (logg > 0) or (logg < 6):
            if (Z > -0.2) or (Z < 1.2):
                return 0
            else:
                return -np.inf
        else:
            return -np.inf
    else:
        return -np.inf
    
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    print(lp)
    print(theta)
    
    if not np.isfinite(lp):
        return -np.inf
    else:
        return lp + lnlike(theta, x, y, yerr)
        

# load clump file that has the fluxes for each clump
clump_path = '/Users/etarantino/Documents/JWST_DATA_PAHS/analysis/con_sub/clumps/'
clump_file = 'clump_props_pah_filt_flux.txt'
clumps = ascii.read(clump_path + clump_file)

clump_num = clumps['clump_num']

# use the first clump for an example
filts = ['F115W', 'F150W', 'F200W']
flux = np.zeros(len(filts))
flux_err = np.zeros(len(filts))

for i, filt in enumerate(filts):
    flux[i] = clumps[i][filt + '_flux']
    flux_err[i] = clumps[i][filt + '_err']
    
    
initial = np.array([5000, 1, -0.1, 1.5263283323384082e-35])
nwalkers = 10
niter = 5
ndim = 4

data = (wave, flux, flux_err)

p0 = [np.array(initial) + 1e-7 * np.random.randn(ndim) for i in range(nwalkers)]

def main(p0,nwalkers,niter,ndim,lnprob,data):
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data)

    print("Running burn-in...")
    p0, _, _ = sampler.run_mcmc(p0, niter)
    sampler.reset()

    print("Running production...")
    pos, prob, state = sampler.run_mcmc(p0, niter)

    return sampler, pos, prob, state

sampler, pos, prob, state = main(p0,nwalkers,niter,ndim,lnprob,data)

# Analyze Results
samples = sampler.get_chain(discard=0,  flat=True)
# samples[:, 0] = np.round(samples[:, 0] / 100) * 100  # Snap temperature back to discrete values

# Plot results
fig = corner.corner(samples)
plt.show()

