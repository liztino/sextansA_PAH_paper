#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculates the TIR from Shi+ 2014 paper and compares to the sum of PAH features

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import astropy.units as u
import scipy.interpolate as inter

# distance to Sextans A
# in Mpc, from McQuinn+ 2017 using TRGB
d = 1.46

# 1 pc to cm 
pc = 3.086e18

d_cm = d * 1e6 * pc

# speed of light in microns
c = 3e14

# luminosity of the sun in erg/s
L_sol = 3.846e33

# from the Shi+ 2014 paper, using "star-forming region 3"
# all in mJy units

# Spitzer table
sf_8 =  0.6
sf_8_unc = 0.05

sf_24 = 6.36
sf_24_unc = 0.65

# herschel table
sf_70 =  265
sf_70_unc = 4

sf_160 = 296
sf_160_unc = 24

# function to convert mjy to a luminosity 
# mjy->erg/s/Hz
def mJy_to_L(mJy):
    return 1e-26 * 4 * np.pi * d_cm**2 * mJy

# function to convert wavelength in microns to frequency in Hz
def um_to_hz(um):
    return c/um

# using wavelength weights from D+L2007
hz_8 = um_to_hz(7.9)
hz_24 = um_to_hz(24)
hz_70 = um_to_hz(71)
hz_160 = um_to_hz(160)

# function to calculate the TIR from Draine + Li 2007 equation 22
# also what Hunt+ 2010 uses to calculate their TIR in BCDs
# unclear what units are?

def L_TIR(flux_8, flux_24, flux_70, flux_160):
    return (0.95 * hz_8 * mJy_to_L(flux_8)) + 1.15*(hz_24 * mJy_to_L(flux_24)) + (hz_70 * mJy_to_L(flux_70)) + (hz_160 * mJy_to_L(flux_160))

print(L_TIR(sf_8, sf_24, sf_70, sf_160)/L_sol)


