#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 15:40:09 2024

@author: etarantino
"""
from astropy.io import fits, ascii
import numpy as np 

def get_pivot_wave(filt):
    # load the pivot wavelengths
    pivot_table_file = '../jwst_pivot_wave.txt'
    pivot_table = ascii.read(pivot_table_file)
    filt_ind = np.where(pivot_table['filter'] == filt)[0][0]
    return pivot_table['pivot_wave'][filt_ind]