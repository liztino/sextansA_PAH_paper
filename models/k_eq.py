#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Lists the equations for the k continuum subtractiom method so they can be used consistently in other scripts


@author: etarantino
"""

import numpy as np
# from get_pivot_wave import get_pivot_wave


def get_pah_low(f1, f2, f3, lam1, lam2, lam3, k):
    # equation calculated through Liz's derivation
    # for when the PAH contaminating filter is the lower wavelength filter
    # for the F560W/F770W/F1000W combo
    # f1 = fc1 + fp1
    # f2 = fc3 + fp2
    # f3 = fc3
    # fp2 = k*fp1
    
    fp1 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(1 - k - ((lam2 - lam1)/(lam3 - lam1)))
    fp2 = fp1*k
    fc1 = f1 - fp1
    fc2 = (((f3 - fc1)/(lam3 - lam1)) * (lam2 - lam1)) + fc1
    slope = (((f3 - fc1)/(lam3 - lam1)))
    
    consub = {}
    consub = {'pah': fp2, 'con': fc2, 'slope': slope}
    
    return consub

def get_pah_up(f1, f2, f3, lam1, lam2, lam3, k):
    # equation calculated through Liz's derivation
    # for when the PAH contaminating filter is the higher wavelength filter
    # for the F300M/F335M/F360M and F1000W/F1130W/F1500W combo
    # f1 = fc1 
    # f2 = fc3 + fp2
    # f3 = fc3 + fp3
    # fp2 = k*fp3
    
    fp3 = (f1 * (1 - ((lam2 - lam1)/(lam3 - lam1))) + f3 * ((lam2 - lam1)/(lam3 - lam1)) - f2)/(((lam2 - lam1)/(lam3 - lam1)) - k)
    fp2 = fp3*k
    fc3 = f3 - fp3
    fc2 = (((fc3 - f1)/(lam3 - lam1)) * (lam2 - lam1)) + f1
    slope = (((fc3 - f1)/(lam3 - lam1)))
    
    consub = {}
    consub = {'pah': fp2, 'con': fc2, 'slope': slope}
    
    return consub


def get_pah_julia(f1, f2, f3, lam1, lam2, lam3, k):
    # equation calculated through Julia's math
    # Liz tested it and it is equivalent to the pah_low function, but is just a slightly different experession
    # here to compare if needed 
    
    fp1 = (f3 * ((lam2 - lam1)/(lam3-lam1)) + f1 * ((lam3 - lam2)/(lam3 - lam1)) - f2)/(((lam3 - lam2)/(lam3 - lam1)) - k)
    fp2 = fp1*k
    
    print('fp1', fp1)
    print('fp2', fp2)
    
    return fp2



