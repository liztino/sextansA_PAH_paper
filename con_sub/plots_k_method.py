#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Quickly plots the relationship between k and flux

@author: etarantino
"""

from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np

# custom functions 
from get_pivot_wave import get_pivot_wave
import k_eq
import k_eq_with_err

plt.ion()

filt_mid = 'F335M'
filt_low = 'F300M'
filt_high = 'F360M'

lam1 = get_pivot_wave(filt_low)
lam2 = get_pivot_wave(filt_mid)
lam3 = get_pivot_wave(filt_high)

f1 = 3
f2 = 10
f3 = 3.5

k = np.linspace(1,10, 1000)

res = k_eq.get_pah_up(f1, f2, f3, lam1, lam2, lam3, k)

plt.figure(1)
plt.clf()
plt.plot(k, res['pah'])
plt.xlabel('k')
plt.ylabel('PAH')

# testing to make the math switch works
# low is the filter that isn't contaminated
filt_pah = 'F770W'
filt_con = 'F1000W'
filt_contam = 'F560W'

lam1 = get_pivot_wave(filt_con)
lam2 = get_pivot_wave(filt_pah)
lam3 = get_pivot_wave(filt_contam)

f1 = 6
f2 = 5
f3 = 3

k = np.linspace(1,10, 1000)

res1 = k_eq.get_pah_up(f1, f2, f3, lam1, lam2, lam3, k)


plt.figure(2)
plt.clf()
plt.plot(k, res1['pah'], label ='Universal', alpha = 0.5)

filt_mid = 'F770W'
filt_low = 'F560W'
filt_high = 'F1000W'

lam1 = get_pivot_wave(filt_low)
lam2 = get_pivot_wave(filt_mid)
lam3 = get_pivot_wave(filt_high)

f1 = 3
f2 = 5
f3 = 6

k = np.linspace(1,10, 1000)

res2 = k_eq.get_pah_low(f1, f2, f3, lam1, lam2, lam3, k)

# now use the real one
plt.plot(k, res2['pah'], label ='Old', alpha = 0.5, ls = '--')

plt.xlabel('k')
plt.ylabel('PAH')

