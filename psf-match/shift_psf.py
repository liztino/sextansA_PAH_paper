#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Single row and column shift for a PSF

@author: Eliz
"""
import numpy as np

def shift_psf(psf, n = 1):
    pad = len(psf) + n
    psf_pad = np.zeros((pad,pad))
    psf_pad[n:,n:] = psf

    return psf_pad