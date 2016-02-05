#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################
##           cross.py
###########################################
## Routines for calculating cross sections
###########################################

import numpy as np
import constants as const

##############################################################
##############################################################

def get_diff_cross(w,wki,sphmatrix,k_ax,duration,A0):

    k_in = np.sqrt(2.0 * const.muelec * w)
    indx = np.argmin(np.abs(k_ax-k_in))
    
    # Get Fourier Transform of the vector
    # potential
    FRWsq = get_fourier_afield(w,wki,duration)

    # The diff cross section
    diff_cross = 4.0 * np.pi**2 * k_ax[indx] * \
                 const.speed_light_au * sphmatrix[:,indx] / \
                 wki / A0 / A0 / FRWsq

    return diff_cross


##############################################################

def get_fourier_afield(w,wki,duration):

    if(w==wki):
        FW = (duration / 4.0)**2
    else:
        FW = 4.0*np.pi**4*np.sin((w-wki)*duration/2.0 )**2 / \
             (duration**2*(w-wki)**2 - 4.0*np.pi**2)**2 / \
             (w-wki)**2

    return FW
