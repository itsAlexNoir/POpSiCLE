#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################
##           cross.py
###########################################
## Routines for calculating cross sections
###########################################

import numpy as np
import constants as const
import numba as nb

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

def get_total_cross(probk,k_ax,w,wki,duration,bandwidth,A0):

    # Get Fourier Transform of the vector
    # potential
    FRWsq = get_fourier_afield(w,wki,duration)
    
    # The total cross section
    total_cross = 4.0 * np.pi**2 * k_ax * \
                 const.speed_light_au * probk / \
                 wki / A0 / A0 / FRWsq

    # Set to zeros everything outside the bandwidth.
    indx = np.where(abs(w-wki)>bandwidth/3.0 )
    total_cross[indx] = np.nan
    
    return total_cross

##############################################################

def get_fourier_afield(w,wki,duration):

    if(w==wki):
        FW = (duration / 4.0)**2
    else:
        FW = 4.0*np.pi**4*np.sin((w-wki)*duration/2.0 )**2 / \
            (duration**2*(w-wki)**2 - 4.0*np.pi**2)**2 / \
            (w-wki)**2
        
    return FW

#################################

def get_fourier_afield(w,wki,duration):

    FW = 4.0*np.pi**4*np.sin((w-wki)*duration/2.0 )**2 / \
         (duration**2*(w-wki)**2 - 4.0*np.pi**2)**2 / \
         (w-wki)**2

    indx = np.where(w==wki)
    if(np.shape(indx)[1]!=0):
        FW = (duration / 4.0)**2

    return FW
