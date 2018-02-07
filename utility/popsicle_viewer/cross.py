#!/usr/bin/env python
# -*- coding: utf-8 -*-

###########################################
##           cross.py
###########################################
## Routines for calculating cross sections
###########################################

import numpy as np
import numba as nb
from . import constants as const

##############################################################
##############################################################

def get_diff_cross(w,wki,sphmatrix,k_ax,duration,A0):

    k_in = np.sqrt(2.0 * const.muelec * w)
    indx = np.argmin(np.abs(k_ax-k_in))
    
    # Get Fourier Transform of the vector
    # potential
    FRWsq = np.zeros((probk.shape[0]))
    get_fourier_afield(w,wki,duration,FRWsq)

    # The diff cross section
    diff_cross = 4.0 * np.pi**2 * k_ax[indx] * \
                 const.speed_light_au * sphmatrix[:,indx] / \
                 wki / A0 / A0 / FRWsq

    return diff_cross


##############################################################

def get_total_cross(probk,k_ax,wki,wl,duration,bandwidth,A0):
    
    # Get Fourier Transform of the vector
    # potential
    FRWsq = np.zeros((probk.shape[0]))
    get_fourier_afield(wl,wki,duration,FRWsq)
    
    # The total cross section
    total_cross = 4.0 * np.pi**2 * const.speed_light_au / \
                  k_ax * probk / wki / (A0**2) / FRWsq
    
    # Set to zeros everything outside the bandwidth.
    indx = np.where(abs(wl-wki)>bandwidth/1.6 )
    total_cross[indx] = 0.0 #np.nan
    
    return total_cross

##############################################################

@nb.guvectorize([(nb.float64[:],nb.float64[:],nb.float64[:],nb.float64[:])],'(),(n),()->(n)')
def get_fourier_afield(wl,wki,duration,FW):
    
    for i in range(wki.shape[0]):
        if(wki[i]==wl[0]):
            FW[i] = (duration[0] / 4.0)**2
        else:
            FW[i] = 4.0*np.pi**4*np.sin((wl[0]-wki[i])*duration[0]/2.0 )**2 / \
                    (duration[0]**2*(wl[0]-wki[i])**2 - 4.0*np.pi**2)**2 / \
                    (wl[0]-wki[i])**2
            
#################################################################################################
#################################################################################################
