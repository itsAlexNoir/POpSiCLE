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
    FRWsq = np.zeros((probk.shape[0]))
    get_fourier_afield(w,wki,duration,FRWsq)

    # The diff cross section
    diff_cross = 4.0 * np.pi**2 * k_ax[indx] * \
                 const.speed_light_au * sphmatrix[:,indx] / \
                 wki / A0 / A0 / FRWsq

    return diff_cross


##############################################################

def get_total_cross(probk,k_ax,w,wki,duration,bandwidth,A0):

    # Get Fourier Transform of the vector
    # potential
    FRWsq = np.zeros((probk.shape[0]))
    get_fourier_afield(w,wki,duration,FRWsq)

    # The total cross section
    total_cross = 4.0 * np.pi**2 * k_ax * \
                 const.speed_light_au * probk / \
                 wki / A0 / A0 / FRWsq

    # Set to zeros everything outside the bandwidth.
    indx = np.where(abs(w-wki)>bandwidth/3.0 )
    total_cross[indx] = 0.0 #np.nan
    
    return total_cross

##############################################################

@nb.guvectorize([(nb.float64[:],nb.float64[:],nb.float64[:],nb.float64[:])],'(n),(),()->(n)')
def get_fourier_afield(w,wki,duration,FW):

    for i in range(w.shape[0]):
        if(w[i]==wki[0]):
            FW[i] = (duration[0] / 4.0)**2
        else:
            FW[i] = 4.0*np.pi**4*np.sin((w[i]-wki[0])*duration[0]/2.0 )**2 / \
                 (duration[0]**2*(w[i]-wki[0])**2 - 4.0*np.pi**2)**2 / \
                 (w[i]-wki[0])**2

#################################################################################################
#################################################################################################
