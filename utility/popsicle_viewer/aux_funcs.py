#!/usr/bin/env python
# -*- coding: utf-8 -*-


########################################
##           aux_funcs.py
########################################
## This file contains the definition
## of the class axes.
########################################

import numpy as np
import numba as nb

#################################################################################
#################################################################################

@nb.jit(nopython=True)
def build_jacobian(a,b,c,jacobian):
    for j in range(b.shape[0]):
        for i in range(a.shape[0]):
            jacobian[i,j] = a[i] * b[j] * c[j]

@nb.jit(nopython=True)
def create_sphjacobian(r_ax,sintheta,jacobian):
    for i in range(sintheta.shape[0]):
        for j in range(r_ax.shape[0]):
            jacobian[i,j] = r_ax[j] * r_ax[j] * sintheta[i]  
            
@nb.jit(nopython=True)
def spherical_azimuthal_integral(sphwave,suma,weights):
    for i in range(sphwave.shape[0]):
        for j in range(sphwave.shape[1]):
            for k in range(sphwave.shape[2]):
                suma[j,k] = suma[j,k] + sphwave[i,j,k] * weights[i]

@nb.jit(nopython=True)
def spherical_angular_integral(sphwave,suma,jacobian,weights):
    for i in range(sphwave.shape[0]):
        for j in range(sphwave.shape[1]):
            for k in range(sphwave.shape[2]):
                suma[k] = suma[k] + sphwave[i,j,k] * jacobian[j,k] * \
                          weights[i,j]

@nb.jit(nopython=True)
def spherical_polar_integral(sphwave,suma,jacobian,weights):
    for i in range(sphwave.shape[0]):
        for j in range(sphwave.shape[2]):
            suma[j] = suma[j] + sphwave[i,j] * jacobian[i,j] * weights[i]
