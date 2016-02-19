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

@nb.jit('void(float64[:],float64[:],float64[:],float64[:,:])',target='cpu')
def buildjacobian(a,b,c,jacobian):
    for j in range(b.shape[0]):
        for i in range(a.shape[0]):
            jacobian[i,j] = a[i] * b[j] * c[j]

@nb.jit('void(float64[:],float64[:,:],float64[:,:])',target='cpu')
def integrate1D_over_2D(intfunc,func,jacobian):
    for i in range(func.shape[1]):
        for j in range(func.shape[0]):
		intfunc[i] = intfunc[i] + func[j,i] * jacobian[j,i]

@nb.jit('void(float64[:],float64[:,:,:],float64[:,:])',target='cpu')
def integrate1D_over_3D(intfunc,func,jacobian):
    for i in range(func.shape[2]):
        for j in range(func.shape[1]):
            for k in range(func.shape[0]):
		intfunc[i] = intfunc[i] + func[k,j,i] * jacobian[j,i]

@nb.jit('void(float64[:],float64[:,:,:],float64[:,:])',target='cpu')
def integrate2D_over_3D(intfunc,func,jacobian):
    for i in range(func.shape[2]):
        for j in range(func.shape[1]):
            for k in range(func.shape[0]):
		intfunc[j,i] = intfunc[j,i] + func[k,j,i] * jacobian[j,i]

# @nb.jit('void(float64[:,:],float64[:,:,:])',target='cpu')
# def integrateOverTheta3D(prob,probpad):
#     for i in range(prob.shape[0]):
#         for j in range(probpad.shape[0]):
#         	for k in range(probpad.shape[1]):
# 	            prob[i,k] = prob[i,k] + probpad[i,j,k]

@nb.jit('void(float64[:,:],float64[:],float64[:])',target='cpu')
def create_sphjacobian(sphjacobian,rpts,sintheta):
    for j in range(rpts.shape[0]):
        for i in range(sintheta.shape[0]):
            sphjacobian[i,j] = rpts[j] * rpts[j] * sintheta[i]  

# @nb.jit('void(float64[:,:],float64[:],float64[:])',target='cpu')
# def create_jacobEnEe(jacobEnEe,elec_ax,nuc_ax):
#     for i in range(nuc_ax.shape[0]):
#         for j in range(elec_ax.shape[0]):
#             jacobEnEe[i,j] = jacobEnEe[i,j] / nuc_ax[j] / elec_ax[i]  
            
# @nb.jit('void(float64[:,:],float64[:,:])',target='cpu')
# def add_jacobEnEe(probEnEe,jacobEnEe):
#     for i in range(probEnEe.shape[0]):
#         for j in range(probEnEe.shape[1]):
#             probEnEe[i,j] = probEnEe[i,j] * jacobEnEe[i,j]

