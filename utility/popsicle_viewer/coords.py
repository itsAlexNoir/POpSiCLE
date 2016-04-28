#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.interpolate as interp
import constants as const
import aux_funcs as aux

########################################
##           coords.py
## This module contains functions for
## transforming coordinates.
########################################


###################################
########### cyl2sph ###############
###################################

def cyl2sph2D(cylfunc,rho,z,r,theta):
    print('Transforming from cylindrical to spherical coordinates...')
    RHO,Z = np.meshgrid(rho,z)
    R_NON = np.sqrt(RHO**2 + Z**2)
    THETA_NON = np.arctan2(RHO, Z)
    #RX,THETAX = np.meshgrid(r,theta)
    
    sphfunc = interp.griddata((np.ravel(R_NON),np.ravel(THETA_NON)),
	    np.ravel(cylfunc),(r[None,:],theta[:,None]),
	    method='cubic',fill_value=1e-30)

    return sphfunc

def cart2sph2Dpolar(cartfunc,x,z,r,theta):
    print('Transforming from cartesian to spherical coordinates...')
    X,Z = np.meshgrid(x,z)
    R_NON = np.sqrt(X**2 + Z**2)
    THETA_NON = np.arctan2(X,Z)
    #R,THETA, PHI = np.meshgrid(r,theta,phi)
    
    sphfunc = interp.griddata((np.ravel(R_NON),np.ravel(THETA_NON)),
    	np.ravel(cartfunc),(r[None,:],theta[:,None],),
    	method='cubic',fill_value=1e-30)

    return sphfunc

def cart2sph2Dazimuthal(cartfunc,x,y,r,phi):
    print('Transforming from cartesian to spherical coordinates...')
    X,Y = np.meshgrid(x,y)
    R_NON = np.sqrt(X**2 + Y**2)
    PHI_NON = np.arctan2(Y,X)
    #R,THETA, PHI = np.meshgrid(r,theta,phi)
    
    sphfunc = interp.griddata((np.ravel(R_NON),np.ravel(PHI_NON)),
		np.ravel(cartfunc),(r[None,:],phi[:,None]),
		method='cubic',fill_value=1e-30)

    return sphfunc

def cart2sph3D(cartfunc,x,y,z,r,theta,phi):
    print('Transforming from cartesian to spherical coordinates...')
    X,Y,Z = np.meshgrid(x,y,z)
    R_NON = np.sqrt(X**2 + Y**2 + Z**2)
    THETA_NON = np.arccos(Z / R_NON)
    PHI_NON = np.arctan2(Y,X)
    #R,THETA, PHI = np.meshgrid(r,theta,phi)
    
    sphfunc = interp.griddata((np.ravel(R_NON),np.ravel(THETA_NON),
		np.ravel(PHI_NON)),np.ravel(cartfunc),(r[None,None,:],theta[None,:,None],
		phi[:,None,None]),method='linear',fill_value=1e-30)

    return sphfunc

def get_radial_prob2D(sphfunc,rpts,theta,theta_weights):
    print('Getting radial probability...')
    # Allocate arrays 
    jacobian = np.zeros(sphfunc.shape)
    probr = np.zeros(sphfunc.shape[1])
    
    aux.create_sphjacobian(jacobian,rpts,np.ones(shape(theta)))
    aux.spherical_polar_integral(sphfunc,probr,jacobian,theta_weights)
    return probr

def get_radial_prob3D(sphfunc,rpts,theta,phi,
                      theta_weights):
    print('Getting radial probability...')
    # Allocate arrays
    jacobian = np.zeros(sphfunc.shape[1:3])
    probr = np.zeros(sphfunc.shape[2])

    weights = np.ones(np.shape(phi))[:,None] * theta_weights[None,:]
    aux.create_sphjacobian(jacobian,rpts,np.ones(shape(theta)))
    aux.spherical_angular_integral(sphfunc,probr,jacobian,weights)
    dphi = phi[1] - phi[0]
    probr[:] *= dphi
    return probr

def get_PAD(sphfunc,rpts,theta,phi):
    print('Getting PAD...')
    # Allocate arrays
    probPAD = np.zeros(sphfunc.shape[1:])
    
    aux.spherical_azimuthal_integral(sphfunc,probPAD,
                                     np.ones(shape(phi)))
    dphi = phi[1] - phi[0]
    probPAD[:] *= dphi
    return probPAD

def probk2probE(probk,ke,Eax,massfactor):
    print('From momentum to Energy probability...')
    probk2E = probk / ke
    Efromk = ke**2 * massfactor
    interpolant_probE = interp.interp1d(Efromk,probk2E,kind='cubic')
    probE = interpolant_probE(Eax)
    return probE

def extract_subset(bigmatrix,xaxis,yaxis,zaxis,limits):
	if(limits[0,0]>xaxis.all() or limits[1,0]<xaxis.all() or 
	limits[0,1]>yaxis.all() or limits[1,1]>yaxis.all() or 
	limits[0,2]>zaxis.all() or limits[1,2]>zaxis.all()):
		print('Limits out of bounds!')
		return
	
	# Locate the limits on the passes axis
	indxmin = np.argmin(np.abs(xaxis+limits[0,0]))
	indxmax = np.argmin(np.abs(xaxis-limits[1,0]))
	indymin = np.argmin(np.abs(yaxis+limits[0,1]))
	indymax = np.argmin(np.abs(yaxis-limits[1,1]))
	indzmin = np.argmin(np.abs(zaxis+limits[0,2]))
	indzmax = np.argmin(np.abs(zaxis-limits[1,2]))
	
	# Create axis for the subset
	subxaxis = xaxis[indxmin:indxmax]
	subyaxis = xaxis[indymin:indymax]
	subzaxis = xaxis[indzmin:indzmax]
	
	# Extract matrix
	submatrix = bigmatrix[indxmin:indxmax,indymin:indymax,indzmin:indzmax]
	
	return subxaxis,subyaxis,subzaxis,submatrix

def func_smoother(ax1,func,ax2):
    interpolant = interp.interp1d(ax1,func,kind='cubic')
    interpfunc = interpolant(ax2)
    return interpfunc

    
