#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
##           axes.py
########################################
## This file contains the definition
## of the class axes.
########################################

# Modules a go-go!!
import numpy as np
import numba as nb

class axes():
    def __init__(self,dEe,Eemaxpts,lmax,mmax,fixed_nuclei=False,
                 dEn=1.0,Enmaxpts=1,mufac=0.5,Mfac=1.0,
                 dE=1.0,Emax=1.0):


        self.dEe      = dEe
        self.Eemaxpts = Eemaxpts
        self.dEn      = dEn
        self.Enmaxpts = Enmaxpts
        self.lmax     = lmax
        self.mmax     = mmax
        self.mmin     = -mmax

        self.dE       = dE
        self.Emax     = Emax
        self.dk       = np.sqrt(2.0 * self.dE)
        self.kmax     = np.sqrt(2.0 * self.Emax)

        self.Eemax = self.dEe * self.Eemaxpts
        self.Enmax = self.dEn * self.Enmaxpts        
        self.dke = np.sqrt(2.0 * self.dEe)
        self.kemax = np.sqrt(2.0 * self.Eemax)

        self.kemaxpts = Eemaxpts
        self.maxthetapts = self.lmax + 1
        self.maxphipts = self.mmax + 1
        self.dphi = 2.0 * np.pi / float(self.maxphipts + 1)
        self.Eemaxpts = int(self.Eemax / self.dEe)+ 1
        self.kmaxpts  = int(self.kmax / self.dk) + 1

        if(fixed_nuclei == False):
            self.dkn = np.sqrt(2.0 * dEn)
            self.knmax = np.sqrt(2.0 * Enmax)
            self.knmaxpts = Enmaxpts
            self.Emaxtotalpts = int(self.Emax / self.dE) + 1

        #######################################################
        ########################################################

        self.ke_ax = np.linspace(self.dke,self.kemax,self.kemaxpts)
        self.k_ax = np.linspace(self.dke,self.kmax,self.kmaxpts)

        ######################################################

        self.costheta_ax , self.theta_weights = \
        np.polynomial.legendre.leggauss(self.maxthetapts)

        self.theta_ax = np.arccos(self.costheta_ax)

        ######################################################

        self.phi_ax = np.linspace(0.0,2.0*np.pi,self.maxphipts)

        ######################################################

        self.Ee_ax = np.linspace(self.dke**2*mufac,self.Eemax,self.Eemaxpts)

        if(fixed_nuclei == False):
            self.kn_ax = np.linspace(self.dkn,self.knmax,self.knmaxpts)
            self.En_ax = np.linspace(0.0,self.Enmax,self.dEn)
            self.E_ax = np.linspace(0.0,self.Eemax,self.dE)

        ########################################################
        ########################################################

        ###########################################################

###########################################################

def gaussleg(x1,x2,degree):

    x,w = np.polynomial.legendre.leggauss(degree)

    xl = (x2-x1) * 0.5
    xm = (x2+x1) * 0.5
    
    xx = 0.5 * (x + 1) * (x2 - x1) + x1
    ww = xl * w

    return xx, ww
