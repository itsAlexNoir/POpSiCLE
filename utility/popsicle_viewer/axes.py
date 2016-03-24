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
    def __init__(self,dke,kemax,lmax,mmax,dk,kmax,
                 dEe,Eemax,fixed_nuclei=False,mufac=0.5,
                 dkn=1.0,knmax=1.0,Mfac=1.0,
                 dEn=1.0,Enmax=1.0,dE=1.0,Emax=1.0):

        self.dke = dke
        self.kemax = kemax
        self.lmax = lmax
        self.mmax = mmax
        self.mmin = -mmax
        self.dEe  = dEe
        self.Eemax = Eemax
        self.dk  = dk
        self.kmax = kmax
        self.dEn   = dEn
        self.Enmax = Enmax
        self.dE   = dE
        self.Emax  = Emax

        self.kemaxpts = int(kemax / dke) + 1
        self.maxthetapts = self.lmax + 1
        self.maxphipts = self.mmax + 1
        self.dphi = 2.0 * np.pi / float(self.maxphipts + 1)
        self.Eemaxpts = int(self.Eemax / self.dEe)+ 1
        self.kmaxpts  = int(self.kmax / self.dk) + 1

        if(fixed_nuclei == False):
            self.dkn = dkn
            self.knmax = knmax
            self.knmaxpts = int(self.knmax / self.dkn) + 1
            self.Enmaxpts = int(self.Enmax / self.dEn) + 1
            self.Emaxtotalpts = int(self.Emax / self.dE) + 1

        #######################################################
        ########################################################

        if(fixed_nuclei):
            self.kn_ax = np.linspace(0.0,self.knmax,self.knmaxpts)

        ######################################################

        self.ke_ax = np.linspace(0.0,self.kemax,self.kemaxpts)
        self.k_ax = np.linspace(0.0,self.kmax,self.kmaxpts)

        ######################################################

        self.costheta_ax , self.theta_weights = \
        np.polynomial.legendre.leggauss(self.maxthetapts)

        self.theta_ax = np.arccos(self.costheta_ax)

        ######################################################

        self.phi_ax = np.linspace(0.0,2.0*np.pi,self.maxphipts)

        ######################################################

        self.Ee_ax = np.linspace(0.0,self.Eemax,self.Eemaxpts)

        if(fixed_nuclei == False):
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
