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
#import numbapro as nb

class axes():
    def __init__(self,dke,kemax,lmax,mmax,fixed_nuclei=False,
    mufac=0.5,dkr=1.0,krmax=1.0,Mfac=1.0):

        self.dke = dke
        self.kemax = kemax
        self.lmax = lmax
        self.mmax = mmax
        self.mmin = -mmax

        self.maxkepts = int(kemax / dke)
        self.maxthetapts = 2 * self.lmax + 1
        self.maxphipts = 2 * self.mmax + 1
        self.dphi = 2.0 * np.pi / float(self.maxphipts)

        if(fixed_nuclei):
            self.dkr = draw_kr
            self.krmax = krmax
            self.maxkrpts = int(self.krmax / self.dkr)

        # Coordinates arrays
        #        self.maxptsEe = maxptsEe
        #        self.maxptsEn = maxptsEn
        #        self.maxptsEtotal = maxptsEtotal
        #        self.dEe = dEe
        #        self.dEn = dEn
        #        self.dE = dE
        #######################################################
        ########################################################

        if(fixed_nuclei):
            self.kr_ax = np.arange(self.dkr,self.krmax,self.dkr)

        # for i in range(maxptskr):
        #     self.kr_ax[i] = float(i+1) * self.dkr

        ######################################################

        self.ke_ax = np.arange(self.dke,self.kemax,self.dke)

        # for i in range(maxptske):
        #     self.ke_ax[i] = float(i+1) * self.dke

        ######################################################

        #for i in range(maxptstheta):
        #    self.theta_ax[i] = float(i+1) * self.dtheta
        self.costheta_ax , self.theta_weights = \
        np.polynomial.legendre.leggauss(self.maxthetapts)

        self.theta_ax = np.arccos(self.costheta_ax)

        ######################################################

        self.phi_ax = np.arange(self.dphi,2.0*np.pi,self.dphi)

        #for i in range(self.maxphipts):
        #    self.phi_ax[i] = float(i+1) * self.dphi

        ######################################################

        self.Ee_ax = self.ke_ax * self.ke_ax * mufac

        # for i in range(maxptsEn):
        #     self.En_ax[i] = float(i+1) * self.dEn
        #
        # for value in variable:
        #pass i in range(maxptsEe):
        #     self.Ee_ax[i] = float(i+1) * self.dEe
        #
        #
        # for value in variable:
        #pass value in variable:
        #pass i in range(maxptsEtotal):
        #     self.E_ax[i] = float(i+1) * self.dE
        #
        ########################################################
        ########################################################

        ###########################################################
