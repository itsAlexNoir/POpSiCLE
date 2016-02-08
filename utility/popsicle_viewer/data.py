#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import h5py as h5py
import constants as const

########################################
##           data.py
########################################

print('Loading data...')

########################################
########### Loading data ###############
########################################

def load_h5matrix(matrix,filename):
    fset = h5py.File(filename+'.h5','r')
    if (np.shape(matrix)==np.shape(fset[fset.keys()[0]])):
        print('Loading data from '+filename+'...')
        dset = fset[fset.keys()[0]]
        matrix = dset[:,:]
    else:
        print('Loading data error. Shapes does not match.')

    return matrix

###########################################################

def load_h5series(filename,frametime=0,domain=None,keaxis=None,
                  thetaaxis=None,phiaxis=None,kraxis=None):

    fset = h5py.File(filename+'.h5','r')
    grp = fset[fset.keys()[0]]
    ntime = len(grp)
    print('Loading selected data frame in group...')
    matrix = grp[grp.keys()[frametime]]

    if(domain is None):
        subset = matrix[:,:]

        return subset, ntime
    else:
        if(len(matrix.shape)==domain.shape[1]):
            if(domain[0,0]>=min(keaxis) and
                   domain[1,0]<=max(keaxis)):
                indxx = np.where((keaxis>=domain[0,0]) &
                                 (keaxis<=domain[1,0]))
                minikeaxis = keaxis[indxke]
                sliceke = slice(min(indxke[0]),max(indxke[0])+1)
            else:
                print('Wrong domain in ke coordinate')
                return

            if(domain.shape[1]>=2):
                if(domain[0,1]>=min(thetaaxis) and
                   domain[1,1]<=max(thetaaxis)):
                    indxtheta = np.where((thetaaxis>=domain[0,1]) &
                                    (thetaaxis<=domain[1,1]))
                    minithetaaxis = thetaaxis[indxtheta]
                    slicetheta = slice(min(indxtheta[0]),max(indxtheta[0])+1)
                else:
                    print('Wrong domain in theta coordinate')
                    return

            if(domain.shape[1]>=3):
                if(domain[0,2]>=min(phiaxis) and
                   domain[1,2]<=max(phiaxis)):
                    indxphi = np.where((phiaxis>=domain[0,2]) &
                                 (phiaxis<=domain[1,2]))
                    miniphiaxis = phiaxis[indxphi]
                    slicephi = slice(min(indxphi[0]),max(indxphi[0])+1)
                else:
                    print('Wrong domain in phi coordinate')
                    return

            if(domain.shape[1]>=4):
                if(domain[0,3]>=min(kraxis) and
                   domain[1,3]<=max(kraxis)):
                    indxkr = np.where((kraxis>=domain[0,3]) &
                                 (kraxis<=domain[1,3]))
                    miniraxis = kraxis[indxkr]
                    slicekr = slice(min(indxkr[0]),max(indxkr[0])+1)
                else:
                    print('Wrong domain in kr coordinate')
                    return
        else:
            print('Shape does not match!')
            return

         # Select subset
        if(domain.shape[1]==1):
            index = tuple([sliceke])
            subset = np.squeeze(matrix[index])
            return subset,ntime,minikeaxis
        elif(domain.shape[1]==2):
            index = tuple([sliceke,sliceheta])
            subset = np.squeeze(matrix[index])
            return subset,ntime,minikeaxis,minithetaaxis
        elif(domain.shape[1]==3):
            index = tuple([sliceke,sliceheta,slicephi])
            subset = np.squeeze(matrix[index])
            return subset,ntime,minikeaxis,minithetaaxis,miniphiaxis
        elif(domain.shape[1]==4):
            index = tuple([sliceke,sliceheta,slicephi,slicekr])
            subset = np.squeeze(matrix[index])
            return subset,ntime,minikeaxis,minithetaaxis,miniphiaxis,minikraxis
        else:
            print('Too much dimensions for our matrix!')
            return

############################################################
