#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import constants as const
import input_reader as inp
import axes as axes
import data as data
import coords as coords
import aux_funcs as aux_funcs
import cross as cross
import graph as graph


########################################
##       popsicle viewer.py
########################################

# Load all the parameters
inp_params = {}
inp.read_input(inp_params,'../input_popview.dat')

params = inp.parameters(inp_params)

### Create axes for simulation
print('Creating axis...')
# For coordinate space
ax = axes.axes(params.dke, params.kemax,
params.lmax,params.mmax)

# Load data from disk
if(params.draw_polar_amplitude):
    polar_prob, ntime = data.load_h5series(params.polar_filename)


# Get ke prob
if(params.draw_mes):
    mes = np.loadtxt(params.mes_filename)

# Get pes
if(params.draw_pes):
    pes = np.loadtxt(params.pes_filename)

#if(params.draw_diff_cross):
    # w = kax.ke_ax[np.argmax(probke)]**2*const.mufac
    # wki = w #kax.ke_ax[np.argmax(probke)]**2*const.mufac + params.Ip
    # FR = cross.get_fourier_afield(w,wki,params.pulse_duration)
    # diff_cross = cross.get_diff_cross(w,wki,probketheta,kax.ke_ax,params.pulse_duration,
    #                             params.A0)

###############################
### Plotting!!!
###############################

print('Plotting selected frames...')

if(params.draw_polar_amplitude):
    print('Plotting polar amplitude...')
    logplot = 1
    clamp = [-9,-6]
    xticks = None #np.arange(0,25,0.2)
    yticks = None #np.arange(-15.,16.,0.2)
    graph.oneframe_surf(polar_prob,ax.ke_ax,ax.theta_ax,clamp,
                        [ax.dke,2.8],[min(ax.theta_ax),max(ax.theta_ax)],
                        xticks,yticks,
                        r'$k_e$(a.u.)',r'$\theta$(rad)',
                        r'$\log_{10}|\Psi(k_e,\theta)|^2$',
                        logplot,params.makeframe,'polar_prob')

if(params.draw_mes):
    print('Plotting photoelectron momentum...')
    xticks = np.arange(0,3.0,0.4)
    yticks = None
    ylim = None #[0.0,max(probke)]
    logplot = 0
    graph.plotprob1d(pes,ax.ke_ax,
                     [0.0,2.8],ylim,
                     xticks,yticks,
                     'Momentum (au)',r'$P(k)$',logplot,
                     params.makeframe,'probke')

if(params.draw_pes):
    print('Plotting photoelectron energy spectra...')
    xticks = np.arange(0,3.0,0.4)
    yticks = None
    ylim = None #[0.0,max(probke)]
    logplot = 0
    graph.plotprob1d(pes,ax.Ee_ax,
                     [0.0,2.8],ylim,
                     xticks,yticks,
                     'Energy (au)',r'$P(E)$',logplot,
                     params.makeframe,'probEe')

# if(params.draw_Etotal):
#     xlim = [0.0,2.0]
#     ylim = None #[0.0,max(probE)]
#     logplot = 0
#     graph.plotprob1d(probE*1e4,kax.E_ax,xlim,ylim,
#                      xticks,yticks,
#                      'Energy (au)',r'$P(E) (10^{-4} au)$',
#                      logplot,params.makeframe,'probE')
#
# if(params.draw_diff_cross):
#     print('Plotting differential cross section...')
#     xticks = None #np.arange(kax.dke,1.0,0.2)
#     ylim = None #[0.0,max(probke)]
#     logplot = 0
#     graph.plotprob1d(diff_cross,kax.theta_ax,
#                      [0.0,np.pi],ylim,
#                      xticks,yticks,
#                      'theta (rad)',r'DCS (au)',logplot,
#                      params.makeframe,'diff_cross')
