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
inp.read_input(inp_params,'./input_popview.dat')

params = inp.parameters(inp_params)

### Create axes for simulation
print('Creating axis...')
# For coordinate space
ax = axes.axes(params.dke, params.kemax,
               params.lmax,params.mmax,
               params.dk,params.kmax,
               params.dEe,params.Eemax)

# Load data from disk
if(params.draw_polar_amplitude):
    polar_prob, ntime = data.load_h5series(params.polar_filename)


# Get ke prob
if(params.draw_mes or params.draw_total_cross):
    params.mes_filename+='.dat'
    mes = np.loadtxt(params.mes_filename)

# Get pes
if(params.draw_pes):
    params.pes_filename+='.dat'
    pes = np.loadtxt(params.pes_filename)

if(params.draw_total_cross):
    #probk = coords.func_smoother(ax.ke_ax,mes,ax.k_ax)
    w   = kax.ke_ax**2 * params.mufac
    wki = w - params.Ip
    total_cross = cross.get_total_cross(probke,ax.ke_ax,wki,
                                        params.w0,
                                        params.pulse_duration,
                                        params.pulse_bandwidth,
                                        params.A0)

#if(params.draw_diff_cross):
    # w = kax.ke_ax[np.argmax(probke)]**2*const.mufac
    # wki = w #kax.ke_ax[np.argmax(probke)]**2*const.mufac + params.Ip
    # FR = cross.get_fourier_afield(w,wki,params.pulse_duration)
    # diff_cross = cross.get_diff_cross(w,wki,probketheta,kax.ke_ax,params.pulse_duration,
    #                             params.A0)

if(params.draw_sampling_pes):
    params.samp_pes_filename+='.dat'
    samp_pes = np.loadtxt(params.samp_pes_filename)
if(params.draw_sampling_pad):
    samp_pad, ntime = data.load_h5series(params.samp_pad_filename)
if(params.draw_wavetime):
    params.samp_wavetime_filename+='.dat'
    samp_wavetime = np.loadtxt(params.samp_wavetime_filename)
    
###############################
### Plotting!!!
###############################

print('Plotting selected frames...')

if(params.draw_polar_amplitude):
    print('Plotting polar amplitude...')
    logplot = 1
    clamp = [-6,-4]
    xticks = None #np.arange(0,25,0.2)
    yticks = None #np.arange(-15.,16.,0.2)
    graph.oneframe_surf(polar_prob,ax.ke_ax,ax.theta_ax,clamp,True,
                        [min(ax.theta_ax),max(ax.theta_ax)],[0.0,3.0],
                        xticks,yticks,
                        r'$\theta$(rad)',r'$k_e$(a.u.)',
                        r'$\log_{10}|\Psi(k_e,\theta)|^2$',
                        logplot,params.makeframe,'polar_prob',
                        params.showframe)

if(params.draw_sph_polar_amplitude):
    print('Plotting polar amplitude in spherical coordinates...')
    logplot = 1
    clamp = [-6,-4]
    xticks = None #np.arange(0,25,0.2)
    yticks = np.arange(0.,5.,0.6)
    graph.oneframe_polarsurf(polar_prob,np.fliplr([ax.theta_ax])[0],ax.ke_ax,
                             clamp,True,
                             [0.0,3.0],
                             xticks,yticks,
                             r'$\theta$ (degrees)',r'$k_e$(a.u.)',
                             r'$\log_{10}|\Psi(k_e,\theta)|^2$',
                             logplot,params.makeframe,'prob_sphpolarPAD',
                             params.showframe)

    
if(params.draw_mes):
    print('Plotting photoelectron momentum...')
    probk = coords.func_smoother(ax.ke_ax,mes[:,1],ax.k_ax)
    xticks = np.arange(0,3.0,0.4)
    yticks = None
    ylim = None #[0.0,max(probke)]
    logplot = 0
    graph.plotprob1d(probk,ax.k_ax,
                     [0.0,2.8],ylim,
                     xticks,yticks,
                     'Momentum (au)',r'$P(k)$',logplot,
                     params.makeframe,'probke',
                     params.showframe)

if(params.draw_pes):
    print('Plotting photoelectron energy spectra...')
    probE = coords.func_smoother(ax.ke_ax**2*params.mufac,pes[:,1],ax.Ee_ax)
    xticks = np.arange(0,3.0,0.4)
    yticks = None
    ylim = None #[0.0,max(probke)]
    logplot = 0
    graph.plotprob1d(probE,ax.Ee_ax,
                     [0.0,2.8],ylim,
                     xticks,yticks,
                     'Energy (au)',r'$P(E)$',logplot,
                     params.makeframe,'probEe',
                     params.showframe)

if(params.draw_wavetime):
    print('Plotting sampling photoelectron energy spectra...')
    xticks = np.arange(0,3.0,0.4)
    yticks = None
    ylim = None #[0.0,max(probke)]
    logplot = 0
    graph.plotprob1d(samp_wavetime[:,1],samp_wavetime[:,0]*const.autime_fs,
                     [0.0,2.8],ylim,
                     xticks,yticks,
                     'Time (au)',r'$P(x)$',logplot,
                     params.makeframe,'prob_wavetime',
                     params.showframe)

        
if(params.draw_sampling_pes):
    print('Plotting sampling photoelectron energy spectra...')
    xticks = np.arange(0,3.0,0.4)
    yticks = None
    ylim = None #[0.0,max(probke)]
    logplot = 0
    graph.plotprob1d(samp_pes[:,1],samp_pes[:,0],
                     [0.0,2.8],ylim,
                     xticks,yticks,
                     'Energy (au)',r'$P(E)$',logplot,
                     params.makeframe,'prob_SPES',
                     params.showframe)

    
if(params.draw_sampling_pad):
    print('Plotting sampling pad...')
    logplot = 1
    clamp = [-9,-6]
    xticks = None #np.arange(0,np.pi,0.4)
    yticks = None #np.arange(0.0,5.0,0.3)
    graph.oneframe_surf(samp_pad,ax.theta_ax,samp_pes[:,0],clamp,False,
                        [min(ax.theta_ax),max(ax.theta_ax)],[0.0,3.0],
                        xticks,yticks,
                        r'$\theta$(rad)',r'Energy (a.u.)',
                        r'$\log_{10}|\Psi(E_e,\theta)|^2$',
                        logplot,params.makeframe,'prob_SPAD',
                        params.showframe)

# if(params.draw_Etotal):
#     xlim = [0.0,2.0]
#     ylim = None #[0.0,max(probE)]
#     logplot = 0
#     graph.plotprob1d(probE*1e4,kax.E_ax,xlim,ylim,
#                      xticks,yticks,
#                      'Energy (au)',r'$P(E) (10^{-4} au)$',
#                      logplot,params.makeframe,'probE',
#                      params.showframe)

if(params.draw_total_cross):
    print('Plotting photoelectron momentum...')
    probk = coords.func_smoother(ax.ke_ax,mes[:,1],ax.k_ax)
    xticks = np.arange(0,3.0,0.4)
    yticks = None
    ylim = None #[0.0,max(probke)]
    logplot = 0
    graph.plotprob1d(total_cross,ax.k_ax**2*0.5,
                     [0.0,2.8],ylim,
                     xticks,yticks,
                     'Ejected electron energy (au)',
                     r'Total cross section (au)',logplot,
                     params.makeframe,'total_cross',
                     params.showframe)


# if(params.draw_diff_cross):
#     print('Plotting differential cross section...')
#     xticks = None #np.arange(kax.dke,1.0,0.2)
#     ylim = None #[0.0,max(probke)]
#     logplot = 0
#     graph.plotprob1d(diff_cross,kax.theta_ax,
#                      [0.0,np.pi],ylim,
#                      xticks,yticks,
#                      'theta (rad)',r'DCS (au)',logplot,
#                      params.makeframe,'diff_cross',
#                      params.showframe)
