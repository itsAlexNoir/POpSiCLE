#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
##           graph.py
########################################
## This file contains the definition
## of the class axes.
########################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
import palettable as brew
#import brewer2mpl as brew
import matplotlib.gridspec as gridspec
import matplotlib.patches as pat

from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import  mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import MaxNLocator

plt.rcParams['font.family']='sans-serif'
plt.ioff()

# Color stuff
gbmap = color.ListedColormap(brew.colorbrewer.sequential.GnBu_9.mpl_colors)
gbcol = brew.colorbrewer.sequential.GnBu_9.mpl_colors
greycol = brew.colorbrewer.sequential.Greys_6.mpl_colors
#quacol = brew.colorbrewer.qualitative.Paired_8.mpl_colors
quacol = brew.colorbrewer.qualitative.Dark2_8.mpl_colors
quacol2 = brew.colorbrewer.qualitative.Set3_3.mpl_colors

################################################################

################################################################

def plotprob1d(func,axes,xlim=None,ylim=None,
        xticks=None,yticks=None,
        xlabel=None,ylabel=None,
               logplot=0,makeframe=0,filename='none',
               showframe=True):

    if(logplot):
        func = np.log10(func)

    fig = plt.figure(figsize=(12,10),facecolor='white')
    gs = gridspec.GridSpec(1,1)
    ax2 = plt.subplot(gs[0,0])
    surf1 = ax2.plot(axes,func,color=quacol[2],linewidth=2)
    ax2.axis('tight')
    if(xticks is not None):
        ax2.set_xticks(np.arange(xlim[0],xlim[1],0.2))
    ax2.tick_params(axis='both',labelsize=16)
    if(xlim is not None):
        ax2.set_xlim([xlim[0],xlim[1]])
    if(ylim is not None):
        ax2.set_ylim([ylim[0],ylim[1]])
    if(xlabel is not None):
        ax2.set_xlabel(xlabel,fontsize=20)
    if(ylabel is not None):
        ax2.set_ylabel(ylabel,fontsize=20)
    fig.set_tight_layout(True)
    # Save the frame
    if(makeframe):
        fig.savefig(filename+'.pdf',format='pdf')
    if(showframe):
        plt.show()
        
################################################################

def oneframe_surf(matrix,ax1,ax2,clamp,rotated=False,
                  xlim=None,ylim=None,xticks=None,yticks=None,
                  xlabel=None,ylabel=None,barlabel=None,logplot=0,
                  makeframe=0,filename='none',showframe=True):

    AX1, AX2 = np.meshgrid(ax1,ax2)

    if(logplot):
        matrix = np.log10(matrix)

    fig = plt.figure(figsize=(10,8),facecolor='white')
    gs = gridspec.GridSpec(1,1)

    ax1 = plt.subplot(gs[0,0])
    if(rotated):
        surf1 = ax1.pcolormesh(AX2,AX1,matrix,vmax=clamp[1],vmin=clamp[0],
                           shading='gouraud',rasterized='True')
    else:    
        surf1 = ax1.pcolormesh(AX1,AX2,matrix,vmax=clamp[1],vmin=clamp[0],
                           shading='gouraud',rasterized='True')
    cb1 = fig.colorbar(surf1,ax=ax1,use_gridspec=True)
    cb1.set_ticks(np.arange(clamp[1],clamp[0]-1,-1))
    cb1.ax.tick_params(labelsize=18)
    if(barlabel is not None):
        cb1.set_label(barlabel,rotation=90,fontsize=20)
    #ax1.set_aspect('equal')
    ax1.axis('tight')
    if(xticks is not None):
        ax1.set_xticks(xticks)
    if(yticks is not None):
        ax1.set_yticks(yticks)
    ax1.tick_params(axis='both',labelsize=16)
    if(xlim is not None):
        ax1.set_xlim(xlim)
    if(ylim is not None):
        ax1.set_ylim(ylim)
    if(xlabel is not None):
        ax1.set_xlabel(xlabel,fontsize=20)
    if(ylabel is not None):
        ax1.set_ylabel(ylabel,fontsize=20)
    #ax2.set_title()
    if(makeframe):
        fig.savefig(filename+'.pdf',format='pdf',rasterized=True)
    if(showframe):
        plt.show()
        
################################################################

def twoframe_surf(matrix1,ax1,ax2,clamp1,matrix2,ax3,ax4,clamp2,
                  xlim1=None,ylim1=None,xticks1=None,yticks1=None,
                  xlabel1=None,ylabel1=None,barlabel1=None,
                  xlim2=None,ylim2=None,xticks2=None,yticks2=None,
                  xlabel2=None,ylabel2=None,barlabel2=None,
                  logplot=0,makeframe=0,filename='none',showframe=True):

    AX1, AX2 = np.meshgrid(ax1,ax2)
    AX3, AX4 = np.meshgrid(ax3,ax4)
    if(logplot):
        matrix1 = np.log10(matrix1)
        matrix2 = np.log10(matrix2)

    fig = plt.figure(figsize=(12,8),facecolor='white')
    gs = gridspec.GridSpec(1,2)

    ax1 = plt.subplot(gs[0,0])
    surf1 = ax1.pcolormesh(AX1,AX2,matrix1,vmax=clamp1[1],vmin=clamp1[0],
                           shading='gouraud',rasterized='True')
    cb1 = fig.colorbar(surf1,ax=ax1,use_gridspec=True)
    cb1.set_ticks(np.arange(clamp1[1],clamp1[0]-1,-1))
    cb1.ax.tick_params(labelsize=18)
    if(barlabel1 is not None):
        cb1.set_label(barlabel1,rotation=90,fontsize=20)
    ax1.set_aspect('equal')
    #ax1.axis('tight')
    if(xticks1 is not None):
        ax1.set_xticks(xticks1)
    if(yticks1 is not None):
        ax1.set_yticks(yticks1)
    ax1.tick_params(axis='both',labelsize=16)
    if(xlim1 is not None):
        ax1.set_xlim(xlim1)
    if(ylim1 is not None):
        ax1.set_ylim(ylim1)
    if(xlabel1 is not None):
        ax1.set_xlabel(xlabel1,fontsize=20)
    if(ylabel1 is not None):
        ax1.set_ylabel(ylabel1,fontsize=20)
    #ax2.set_title()

    ax2 = plt.subplot(gs[0,1])
    surf2 = ax1.pcolormesh(AX3,AX4,matrix2,vmax=clamp2[1],vmin=clamp2[0],
                           shading='gouraud',rasterized='True')
    cb2 = fig.colorbar(surf2,ax=ax2,use_gridspec=True)
    cb2.set_ticks(np.arange(clamp2[1],clamp2[0]-1,-1))
    cb2.ax.tick_params(labelsize=18)
    if(barlabel2 is not None):
        cb1.set_label(barlabel2,rotation=90,fontsize=20)
    ax2.set_aspect('equal')
    #ax1.axis('tight')
    if(xticks2 is not None):
        ax2.set_xticks(xticks2)
    if(yticks2 is not None):
        ax2.set_yticks(yticks2)
    ax2.tick_params(axis='both',labelsize=16)
    if(xlim2 is not None):
        ax2.set_xlim(xlim2)
    if(ylim2 is not None):
        ax2.set_ylim(ylim2)
    if(xlabel2 is not None):
        ax2.set_xlabel(xlabel2,fontsize=20)
    if(ylabel2 is not None):
        ax2.set_ylabel(ylabel2,fontsize=20)
    #ax2.set_title()
    if(makeframe):
        fig.savefig(filename+'.pdf',format='pdf',rasterized=True)
    if(showframe):
        plt.show()

################################################################

def threeframe_surf(matrix1,ax1,ax2,clamp1,matrix2,ax3,ax4,clamp2,
                    matrix3,ax5,ax6,clamp3,
                    xlim1=None,ylim1=None,xticks1=None,yticks1=None,
                    xlabel1=None,ylabel1=None,barlabel1=None,
                    xlim2=None,ylim2=None,xticks2=None,yticks2=None,
                    xlabel2=None,ylabel2=None,barlabel2=None,
                    xlim3=None,ylim3=None,xticks3=None,yticks3=None,
                    xlabel3=None,ylabel3=None,barlabel3=None,
                    logplot=0,makeframe=0,filename='none',showframe=True):

    AX1, AX2 = np.meshgrid(ax1,ax2)
    AX3, AX4 = np.meshgrid(ax3,ax4)
    AX5, AX6 = np.meshgrid(ax5,ax6)
    if(logplot):
        matrix1 = np.log10(matrix1)
        matrix2 = np.log10(matrix2)
        matrix3 = np.log10(matrix3)

    fig = plt.figure(figsize=(20,8),facecolor='white')
    gs = gridspec.GridSpec(1,3)

    ax1 = plt.subplot(gs[0,0])
    surf1 = ax1.pcolormesh(AX1,AX2,matrix1,vmax=clamp1[1],vmin=clamp1[0],
                           shading='gouraud',rasterized='True')
    cb1 = fig.colorbar(surf1,ax=ax1,use_gridspec=True)
    cb1.set_ticks(np.arange(clamp1[1],clamp1[0]-1,-1))
    cb1.ax.tick_params(labelsize=18)
    if(barlabel1 is not None):
        cb1.set_label(barlabel1,rotation=90,fontsize=20)
    ax1.set_aspect('equal')
    #ax1.axis('tight')
    if(xticks1 is not None):
        ax1.set_xticks(xticks1)
    if(yticks1 is not None):
        ax1.set_yticks(yticks1)
    ax1.tick_params(axis='both',labelsize=16)
    if(xlim1 is not None):
        ax1.set_xlim(xlim1)
    if(ylim1 is not None):
        ax1.set_ylim(ylim1)
    if(xlabel1 is not None):
        ax1.set_xlabel(xlabel1,fontsize=20)
    if(ylabel1 is not None):
        ax1.set_ylabel(ylabel1,fontsize=20)
    #ax2.set_title()

    ax2 = plt.subplot(gs[0,1])
    surf2 = ax2.pcolormesh(AX3,AX4,matrix2,vmax=clamp2[1],vmin=clamp2[0],
                           shading='gouraud',rasterized='True')
    cb2 = fig.colorbar(surf2,ax=ax2,use_gridspec=True)
    cb2.set_ticks(np.arange(clamp2[1],clamp2[0]-1,-1))
    cb2.ax.tick_params(labelsize=18)
    if(barlabel2 is not None):
        cb2.set_label(barlabel2,rotation=90,fontsize=20)
    ax2.set_aspect('equal')
    #ax1.axis('tight')
    if(xticks2 is not None):
        ax2.set_xticks(xticks2)
    if(yticks2 is not None):
        ax2.set_yticks(yticks2)
    ax2.tick_params(axis='both',labelsize=16)
    if(xlim2 is not None):
        ax2.set_xlim(xlim2)
    if(ylim2 is not None):
        ax2.set_ylim(ylim2)
    if(xlabel2 is not None):
        ax2.set_xlabel(xlabel2,fontsize=20)
    if(ylabel2 is not None):
        ax2.set_ylabel(ylabel2,fontsize=20)
    #ax2.set_title()


    ax3 = plt.subplot(gs[0,2])
    surf3 = ax3.pcolormesh(AX5,AX6,matrix3,vmax=clamp3[1],vmin=clamp3[0],
                           shading='gouraud',rasterized='True')
    cb3 = fig.colorbar(surf3,ax=ax3,use_gridspec=True)
    cb3.set_ticks(np.arange(clamp3[1],clamp3[0]-1,-1))
    cb3.ax.tick_params(labelsize=18)
    if(barlabel3 is not None):
        cb3.set_label(barlabel3,rotation=90,fontsize=20)
    ax3.set_aspect('equal')
    #ax1.axis('tight')
    if(xticks3 is not None):
        ax3.set_xticks(xticks3)
    if(yticks3 is not None):
        ax3.set_yticks(yticks3)
    ax3.tick_params(axis='both',labelsize=16)
    if(xlim3 is not None):
        ax3.set_xlim(xlim3)
    if(ylim3 is not None):
        ax3.set_ylim(ylim3)
    if(xlabel3 is not None):
        ax3.set_xlabel(xlabel3,fontsize=20)
    if(ylabel3 is not None):
        ax3.set_ylabel(ylabel3,fontsize=20)
    #ax3.set_title()

    if(makeframe):
        fig.savefig(filename+'.pdf',format='pdf',rasterized=True)
    if(showframe):
        plt.show()

################################################################

def fourframe_surf(matrix1,ax1,ax2,clamp1,matrix2,ax3,ax4,clamp2,
                   matrix3,ax5,ax6,clamp3,matrix4,ax7,ax8,clamp4,
                   xlim1=None,ylim1=None,xticks1=None,yticks1=None,
                   xlabel1=None,ylabel1=None,barlabel1=None,
                   xlim2=None,ylim2=None,xticks2=None,yticks2=None,
                   xlabel2=None,ylabel2=None,barlabel2=None,
                   xlim3=None,ylim3=None,xticks3=None,yticks3=None,
                   xlabel3=None,ylabel3=None,barlabel3=None,
                   xlim4=None,ylim4=None,xticks4=None,yticks4=None,
                   xlabel4=None,ylabel4=None,barlabel4=None,
                   logplot=0,makeframe=0,filename='none',showframe=True):

    AX1, AX2 = np.meshgrid(ax1,ax2)
    AX3, AX4 = np.meshgrid(ax3,ax4)
    AX5, AX6 = np.meshgrid(ax5,ax6)
    AX7, AX8 = np.meshgrid(ax7,ax8)
    if(logplot):
        matrix1 = np.log10(matrix1)
        matrix2 = np.log10(matrix2)
        matrix3 = np.log10(matrix3)
        matrix4 = np.log10(matrix4)

    fig = plt.figure(figsize=(20,10),facecolor='white')
    gs = gridspec.GridSpec(2,2)

    ax1 = plt.subplot(gs[0,0])
    surf1 = ax1.pcolormesh(AX1,AX2,matrix1,vmax=clamp1[1],vmin=clamp1[0],
                           shading='gouraud',rasterized='True')
    cb1 = fig.colorbar(surf1,ax=ax1,use_gridspec=True)
    cb1.set_ticks(np.arange(clamp1[1],clamp1[0]-1,-1))
    cb1.ax.tick_params(labelsize=18)
    if(barlabel1 is not None):
        cb1.set_label(barlabel1,rotation=90,fontsize=20)
    ax1.set_aspect('equal')
    #ax1.axis('tight')
    if(xticks1 is not None):
        ax1.set_xticks(xticks1)
    if(yticks1 is not None):
        ax1.set_yticks(yticks1)
    ax1.tick_params(axis='both',labelsize=16)
    if(xlim1 is not None):
        ax1.set_xlim(xlim1)
    if(ylim1 is not None):
        ax1.set_ylim(ylim1)
    if(xlabel1 is not None):
        ax1.set_xlabel(xlabel1,fontsize=20)
    if(ylabel1 is not None):
        ax1.set_ylabel(ylabel1,fontsize=20)
    #ax2.set_title()

    ax2 = plt.subplot(gs[0,1])
    surf2 = ax2.pcolormesh(AX3,AX4,matrix2,vmax=clamp2[1],vmin=clamp2[0],
                           shading='gouraud',rasterized='True')
    cb2 = fig.colorbar(surf2,ax=ax2,use_gridspec=True)
    cb2.set_ticks(np.arange(clamp2[1],clamp2[0]-1,-1))
    cb2.ax.tick_params(labelsize=18)    
    if(barlabel2 is not None):
        cb2.set_label(barlabel2,rotation=90,fontsize=20)
    ax2.set_aspect('equal')
    #ax2.axis('tight')
    if(xticks2 is not None):
        ax2.set_xticks(xticks2)
    if(yticks2 is not None):
        ax2.set_yticks(yticks2)
    ax2.tick_params(axis='both',labelsize=16)
    if(xlim2 is not None):
        ax2.set_xlim(xlim2)
    if(ylim2 is not None):
        ax2.set_ylim(ylim2)
    if(xlabel2 is not None):
        ax2.set_xlabel(xlabel2,fontsize=20)
    if(ylabel2 is not None):
        ax2.set_ylabel(ylabel2,fontsize=20)
    #ax2.set_title()


    ax3 = plt.subplot(gs[1,0])
    surf3 = ax3.pcolormesh(AX5,AX6,matrix3,vmax=clamp3[1],vmin=clamp3[0],
                           shading='gouraud',rasterized='True')
    cb3 = fig.colorbar(surf3,ax=ax3,use_gridspec=True)
    cb3.set_ticks(np.arange(clamp3[1],clamp3[0]-1,-1))
    cb3.ax.tick_params(labelsize=18)
    if(barlabel3 is not None):
        cb3.set_label(barlabel3,rotation=90,fontsize=20)
    ax3.set_aspect('equal')
    #ax3.axis('tight')
    if(xticks3 is not None):
        ax3.set_xticks(xticks3)
    if(yticks3 is not None):
        ax3.set_yticks(yticks3)
    ax3.tick_params(axis='both',labelsize=16)
    if(xlim3 is not None):
        ax3.set_xlim(xlim3)
    if(ylim3 is not None):
        ax3.set_ylim(ylim3)
    if(xlabel3 is not None):
        ax3.set_xlabel(xlabel3,fontsize=20)
    if(ylabel3 is not None):
        ax3.set_ylabel(ylabel3,fontsize=20)
    #ax3.set_title()

    ax4 = plt.subplot(gs[1,1])
    surf4 = ax4.pcolormesh(AX8,AX7,matrix4,vmax=clamp4[1],vmin=clamp4[0],
                           shading='gouraud',rasterized='True')
    cb4 = fig.colorbar(surf4,ax=ax4,use_gridspec=True)
    cb4.set_ticks(np.arange(clamp4[1],clamp4[0]-1,-1))
    cb4.ax.tick_params(labelsize=18)
    if(barlabel4 is not None):
        cb4.set_label(barlabel4,rotation=90,fontsize=20)
    #ax4.set_aspect('equal')
    ax4.axis('tight')
    if(xticks4 is not None):
        ax4.set_xticks(xticks4)
    if(yticks4 is not None):
        ax4.set_yticks(yticks4)
    ax4.tick_params(axis='both',labelsize=16)
    if(xlim4 is not None):
        ax4.set_xlim(xlim4)
    if(ylim4 is not None):
        ax4.set_ylim(ylim4)
    if(xlabel3 is not None):
        ax4.set_xlabel(xlabel4,fontsize=20)
    if(ylabel4 is not None):
        ax4.set_ylabel(ylabel4,fontsize=20)
    #ax4.set_title()

    if(makeframe):
        fig.savefig(filename+'.pdf',format='pdf',rasterized=True)
    if(showframe):
        plt.show()

################################################################

def oneframe_polarsurf(matrix,ax1,ax2,clamp,axes_aspect=True,
                       ylim=None,xticks=None,yticks=None,
                       xlabel=None,ylabel=None,barlabel=None,logplot=0,
                       makeframe=0,filename='none',showframe=True,Up=0.0):

    matrix_twice = np.concatenate((matrix,matrix[::-1,:],matrix[0:1,:]),axis=0)
    ax1_twice = np.concatenate((ax1,2.0*np.pi - ax1[::-1],ax1[0:1]),axis=0)
    
    AX1, AX2 = np.meshgrid(ax2,ax1_twice)
    if(logplot):
        matrix_twice = np.log10(matrix_twice)

    fig = plt.figure(figsize=(10,8),facecolor='white')
    #fig = plt.figure(figsize=(8,13),facecolor='white')
    gs = gridspec.GridSpec(1,1)
    
    ax1 = plt.subplot(gs[0,0],polar=True)
    surf1 = ax1.pcolormesh(AX2,AX1,matrix_twice,vmax=clamp[1],vmin=clamp[0],
                               shading='gouraud',rasterized='True')
    ax1.plot(np.linspace(0, 2*np.pi, 100), np.ones(100)*2.0*Up,
            lw=2,color='white', linestyle='-')
    ax1.plot(np.linspace(0, 2*np.pi, 100), np.ones(100)*10.0*Up,
            lw=2,color='white', linestyle='-')
    #ax1.text(np.pi/1.9, 2.8*Up, r'$2U_p$', fontsize=20,color='white')
    ax1.text(np.pi/1.9, 4.5*Up, r'$2U_p$', fontsize=20,color='white')
    ax1.text(np.pi/2.0, 9.2*Up, r'$10U_p$', fontsize=20,color='white')

    cb1 = fig.colorbar(surf1,ax=ax1,use_gridspec=True,pad=0.1)
    cb1.set_ticks(np.arange(clamp[1],clamp[0]-1,-1))
    cb1.ax.tick_params(labelsize=18)
    if(barlabel is not None):
        cb1.set_label(barlabel,rotation=90,fontsize=20)
    if(axes_aspect is True):
        ax1.set_aspect('equal')
    #ax1.axis('tight')
    ax1.set_theta_offset(np.pi/2.0)
    if(xticks is not None):
        ax1.set_xticks(xticks)
    ax1.set_rlabel_position(0)
    if(yticks is not None):        
        ax1.set_yticks([])
        #ax1.set_yticks(yticks)
    ax1.tick_params(axis='both',labelsize=16)
    ax1.tick_params(axis='y',colors='white')
    if(ylim is not None):
        ax1.set_ylim(ylim)
    if(xlabel is not None):
        ax1.set_xlabel(xlabel,fontsize=20)
    #if(ylabel is not None):
    #    ax1.set_ylabel(ylabel,fontsize=20)
    if(makeframe):
        fig.savefig(filename+'.pdf',format='pdf',rasterized=True)
    if(showframe):
        plt.show()

################################################################

def setup_axes3(fig, rect,r,theta,rticks,rlabel,thetalabel):
    """
    Sometimes, things like axis_direction need to be adjusted.
    """

    # rotate a bit for better orientation
    tr_rotate = Affine2D().translate(-90, 0)

    # scale degree to radians
    tr_scale = Affine2D().scale(-np.pi/180., 1.)

    tr = tr_rotate + tr_scale + PolarAxes.PolarTransform()

    #grid_locator1 = angle_helper.LocatorHMS(4)
    tick_formatter1 = angle_helper.FormatterHMS()

    numticks = len(rticks)
    grid_locator2 = MaxNLocator(numticks)

    ra0, ra1 = 0, 180
    cz0, cz1 = 0, max(r)
    grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                        extremes=(ra0, ra1, cz0, cz1),
                                        #grid_locator1=grid_locator1,
                                        grid_locator2=grid_locator2,
                                        #tick_formatter1=tick_formatter1,
                                        #tick_formatter2=None,
                                       )

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    # adjust axis
    ax1.axis["right"].set_axis_direction("bottom")
    ax1.axis["left"].set_axis_direction("top")
    ax1.axis["right"].toggle(ticklabels=True, label=True)
    ax1.axis["left"].toggle(ticklabels=True, label=True)
    ax1.axis["left"].major_ticklabels.set_axis_direction("left")
    ax1.axis["right"].major_ticklabels.set_axis_direction("left")    
    ax1.axis["right"].major_ticklabels.set_size(18)
    ax1.axis["left"].major_ticklabels.set_size(18)
    
    ax1.axis["bottom"].set_visible(False)
    ax1.axis["top"].set_axis_direction("top")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("bottom")
    ax1.axis["top"].label.set_axis_direction("bottom")
    ax1.axis["top"].major_ticklabels.set_size(18)
    
    
    #ax1.axis["right"].label.set_fontsize(18)
    ax1.axis["left"].label.set_fontsize(18)
    ax1.axis["top"].label.set_fontsize(18)
    
    ax1.axis["left"].label.set_text(rlabel)
    ax1.axis["top"].label.set_text(thetalabel)

    # create a parasite axes whose transData in RA, cz
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch # for aux_ax to have a clip path as in ax
    ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                        # drawn twice, and possibly over some other
                        # artists. So, we decrease the zorder a bit to
                        # prevent this.

    return ax1, aux_ax

################################################################

def setup_pcolormesh_semipolar(aux_ax, aux_X, aux_Y, Z, **kwargs):
    ax = aux_ax._parent_axes
    shape_orig = aux_X.shape
    TR = np.array([aux_X, aux_Y]).reshape((2, -1)).transpose() # coordinates in aux_ax
    XY = aux_ax.transAux.transform(TR) # coordinates in ax
    X, Y = XY.transpose().reshape((2, shape_orig[0], shape_orig[1]))

    CS = ax.pcolormesh(X, Y, Z, **kwargs)


    return CS


################################################################

def semipolar_surf(matrix,ax1,ax2,clamp,
                       ylim=None,xticks=None,yticks=None,
                       xlabel=None,ylabel=None,barlabel=None,logplot=0,
                       makeframe=0,filename='none',showframe=True):

    ax1 = ax1 * 180. / np.pi
    maxindx = np.argmin(np.abs(ax2-ylim[1])) + 1
    ax2_red = ax2[:maxindx]
    AX1, AX2 = np.meshgrid(ax1,ax2_red)
    #AX1, AX2 = np.meshgrid(ax2_red,ax1)

    matrix = matrix[:,:maxindx]

    if(logplot):
        matrix = np.log10(matrix)
        
    matrix = np.transpose(matrix)

    # Create frame
    fig = plt.figure(figsize=(10,12),facecolor='white')

    # Create axis
    ax3, aux_ax3 = setup_axes3(fig,111,ax2_red,ax1,yticks,
                               ylabel, xlabel)

    # Create actual plot
    surf1 = setup_pcolormesh_semipolar(aux_ax3,AX1,AX2,matrix,
                         vmax=clamp[1],vmin=clamp[0],
                         shading='gouraud',rasterized='True')

    cb1 = fig.colorbar(surf1,ax=ax3,use_gridspec=True,pad=0.15)
    cb1.set_ticks(np.arange(clamp[1],clamp[0]-1,-1))
    cb1.ax.tick_params(labelsize=18)
    if(barlabel is not None):
        cb1.set_label(barlabel,rotation=90,fontsize=20)

    if(makeframe):
        fig.savefig(filename+'.pdf',format='pdf',rasterized=True)
    if(showframe):
        plt.show()


################################################################
