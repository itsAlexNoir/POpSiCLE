#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
########################################
########################################
######### H2p moviemaker_hdf5 ##########
########################################
########################################


### Import modules
import h5py as h5py
import numpy as np
import scipy.interpolate as intp
#import numba as nb
import numbapro as nb
import matplotlib.pyplot as plt
import brewer2mpl as brew
import matplotlib.gridspec as gridspec
import matplotlib.patches as pat

plt.rcParams['font.family']='sans-serif'

########################################
####### Options for this script ########
########################################

 ## Load and calculate different data slices.
calc_krkrhokz = 0
calc_krhokz = 1
calc_krkz = 0
calc_PAD   = 1
calc_cross_section = 1

 ## Which frame do you want to plot?
draw_krkz = 0
draw_krhokz = 1
draw_EnEe = 0
draw_PAD = 0
draw_Etotal = 1
draw_cross_section = 1

 ## If you want to save the frame to a file
makeframe = 0
 ## Time of the final snapshot
frametime = 0
 ## If you want to draw over the plot 
drawover = 0

########################################
#### Parameters of the simulation ######
########################################

Nprocr =  1
Nprocrho = 4
Nprocz = 10

maxptsprockr = 1
maxptsprockrho = 90
maxptsprockz = 201

maxptsprockn = 1
maxptsprocke = 201
maxptsprocktot = 201

maxptskr = Nprocr * maxptsprockr 
maxptskrho = Nprocrho * maxptsprockrho
maxptskz = Nprocz * maxptsprockz

maxptskn = Nprocr * maxptsprockn 
maxptske = Nprocz * maxptsprocke
maxptsktot = Nprocz * maxptsprocktot

maxptsEn = 1
maxptsEe = 5000
maxptsEtotal = 5000

dkr = 1.0
dkrho = 0.0290
dkz = 0.01563759409452361 #0.0157

dkn = 1.0
dke = 0.01563759409452361 #0.0157
dtheta = 0.01563759409452361 #0.0157

dEn = 1.0
dEe = 0.001
dE = 0.001

krho_scaled = 0
krho_flat = 0

dt = 0.05
deltabench = 0.8
clamp_max = -6
clamp_min = - 12

maxkr = maxptskr * dkr
maxkrho = maxptskrho * dkrho
maxkz = (maxptskz - 1)/2 * dkz
maxkn = (maxptskn - 1)/2 * dkn
maxke = (maxptske - 1)/2 * dke

maxptstheta = int(np.pi / dtheta)

wavelength = 20.0 # nm
intensity  = 1.0e12 # W/cm2
no_cycles = 8
Eg = -0.595258373095276
Ip = 0.595258373095623 + 0.5 #1.1

########################################
############## Constants ###############
########################################

M1 = 1836.0
M2 = 1836.0
femto = 1e-15
autime = 2.418884326505e-17
aulength = 5.291e-11 # m
intensityau = 3.5e16 # W/cm2
speed_light = 137.0
fine_struc = 1.0 / speed_light

munuc = M1 * M2 / (M1 + M2)
muelec = ( M1 + M2 ) / (M1 + M2 + 1.0)

mufac = 0.5 / muelec
Mfac  = 0.5 / munuc

########################################

 # Let's convert the frequency to atomic units
wavelength0 = wavelength * 1e-9 / aulength
w0 = 2 * np.pi * 137 / wavelength0
period = 2.0 * np.pi / w0
pulse_duration = period * float(no_cycles)
pulse_bandwidth = 4.0 * np.pi / pulse_duration 
quiver = np.sqrt(intensity / intensityau) / w0 / w0
Up = (intensity / intensityau) / 4. / w0 / w0
E0 = np.sqrt(intensity / intensityau)
A0 = E0 / w0

########################################
############ Creating axes #############
########################################

print('Creating axis...')

kz_ax = np.zeros((maxptskz))
khpts = np.zeros((maxptskz))
khp = np.zeros((maxptskz))

kr_ax = np.zeros((maxptskr))
kfpts = np.zeros((maxptskr))
kfp = np.zeros((maxptskr))

krho_ax = np.zeros((maxptskrho))
kgpts = np.zeros((maxptskrho))
kgp = np.zeros((maxptskrho))

 # Scaled coordinates
sqrtkrho = np.zeros((maxptskrho))
sqrt1krho = np.zeros((maxptskrho))

if(draw_EnEe):
	kn_ax = np.zeros((maxptskn))
	knpts = np.zeros((maxptskn))
	knp = np.zeros((maxptskn))
	
if(draw_EnEe or calc_PAD):
	ke_ax = np.zeros((maxptske))
	kepts = np.zeros((maxptske))
	kep = np.zeros((maxptske))

if(calc_PAD or draw_EnEe):
	theta_ax = np.zeros((maxptstheta))
	En_ax = np.zeros((maxptsEn))
	Ee_ax = np.zeros((maxptsEe))
	
if(draw_Etotal):
	ktot_ax = np.zeros((maxptsktot))
	ktotpts = np.zeros((maxptsktot))
	ktotp = np.zeros((maxptsktot))
	E_ax = np.zeros((maxptsEtotal))


if(calc_krkz):
	jacobkrkz = np.ones((maxptskz,maxptskr))
if(calc_krhokz):
	jacobkrhokz = np.ones((maxptskz,maxptskrho))
#if(draw_knke):
#	jacobknke = np.zeros((maxptskn,maxptske))

######################################################

for i in range(maxptskr):
    kr_ax[i] = (i+1) * dkr

kfpts = kr_ax

kfp[:] = 1.0
    
######################################################

if(krho_scaled):
	if(krho_flat):
		krhoalpha = 1.0
		krhobeta = 1.0
	else:
		krhoalpha = 1.0
		krhobeta = 0.0
    
	for i in range(maxptskrho):
		krho_ax[i] = (i+1) * dkrho
	
	sqrtkrho = np.sqrt(krho_ax)
	sqrt1krho[:] = np.sqrt(krhoalpha + krhobeta * krho_ax[:])

	kgpts = krho_ax * sqrtkrho / sqrt1krho
	kgp[:] = 0.5 * sqrtkrho[:] *  (2.0 * krhobeta * krho_ax[:] + 3.0 * krhoalpha)

	kgp = kgp / (sqrt1krho ** 3)
else:
	for i in range(maxptskrho):
		krho_ax[i] = (i+1) * dkrho
	kgpts = krho_ax
	kgp[:] = 1.0
	
#######################################################

for i in range(maxptskz):
    kz_ax[i] = -maxkz + (i+1) * dkz


khpts = kz_ax
khp[:] = 1.0

#######################################################

#######################################################

if(draw_EnEe):
	for i in range(maxptskn):
		kn_ax[i] = (i+1) * dkn


	knpts = kn_ax
	knp[:] = 1.0

#######################################################

#######################################################

if(calc_PAD or draw_EnEe):
	for i in range(maxptske):
		ke_ax[i] = (i+1) * dke


	kepts = ke_ax
	kep[:] = 1.0

########################################################

if(calc_PAD or draw_EnEe):
	for i in range(maxptstheta):
		theta_ax[i] = (i+1) * dtheta

	for i in range(maxptsEn):
		En_ax[i] = (i+1) * dEn

	for i in range(maxptsEe):
		Ee_ax[i] = (i+1) * dEe


#########################################################

if(draw_Etotal):
	for i in range(maxptsEtotal):
		E_ax[i] = (i+1) * dE

############################

@nb.jit('void(float64[:],float64[:],float64[:],float64[:,:])',target='cpu')
def buildjacobian(a,b,c,jacobian):
    for j in range(b.shape[0]):
        for i in range(a.shape[0]):
            jacobian[i,j] = a[i] * b[j] * c[j]

############################

# if(calc_krkz):
# 	buildjacobian(khp,kfp,kfp,jacobkrkz)
# 	jacobkrkz = 1.0 / jacobkrkz
# if(calc_krhokz):
# 	buildjacobian(khp,kgpts,kgp,jacobkrhokz)
# 	jacobkrhokz = 1.0 / jacobkrhokz
#if(draw_knke):
#	buildjacobian(kep,kfp,kfp,jacobknke)
#	jacobknke = 1.0 / jacobknke


##########################################

 # Create meshgrid
if(calc_krkz):
	KR,KZ1 = np.meshgrid(kfpts,khpts)
if(calc_krhokz or draw_EnEe):
	KRHO,KZ2 = np.meshgrid(kgpts,khpts)
if(calc_PAD or draw_EnEe):
	EEpad,THETA = np.meshgrid(Ee_ax,theta_ax)
if(draw_EnEe):
	EE,EN = np.meshgrid(Ee_ax,En_ax)
	
########################################
########### Loading data ###############
########################################

print('Loading data...')

RedToBlue = brew.get_map('RdBu', 'Diverging', 9,reverse=True).mpl_colormap
blues = brew.get_map('YlGnBu', 'Sequential', 9,reverse=True).mpl_colormap
greys = brew.get_map('Greys', 'Sequential', 9,reverse=True).mpl_colormap
bmap = brew.get_map('Paired','Qualitative',10).mpl_colors


if(calc_krkrhokz):
	fkrkrhokz = h5py.File('probkrkrhokz.h5','r')
	probkrkz = np.zeros((maxptskz,maxptskrho,maxptskr))
	print('Loading selected data frame...')
	probkrkrhokz = fkrkrhokz[fkrkrhokz.keys()[0]]
	probkrkrhokz = np.fft.fftshift(probkrkrhokz,axes=0)
if(calc_krkz):
	fkrkz = h5py.File('probkrkz.h5','r')
	probkrkz = np.zeros((maxptskz,maxptskr))
	print('Loading selected data frame...')
	probkrkz = fkrkz[fkrkz.keys()[0]]
	#probkrkz *= dkrho
	probkrkz *= jacobkrkz * dkrho
	probkrkz = np.fft.fftshift(probkrkz,axes=0)
if(calc_krhokz):
	fkrhokz = h5py.File('probkrhokz.h5','r')
	probkrhokz = np.zeros((maxptskz,maxptskrho))
	print('Loading selected data frame...')
	probkrhokz = fkrhokz[fkrhokz.keys()[0]]
	#probkrhokz *= dkr
	probkrhokz *= jacobkrhokz * dkr
	probkrhokz = np.fft.fftshift(probkrhokz,axes=0)
	#probkrhokz = np.fft.fftshift(probkrhokz)
# if(draw_knke):
# 	fknke = h5py.File('probknke.h5','r')
# 	probknke = np.zeros((maxptskz,maxptskrho))
# 	print('Loading selected data frame...')
# 	probknke = fknke[fknke.keys()[0]]
#  	probknke *= jacobknke #* dkr
#  	probknke = np.fft.fftshift(probknke,axes=0)
# if(draw_ktotal and draw_Etotal):
# 	fktotal = h5py.File('probktotal.h5','r')
# 	probktotal = np.zeros((maxptsktot))
# 	print('Loading selected data frame...')
# 	probktotal = fktotal[fktotal.keys()[0]]
# 	probknke *= jacobknke * dkr
# 	probknke = np.fft.fftshift(probknke,axes=0)
# if(draw_Etotal):
# 	probEtotal = probktotal / ktotpts


#########################################
########### Calculate PAD ###############
#########################################

if(calc_PAD or draw_EnEe):
	
	print('Calculating PAD...')
	#KE_NON = np.sqrt(KRHO**2 + KZ2**2)
	E_NON = (KRHO**2 + KZ2**2) * mufac
	THETA_NON = np.arctan2(KRHO, KZ2)
	
	print('Interpolating grid to get PAD...')
	if(calc_krhokz):
		probpad2D = intp.griddata((np.ravel(E_NON),np.ravel(THETA_NON)), np.ravel(probkrhokz),
		(Ee_ax[None,:],theta_ax[:,None]),method='cubic',fill_value=1e-30)
	
	if(calc_krkrhokz):
		#probpad3D = np.zeros((maxptskr,maxptske,maxptstheta))
		for ikr in range(maxptskr):
			probpad3D[ikr,:,:] = intp.griddata((np.ravel(E_NON),np.ravel(THETA_NON)),
			 np.ravel(probkrkrhokz[ikr,:,:]),(Ee_ax[None,:],theta_ax[:,None]),
			 method='cubic',fill_value=1e-30)
	
###############################################################################
############################## Calculate EnEe #################################
###############################################################################

@nb.jit('void(float64[:],float64[:,:],float64[:])',target='cpu')
def integrateOverTheta2D(prob,probpad,jacob):
    for i in range(prob.shape[0]):
        for j in range(probpad.shape[0]):
        	for k in range(probpad.shape[1]):
	            prob[k] = prob[k] + probpad[j,k] * jacob[j]

@nb.jit('void(float64[:,:],float64[:,:,:])',target='cpu')
def integrateOverTheta3D(prob,probpad):
    for i in range(prob.shape[0]):
        for j in range(probpad.shape[0]):
        	for k in range(probpad.shape[1]):
	            prob[i,k] = prob[i,k] + probpad[i,j,k]

@nb.jit('void(float64[:,:],float64[:,:],float64[:],float64[:])',target='cpu')
def buildEnEe(probEnEe,probknke,kn,ke):
    for i in range(probEnEe.shape[0]):
        for j in range(probEnEe.shape[1]):
            probEnEe[i,j] = probknke[i,j] / kn[j] / ke[i]
            
#################################################################################

if(draw_EnEe):
	# Integrate over theta
	probknke = np.zeros((maxptske,maxptskn))
	
	if(calc_krkz):
		KE1_NON = np.sqrt(KZ1**2 + KZ1**2)
		
		probknke = intp.griddata((np.ravel(KE1_NON),np.ravel(KR)), np.ravel(probkrkz),
		(kepts[:,None],knpts[None,:]),method='cubic',fill_value=1e-30)
	
	# Transform to energies
	probEnEe = np.zeros((maxptske,maxptskn))
	buildEnEe(probEnEe,probknke,knpts,kepts)
	probEnEe[:,:] *= muelec * munuc

#################################################################################

if(draw_Etotal):
	
	if(calc_krhokz):
		probEtotal = np.zeros((maxptsEtotal))
		jacobketheta = np.sin(theta_ax)
		integrateOverTheta2D(probEtotal,probpad2D,jacobketheta)
		probEtotal[:] *= dtheta
		
		# Interpolate points for total energy
		probE = intp.interp1d(Ee_ax,probEtotal,kind='cubic')
		probEtotal = probE(E_ax)
		
	if(calc_krkz):
		ETOT_NON = EN_NON + EE_NON
		probE = intp.interp1d(np.ravel(ETOT_NON),np.ravel(probEnEe),kind='cubic')
		probEtotal = probE(E_ax)

	wkmax = E_ax[np.argmax(probEtotal)]

#################################################################################

 ### Calculate one-photon cross section
if(calc_cross_section):
	wi = Eg
	FRWsq = np.zeros((maxptsEtotal))
       	cross = np.zeros((maxptsEtotal))
	 ### The Fourier transform of the pulse.
	 ### (Analytical expression for sin2 shape)
	FRWsq = 4.0 * np.pi**4 * np.sin( ((w0+wi) - E_ax) * pulse_duration * 0.5)**2 / \
	    ((pulse_duration**2 * ((w0+wi) - E_ax )**2 - 4.0 * np.pi**2)**2 *((w0+wi) - E_ax )**2)
	# Finally, the cross section
					 
	cross = 4.0 * np.pi**2 * speed_light * probEtotal / (E_ax - wi) / FRWsq / A0 / A0
	#cross = probEtotal / FRWsq / A0 / A0
	 ### Remove function out of the bandwidth
	mm = np.where(E_ax < (w0+wi) - pulse_bandwidth/2.0)
	pp = np.where(E_ax > (w0+wi) + pulse_bandwidth/2.0)
	cross[mm] = 0.0
	cross[pp] = 0.0
	
####################################################################################
##########################  Movie of the evolution  ################################
####################################################################################

print('Plotting selected frames...')

###################################################################### 
######################### For vertical plots #########################
############################# First frame ############################
###################################################################### 	


### Set the window
if(draw_krkz or draw_krhokz):
	#fig = plt.figure(figsize=(20,12),facecolor='white')
	fig = plt.figure(figsize=(8,13),facecolor='white')
	#fig = plt.figure(figsize=(12,12),facecolor='white')

if(draw_krkz and draw_krhokz):
	gs = gridspec.GridSpec(1,2)
elif(draw_krkz):
		gs = gridspec.GridSpec(1,1)
elif(draw_krhokz):
		gs = gridspec.GridSpec(1,1)

# ### R vs. Z
if(draw_krkz):
	print('Plotting KR vs KZ...')
	ax2 = plt.subplot(gs[0,0])
	surf1 = ax2.pcolormesh(KR,KZ1,np.log10(probkrkz),vmax=clamp_max,vmin=clamp_min,shading='gouraud')
	cb1 = fig.colorbar(surf1,ax=ax2,use_gridspec=True)
	cb1.set_ticks(np.arange(clamp_max,clamp_min-1,-1))
	cb1.set_label(r'$\log_{10}|\Psi(kR,kz)|^2$',rotation=90,fontsize=20)
	#ax2.set_aspect('equal')
	ax2.axis('tight')
	#ax2.set_xticks(np.arange(0,max(r_ax),5.))
	ax2.set_xticks(np.arange(0,30,2.))
	ax2.set_yticks(np.arange(-30,30,0.5))
	#ax2.set_yticks(np.arange(-maxptsz*dz/2,maxptsz*dz/2+10,20.))
	#ax2.set_yticks(np.linspace(-160.,160.,11.,endpoint=True))
	ax2.tick_params(axis='both',labelsize=16)
	ax2.set_xlim([0,20])
	ax2.set_ylim([-4,4])
	ax2.set_xlabel(r'$k_R$(a.u.)',fontsize=20)
	ax2.set_ylabel(r'$k_z$(a.u.)',fontsize=20)
	#ax2.set_title()

	if(drawover):
		ax2.plot([0,0],[0,35],'r',lw=2)
		ax2.plot([0,-16.6],[0,35],'r',lw=2)
		ax2.plot([0,16.6],[0,35],'r',lw=2)
		ax2.add_artist(pat.Arc((0,0), 10, 10, angle=0, theta1=65, theta2=115,color='r',lw=2))
		ax2.text(-3,15,r'$25^{\circ}$',color='r',fontsize=18)
		ax2.text(1.5,15,r'$25^{\circ}$',color='r',fontsize=18)



#  ### RHO vs. Z
if(draw_krkz & draw_krhokz):
	ax3 = plt.subplot(gs[0,1])
elif(draw_krhokz & (draw_krkz==0)):
	ax3 = plt.subplot(gs[0,0])

if(draw_krhokz):
	print('Plotting KRHO vs KZ...')
	surf2 = ax3.pcolormesh(KRHO,KZ2,np.log10(probkrhokz),vmax=clamp_max,vmin=clamp_min,shading='gouraud')
	cb2 = fig.colorbar(surf2,ax=ax3,use_gridspec=True)
	cb2.set_ticks(np.arange(clamp_max,clamp_min-1,-1))
	cb2.set_label(r'$\log_{10}|\Psi(k\rho,kz)|^2$',rotation=90,fontsize=20)
	#ax3.set_aspect('equal')
	ax3.axis('tight')
	#ax3.set_xticks(np.arange(0,max(krho_ax),20.))
	ax3.set_xticks(np.arange(0.,25.,1.))
	#ax3.set_yticks(np.arange(-max(kz_ax),max(kz_ax),1.))
	ax3.set_yticks(np.arange(-15.,16.0,1.))
	ax3.tick_params(axis='both',labelsize=16)
	ax3.set_xlim([dkrho,4])
	ax3.set_ylim([-4,4])
	ax3.set_xlabel(r'$k\rho$(a.u.)',fontsize=20)
	ax3.set_ylabel(r'$k_z$(a.u.)',fontsize=20)
	#ax3.set_title()
 

if(draw_krkz or draw_krhokz): 
	gs.tight_layout(fig)
	fig.set_tight_layout(True)

	 # Save the frame
	if(makeframe):
		filename = 'kimage%0.4d.pdf' % (frametime)
		fig.savefig(filename,format='pdf')

	
#########################################################################
#########################################################################

### Set the window
if(draw_PAD):
	#fig = plt.figure(figsize=(20,12),facecolor='white')
	fig = plt.figure(figsize=(13,8),facecolor='white')
	#fig = plt.figure(figsize=(12,12),facecolor='white')

	# ### KE vs. THETA

	print('Plotting THETA vs KE...')
	gs = gridspec.GridSpec(1,1)
	ax2 = plt.subplot(gs[0,0])
	surf1 = ax2.pcolormesh(THETA,EEpad,np.log10(probpad2D),vmax=clamp_max,vmin=clamp_min,shading='gouraud')
	cb1 = fig.colorbar(surf1,ax=ax2,use_gridspec=True)
	cb1.set_ticks(np.arange(clamp_max,clamp_min-1,-1))
	cb1.set_label(r'$\log_{10}|\Psi(E_e,\theta)|^2$',rotation=90,fontsize=20)
	#ax2.set_aspect('equal')
	ax2.axis('tight')
	#ax2.set_xticks(np.arange(0,20,2.))
	#ax2.set_yticks(np.arange(0,20,0.5))
	#ax2.set_yticks(np.linspace(-160.,160.,11.,endpoint=True))
	ax2.tick_params(axis='both',labelsize=16)
	ax2.set_xlim([0,np.pi])
	ax2.set_ylim([dEe,4])
	ax2.set_xlabel(r'$\theta$ (rad)',fontsize=20)
	ax2.set_ylabel(r'$E_e$(a.u.)',fontsize=20)
	
# 	if(drawover):
# 		ax2.plot([0,0],[0,35],'r',lw=2)
# 		ax2.plot([0,-16.6],[0,35],'r',lw=2)
# 		ax2.plot([0,16.6],[0,35],'r',lw=2)
# 		ax2.add_artist(pat.Arc((0,0), 10, 10, angle=0, theta1=65, theta2=115,color='r',lw=2))
# 		ax2.text(-3,15,r'$25^{\circ}$',color='r',fontsize=18)
# 		ax2.text(1.5,15,r'$25^{\circ}$',color='r',fontsize=18)


	#gs.tight_layout(fig)
	fig.set_tight_layout(True)
	
	 # Save the frame
	if(makeframe):
		filename = 'PAD%0.4d.pdf' % (frametime)
		fig.savefig(filename,format='pdf')
	
#########################################################################################


### Set the window
if(draw_EnEe):
	EN_NON, EE_NON = np.meshgrid(knpts*knpts*Mfac,kepts*kepts*mufac)
	fig = plt.figure(figsize=(8,13),facecolor='white')
	#fig = plt.figure(figsize=(12,12),facecolor='white')

	# ### EN vs. EE

	print('Plotting KN vs KE...')
	gs = gridspec.GridSpec(1,1)
	ax2 = plt.subplot(gs[0,0])
	surf1 = ax2.pcolormesh(EN_NON,EE_NON,np.log10(probEnEe),vmax=clamp_max,vmin=clamp_min,shading='gouraud')
	cb1 = fig.colorbar(surf1,ax=ax2,use_gridspec=True)
	cb1.set_ticks(np.arange(clamp_max,clamp_min-1,-1))
	cb1.set_label(r'$\log_{10}|\Psi(E_N,E_e)|^2$',rotation=90,fontsize=20)
	#ax2.set_aspect('equal')
	ax2.axis('tight')
	ax2.set_xticks(np.arange(0,20,0.5))
	ax2.set_yticks(np.arange(0,20,0.5))
	ax2.tick_params(axis='both',labelsize=16)
	ax2.set_xlim(0,2)
	ax2.set_ylim(0,4)
	ax2.set_xlabel(r'$E_N$(a.u.)',fontsize=20)
	ax2.set_ylabel(r'$E_e$(a.u.)',fontsize=20)
	#ax2.set_title()

# 	if(drawover):
# 		ax2.plot([0,0],[0,35],'r',lw=2)
# 		ax2.plot([0,-16.6],[0,35],'r',lw=2)
# 		ax2.plot([0,16.6],[0,35],'r',lw=2)
# 		ax2.add_artist(pat.Arc((0,0), 10, 10, angle=0, theta1=65, theta2=115,color='r',lw=2))
# 		ax2.text(-3,15,r'$25^{\circ}$',color='r',fontsize=18)
# 		ax2.text(1.5,15,r'$25^{\circ}$',color='r',fontsize=18)


	#gs.tight_layout(fig)
	fig.set_tight_layout(True)
	
	# Save the frame
	if(makeframe):
		filename = 'jointEnEe%0.4d.pdf' % (frametime)
		fig.savefig(filename,format='pdf')

#############################################################################

if(draw_Etotal):
	fig = plt.figure(figsize=(12,10),facecolor='white')
	
	# ### Etotal

	print('Plotting E TOTAL...')
	gs = gridspec.GridSpec(1,1)
	ax2 = plt.subplot(gs[0,0])
	surf1 = ax2.plot(E_ax,(probEtotal),color=bmap[5],linewidth=2)
	ax2.axis('tight')
	#ax2.set_xticks(np.arange(0,max(r_ax),5.))
	#ax2.set_xticks(np.arange(0,45,10.))
	ax2.tick_params(axis='both',labelsize=16)
	ax2.set_xlim([0,4])
	#ax2.set_ylim(clamp_min,clamp_max)
	ax2.set_xlabel('Energy(a.u.)',fontsize=20)
	ax2.set_ylabel(r'$|\Psi(E)|^2$(a.u.)',fontsize=20)
	#ax2.set_title()
	
	#gs.tight_layout(fig)
	fig.set_tight_layout(True)
	
	# Save the frame
	if(makeframe):
		filename = 'Etotal%0.4d.pdf' % (frametime)
		fig.savefig(filename,format='pdf')

#############################################################################

if(draw_cross_section):
	fig = plt.figure(figsize=(12,10),facecolor='white')
	
	# ### Etotal

	print('Plotting Cross section...')
	gs = gridspec.GridSpec(1,1)
	ax2 = plt.subplot(gs[0,0])
	surf1 = ax2.plot(E_ax,(cross),color=bmap[5],linewidth=2)
	ax2.axis('tight')
	#ax2.set_xticks(np.arange(0,max(r_ax),5.))
	ax2.set_xticks(np.arange(0,10,0.2))
	ax2.tick_params(axis='both',labelsize=16)
	ax2.set_xlim([0,4])
	#ax2.set_ylim(clamp_min,clamp_max)
	ax2.set_xlabel('Energy (a.u.)',fontsize=20)
	ax2.set_ylabel(r'Cross Section (a.u.)',fontsize=20)
	#ax2.set_title()
	
	#gs.tight_layout(fig)
	fig.set_tight_layout(True)
	
	# Save the frame
	if(makeframe):
		filename = 'Cross_section%0.4d.pdf' % (frametime)
		fig.savefig(filename,format='pdf')

########################################################################
	
if(makeframe==0):
	plt.show()

print('Script finished!!!')
