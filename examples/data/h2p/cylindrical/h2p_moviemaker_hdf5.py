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
#import numba as nb
import numbapro as nb
import matplotlib.pyplot as plt
import brewer2mpl as brew
import matplotlib.gridspec as gridspec
import matplotlib.patches as pat
import matplotlib.animation as animation

plt.rcParams['font.family']='sans-serif'

########################################
####### Options for this script ########
########################################

 # Which frame do you want to plot?
draw_rz = 0
draw_rhoz = 1

 # If you want a movie
makemovie = 0       
 # If you want to save the frame to a file
makeframe = 0
 # Select the time step you want to plot
frametime = 0
# If you want to draw over the plot 
drawover = 0

########################################
#### Parameters of the simulation ######
########################################

Nprocr =  1
Nprocrho = 4
Nprocz = 10

maxptsprocr = 1
maxptsprocrho = 90
maxptsprocz = 201

maxptsr = Nprocr * maxptsprocr 
maxptsrho = Nprocrho * maxptsprocrho
maxptsz = Nprocz * maxptsprocz


dr = 1.0
drho = 0.1
dz = 0.2

rho_scaled = 1

dt = 0.05
deltabench = 0.2
graphicstime = 0.15
clamp_val = - 15

maxr = maxptsr * dr
maxrho = maxptsrho * drho
maxz = (maxptsz - 1)/2 * dz


########################################
############## Constants ###############
########################################

laserfrequency = 3e8 / 780.0e-9
femto = 1e-15
autime = 2.418884326505e-17

########################################
############ Creating axes #############
########################################

print('Creating axis...')

rho_ax = np.zeros((maxptsrho))
r_ax = np.zeros((maxptsr))
z_ax = np.zeros((maxptsz))

 # Scaled coordinates
gpts = np.zeros((maxptsrho))
fpts = np.zeros((maxptsr))
hpts = np.zeros((maxptsz))

gp = np.zeros((maxptsrho))
fp = np.zeros((maxptsr))
hp = np.zeros((maxptsz))

sqrtrho = np.zeros((maxptsrho))
sqrt1rho = np.zeros((maxptsrho))

jacobrz = np.zeros((maxptsz,maxptsr))
jacobrhoz = np.zeros((maxptsz,maxptsrho))

###########################

for i in range(maxptsr):
    r_ax[i] = (i+1) * dr

fpts = r_ax

fp[:] = 1.0
    
###########################

if(rho_scaled):
	rhoalpha = 1.0
	rhobeta = 0.0
else:
	rhoalpha = 1.0
	rhobeta = 1.0
    
for i in range(maxptsrho):
    rho_ax[i] = (i+1) * drho
    
sqrtrho = np.sqrt(rho_ax)
sqrt1rho[:] = np.sqrt(rhoalpha + rhobeta * rho_ax[:])

gpts = rho_ax * sqrtrho / sqrt1rho
gp[:] = 0.5 * sqrtrho[:] *  (2.0 * rhobeta * rho_ax[:] + 3.0 * rhoalpha)

gp = gp / (sqrt1rho ** 3)

###########################

for i in range(maxptsz):
    z_ax[i] = -maxz + (i+1) * dz

hpts = z_ax

hp[:] = 1.0

###########################

@nb.jit('void(float64[:],float64[:],float64[:],float64[:,:])',target='cpu')
def buildjacobian(a,b,c,jacobian):
    for j in range(b.shape[0]):
        for i in range(a.shape[0]):
            jacobian[i,j] = a[i] * b[j] * c[j]


buildjacobian(hp,fp,fp,jacobrz)
jacobrz = 1.0 / jacobrz
buildjacobian(hp,gpts,gp,jacobrhoz)
jacobrhoz = 1.0 / jacobrhoz

########################################
########### Loading data ###############
########################################

print('Loading data...')

RedToBlue = brew.get_map('RdBu', 'Diverging', 9,reverse=True).mpl_colormap
blues = brew.get_map('YlGnBu', 'Sequential', 9,reverse=True).mpl_colormap
greys = brew.get_map('Greys', 'Sequential', 9,reverse=True).mpl_colormap

time_data = np.loadtxt('./sim_h2plus_ce/field.dat')
numtime = len(time_data)
time_data[:,0] = time_data[:,0] * autime / femto

if(draw_rz):
	frz = h5py.File('probrz.h5','r')
	grp_rz = frz[frz.keys()[0]]
	ntime = len(grp_rz)
	probrz = np.zeros((maxptsz,maxptsr))

if(draw_rhoz):
	frhoz = h5py.File('probrhoz.h5','r')
	grp_rhoz = frhoz[frhoz.keys()[0]]
	if(draw_rz==0):
            ntime = len(grp_rhoz)
	probrhoz = np.zeros((maxptsz,maxptsrho))


graphicstep = int(graphicstime / deltabench)

time = np.zeros((ntime))
efield = np.zeros((ntime,2))

# Build temporal axis of the movie
for it in range(ntime):
    itime = it * graphicstep
    time[it] = time_data[itime,0]
    efield[it,0] = time_data[itime,1]
    efield[it,1] = time_data[itime,2]

################################################
# 

if(makemovie):
    print('Loading first data frame...')
    if(draw_rz):
    	probrz = grp_rz[grp_rz.keys()[0]]
    	probrz *= jacobrz * drho
    if(draw_rhoz):
    	probrhoz = grp_rhoz[grp_rhoz.keys()[0]]
    	probrhoz *= jacobrhoz * dr
   	
else:
    print('Loading selected data frame...')
    if(draw_rz):
    	probrz = grp_rz[grp_rz.keys()[frametime]]
    	probrz *= jacobrz * drho
    if(draw_rhoz):
    	probrhoz = grp_rhoz[grp_rhoz.keys()[frametime]]
    	probrhoz *= jacobrhoz * dr

#########################################################################
######################  Movie of the evolution  #########################
#########################################################################

if (makemovie):
	print('Starting movie...')
	

###################################################################### 
######################### For vertical plots #########################
############################# First frame ############################
###################################################################### 	

 # Create meshgrid
R,Z1 = np.meshgrid(r_ax,z_ax)
RHO,Z2 = np.meshgrid(gpts,z_ax)

### Set the window
#fig = plt.figure(figsize=(20,12),facecolor='white')
#fig = plt.figure(figsize=(15,12),facecolor='white')
#fig = plt.figure(figsize=(12,12),facecolor='white')
fig = plt.figure(figsize=(7,16),facecolor='white')
if(draw_rz & draw_rhoz):
	gs = gridspec.GridSpec(2,2,height_ratios=[1,8])
elif(draw_rz):
	gs = gridspec.GridSpec(2,1,height_ratios=[1,8])
elif(draw_rhoz):
	gs = gridspec.GridSpec(2,1,height_ratios=[1,8])
	
### Laser pulse
ax1 = plt.subplot(gs[0,:])
pulse1 = ax1.plot(time_data[:,0],time_data[:,2],'r',linewidth=2)
ax1.axis('tight')
#ax1.set_xlabel('time (fs)',fontsize=20)
ax1.set_ylabel('E field',fontsize=20)
#time_text = ax1.set_title('Time: ',fontsize=20)
time_text = ax1.set_title('Time: %3.2f fs' %(time[frametime]),fontsize=20)
ax1.set_xticks(np.arange(min(time_data[:,0]),max(time_data[:,0]),5.))
ax1.set_yticks([min(time_data[:,2]),0.0,max(time_data[:,2])])
#ax1.set_yticks([-0.03,0.0,0.03])
ax1.tick_params(axis='both',labelsize=16)
#pulse2, = ax1.plot([],[],'r',marker='o',markersize=10)
pulse2, = ax1.plot(time[frametime],efield[frametime,1],'r',marker='o',markersize=10)

### R vs. Z
if(draw_rz):
	ax2 = plt.subplot(gs[1,0])
	surf1 = ax2.pcolormesh(R,Z1,np.log10(probrz),vmax=0,vmin=clamp_val,shading='gouraud',rasterized=True)
	cb1 = fig.colorbar(surf1,ax=ax2,use_gridspec=True)
	cb1.set_ticks(np.arange(0,clamp_val-1,-1))
	cb1.set_label(r'$\log_{10}|\Psi(R,z)|^2$',rotation=90,fontsize=20)
	#ax2.set_aspect('equal')
	ax2.axis('tight')
	#ax2.set_xticks(np.arange(0,max(r_ax),10.))
	ax2.set_xticks(np.arange(0.,50.,10.))
	#ax2.set_yticks(np.arange(-200,210,20.))
	ax2.set_yticks(np.linspace(-50.,50.,11.,endpoint=True))
	ax2.tick_params(axis='both',labelsize=16)
	#ax2.set_xlim([0,32])
	#ax2.set_xlim([0,40])
	#ax2.set_ylim([-30,30])
	#ax2.set_ylim([-55,55])
	#ax2.set_ylim([-40,40])
	#ax2.set_ylim([-50,50])
	ax2.set_xlabel('R(a.u.)',fontsize=20)
	ax2.set_ylabel('z(a.u.)',fontsize=20)
	#ax2.set_title()

	if(drawover):
		ax2.plot([0,0],[0,35],'r',lw=2)
		ax2.plot([0,-16.6],[0,35],'r',lw=2)
		ax2.plot([0,16.6],[0,35],'r',lw=2)
		ax2.add_artist(pat.Arc((0,0), 10, 10, angle=0, theta1=65, theta2=115,color='r',lw=2))
		ax2.text(-3,15,r'$25^{\circ}$',color='r',fontsize=18)
		ax2.text(1.5,15,r'$25^{\circ}$',color='r',fontsize=18)



#  ### RHO vs. Z
if(draw_rz & draw_rhoz):
	ax3 = plt.subplot(gs[1,1])
elif(draw_rhoz & (draw_rz==0)):
	ax3 = plt.subplot(gs[1,0])
	
if(draw_rhoz):
	surf2 = ax3.pcolormesh(RHO,Z2,np.log10(probrhoz),vmax=0,vmin=clamp_val,shading='gouraud',rasterized=True)
	cb2 = fig.colorbar(surf2,ax=ax3,use_gridspec=True)
	cb2.set_ticks(np.arange(0,clamp_val-1,-1))
	cb2.set_label(r'$\log_{10}|\Psi(\rho,z)|^2$',rotation=90,fontsize=20)
	ax3.set_aspect('equal')
	#ax3.axis('tight')
	#ax3.set_xticks(np.arange(0,max(gpts),20.))
	#ax3.set_xticks(np.arange(0,70,10.))
	#ax3.set_yticks(np.arange(-max(z_ax),max(z_ax),10.))
	#ax3.set_yticks(np.arange(-200,210,20.))
	#ax3.set_yticks(np.linspace(-50.,50.,11.,endpoint=True))
	ax3.tick_params(axis='both',labelsize=16)
	#ax3.set_xlim([0,30])
	#ax3.set_xlim([0,50])
	#ax3.set_xlim([0,40])
	#ax3.set_xlim([0,100])
	#ax3.set_ylim([-55,55])
	#ax3.set_ylim([-40,40])
	#ax3.set_ylim([-50,50])
	ax3.set_xlabel(r'$\rho$(a.u.)',fontsize=20)
	ax3.set_ylabel('z(a.u.)',fontsize=20)
	#ax3.set_title()
 
 
gs.tight_layout(fig)

 # Save the frame
if(makeframe):
    filename = 'image%0.4d.pdf' % (frametime)
    fig.savefig(filename,format='pdf',dpi=400,rasterized=True)

#  ######################################################################
 
 
# ####################################################################
# ###################### Animating functions #########################
 
def init():
    pulse2.set_data([],[])
    time_text.set_text('')
    #surf1.set_array([])
    #surf2.set_array([])
    return pulse2,time_text


# ######################################################
 
def animate(it):
      
    print'time: %0.4f fs' %(time[it])
     
    if(draw_rz): 
    	global probrz
    	probrz[:,:] = 0.0
    	probrz = grp_rz[grp_rz.keys()[it]]
    	probrz *= jacobrz * drho
    	surf1.set_array(np.ravel(np.log10(probrz)))
    if(draw_rhoz):
        global probrhoz
        probrhoz[:,:] = 0.0
        probrhoz = grp_rhoz[grp_rhoz.keys()[it]]
        probrhoz *= jacobrhoz * dr
        surf2.set_array(np.ravel(np.log10(probrhoz)))
    
    pulse2.set_data(time[it],efield[it,1])
    time_text.set_text('Time: %3.2f fs' %(time[it]))
    
    if(draw_rz & draw_rhoz): 
        return pulse2,time_text,surf1,surf2
    elif(draw_rz):
    	return pulse2,time_text,surf1
    elif(draw_rhoz):
    	return pulse2,time_text,surf2
    	
# ######################################################

gs.tight_layout(fig)

framerate = 24 #int(10./graphicstime)

if(makemovie):  
    anim = animation.FuncAnimation(fig,animate,init_func=init,frames=ntime,interval=200,blit=True,repeat=False)
    anim.save('h2p_vib.mp4',fps=framerate,extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])

if((makemovie==0)&(makeframe==0)): 
	plt.show()


print('Script finished!!!')
