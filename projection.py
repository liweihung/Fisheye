#-----------------------------------------------------------------------------#
#projection.py
#
#NPS Night Skies Program
#
#Last updated: 2020/07/22
#
#This script 
#Input: 
#   (1) filepath.mask to get the x,y center and the fisheye view radius
#	(2) all the processed fisheye images
#
#Output:
#   (1) 
#   (2) 
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#
import matplotlib.projections as mprojections
import numpy as n

from astropy.io import fits
from glob import glob
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import Wedge
from matplotlib.projections.geo import HammerAxes
import matplotlib.spines as mspines

# Local Source
import colormaps
import filepath
import upper_hammer

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#						  Define Plot settings							   	   #
#------------------------------------------------------------------------------#
def fisheye_plot(img):
	"""
	Plots the input image in a fisheye view. Saves the output as a png image. 
	"""
	fig = plt.figure() 
	ax = plt.subplot(111, projection='polar')
	ax.pcolormesh(theta,r_deg,img,vmin=14,vmax=24)
	ax.set_rlim(0,90)
	ax.set_yticklabels([])
	ax.tick_params(colors='darkgray')
	ax.set_theta_zero_location('N') 
	ax.grid(True, color='gray', linestyle='dotted', linewidth=.5)
	plt.savefig(f[:-4]+'_fisheye.png', dpi=300)
	
	
def hammer_plot(img):
	"""
	Plots the input image in Hammer equal area projection. Saves the output as 
	a png image.
	"""
	plt.figure(figsize=(15,5.2))
	ax = plt.subplot(111, projection="upper_hammer")
	ax.pcolormesh(theta_sorted,r_sorted,img,vmin=14,vmax=24)
	ax.grid(True)
	#plt.title("Hammer")
	plt.tight_layout(rect=(0.03,-0.6,0.98,0.97))
	plt.savefig(f[:-4]+'_hammer.png')


#General plot settings
plt.close('all')
plt.style.use('dark_background')
plt.rcParams['image.cmap'] = 'NPS_mag'
plt.cm.get_cmap().set_bad(color='black')	

#Mask - read in the fisheye mask center coordinates and radius
mask = fits.open(filepath.mask,uint=False)[0].header
xc, yc, r0 = mask['CENTERX'], mask['CENTERY'], mask['RADIUS']
X, Y = n.meshgrid(n.arange(-r0,r0),n.arange(-r0,r0))

#Polar coordinates
r = n.sqrt(X**2+Y**2) / r0
theta = n.arctan2(Y,X)

#Fisheye takes r in degree
r_deg = 90*r

#Hammer plot requires the values to be sorted
r_str = n.pi/2 - r * n.pi/2
inds = n.argsort(theta[:,0])
theta_sorted = theta[inds,:] 
r_sorted = r_str[inds,:]

#------------------------------------------------------------------------------#
#				Plot the image in fisheye and Hammer projections			   #
#------------------------------------------------------------------------------#
for f in glob(filepath.data_cal+'*light*.fit'):
	print(f)
	img = fits.open(f,uint=False)[0].data[yc-r0:yc+r0,xc-r0:xc+r0] 
	#fisheye_plot(img)
	hammer_plot(img[inds,:])

#plt.show(block=False)
#plt.savefig(f+'_hammer4.jpg', dpi=300)

'''
#----------------old code below------------------------------------------------#

plt.figure()
plt.rcParams['image.cmap'] = 'NPS_mag'
current_cmap = plt.cm.get_cmap()
current_cmap.set_bad(color='black')
plt.imshow(fits.open(f,uint=False)[0].data,vmin=14,vmax=24)

plt.colorbar()
plt.show(block=False)

'''