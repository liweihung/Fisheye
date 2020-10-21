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
#-----------------------------------------------------------------------------#
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

#------------------------------------------------------------------------------#
#					  Create Upper Hammer Projection Class				   	   #
#------------------------------------------------------------------------------#
class UpperHammerAxes(HammerAxes):
	name = 'upper_hammer'
	def cla(self):
		HammerAxes.cla(self)
		Axes.set_xlim(self, -n.pi, n.pi)
		Axes.set_ylim(self, 0, n.pi / 2.0)
		self.xaxis.set_ticks_position('top')
		self.tick_params(axis='x', length=0)
		self.tick_params(colors='darkgray')
		self.grid(color='gray', linestyle='dotted', linewidth=.5)
	
	def _gen_axes_patch(self):
		return Wedge((0.5, 0.5), 0.5, 0, 180)
	
	def _gen_axes_spines(self):
		pass
		path = Wedge((0, 0), 1.0, 0, 180).get_path()
		spine = mspines.Spine(self, 'circle', path)
		spine.set_color('gray')
		spine.set_patch_circle((0.5, 0.5), 0.5)
		return {'wedge':spine}
		
mprojections.register_projection(UpperHammerAxes)

#------------------------------------------------------------------------------#
#						  Define Plot settings							   	   #
#------------------------------------------------------------------------------#
plt.close('all')
plt.style.use('dark_background')
plt.rcParams['image.cmap'] = 'NPS_mag'
plt.cm.get_cmap().set_bad(color='black')

def fisheye_plot(theta, r, img):
	"""
	Plots the input image in a fisheye view. Saves the output as a png image. 
	"""
	fig = plt.figure() 
	ax = plt.subplot(111, projection='polar')
	ax.pcolormesh(theta,r,img,vmin=14,vmax=24)
	ax.set_rlim(0,90)
	ax.set_yticklabels([])
	ax.tick_params(colors='darkgray')
	ax.set_theta_zero_location('N') 
	ax.grid(True, color='gray', linestyle='dotted', linewidth=.5)
	#plt.savefig('.png', dpi=200)
	
	
def hammer_plot(theta,r,img):
	"""
	Plots the input image in Hammer equal area projection. Saves the output as 
	a png image.
	"""
	plt.figure(figsize=(11,3.8))
	ax = plt.subplot(111, projection="upper_hammer")
	ax.pcolormesh(theta,r,img,vmin=14,vmax=24)
	ax.grid(True)
	#plt.title("Hammer")
	plt.tight_layout(rect=(0.03,-0.6,0.98,0.97))
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#			  Input files												   	   #
#------------------------------------------------------------------------------#

#Mask - read in the fisheye mask center coordinates and radius
mask = fits.open(filepath.mask,uint=False)[0].header
xc, yc, r0 = mask['CENTERX'], mask['CENTERY'], mask['RADIUS']
X, Y = n.meshgrid(n.arange(-r0,r0),n.arange(-r0,r0))

#Polar coordinates
r = n.sqrt(X**2+Y**2) / r0
theta = n.arctan2(Y,X)

for f in glob(filepath.data_cal+'*light*'):
	img = fits.open(f,uint=False)[0].data[yc-r0:yc+r0,xc-r0:xc+r0] 

#-----------------------------------------------------------------------------#
#							Plot the fitting results						  #
#-----------------------------------------------------------------------------#
#Fisheye takes r in degree
r_deg = 90*r
fisheye_plot(theta, r_deg, img)

#Hammer plot requires the values to be sorted
r_str = n.pi/2 - r * n.pi/2
inds = n.argsort(theta[:,0])
theta_sorted = theta[inds,:] 
r_sorted = r_str[inds,:]
img_sorted = img[inds,:]
hammer_plot(theta_sorted, r_sorted, img_sorted)

plt.show(block=False)
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