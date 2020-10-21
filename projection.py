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

def fisheye_plot(theta, r, img):
	"""
	Plots the input image in a fisheye view. Saves the output as a png image. 
	"""
	plt.style.use('dark_background')
	fig = plt.figure()#figsize=(6,5)) 
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

#Mask - read in the fisheye mask to find the center of the view
mask = fits.open(filepath.mask,uint=False)[0].data
view = n.where(mask==1)
yc, xc = n.round(n.mean(view,axis=1))	#center of fisheye view
radius = n.min(n.floor((n.max(view,axis=1)-n.min(view,axis=1))/2)) #view radius

y = n.arange(mask.shape[0])- yc
x = n.arange(mask.shape[1])- xc
X, Y = n.meshgrid(x,y)

r = n.sqrt(X**2+Y**2) / radius
theta = n.arctan2(Y,X)

y1, x1 = n.max(view,axis=1)
y0, x0 = n.min(view,axis=1)

for f in glob(filepath.data_cal+'*light*'):
	img = fits.open(f,uint=False)[0].data 
	w = n.where(r<0)
	img[w] = 0#n.nan



r = r[y0:y1,x0:x1]
theta = theta[y0:y1,x0:x1]
img = img[y0:y1,x0:x1]

#-----------------------------------------------------------------------------#
#							Plot the fitting results						  #
#-----------------------------------------------------------------------------#
plt.rcParams['image.cmap'] = 'NPS_mag'
current_cmap = plt.cm.get_cmap()
current_cmap.set_bad(color='black')

plt.close('all')

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