#-----------------------------------------------------------------------------#
#projection.py
#
#NPS Night Skies Program
#
#Last updated: 2020/11/05
#
#This script plots the fits images in fisheye and Hammer projections. 
#
#Input: 
#   (1) filepath.mask to get the x,y center and the fisheye view radius
#	(2) all the processed fisheye fit images
#
#Output:
#   (1) *fisheye.png
#   (2) *hammer.png
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#
import copy
import matplotlib as mpl
import numpy as n
import warnings

from astropy.io import fits
from glob import glob
from matplotlib import pyplot as plt
from skimage.transform import rotate

# Local Source
import colormaps
import filepath
import upper_hammer

#------------------------------------------------------------------------------#
def main():
	"""
	This script plots fits images in fisheye and Hammer projection. See the 
	script description for detail.
	"""
	#--------------------------------------------------------------------------#
	#						  Generate Polar Coordinates				   	   #
	#--------------------------------------------------------------------------#
	#Mask - read in the fisheye mask center coordinates and radius
	mask = fits.open(filepath.mask,uint=False)[0].header
	xc, yc, r0 = mask['CENTERX'], mask['CENTERY'], mask['RADIUS']
	X, Y = n.meshgrid(n.arange(-r0,r0),n.arange(-r0,r0))
	
	#Polar coordinates
	r = n.sqrt(X**2+Y**2) / r0
	theta = -n.arctan2(Y,X)
	
	#Fisheye takes r in degree
	r_deg = 90 * r
	theta_f = theta + n.pi/2
	
	#Hammer plot requires the values to be sorted
	r_str = n.pi/2 - r * n.pi/2
	inds = n.argsort(theta[:,0])
	theta_s = theta[inds,:] 
	r_s = r_str[inds,:]
	
	
	#--------------------------------------------------------------------------#
	#						  Define Plot settings						   	   #
	#--------------------------------------------------------------------------#
	#General plot settings
	plt.close('all')
	plt.style.use('dark_background')
	plt.rcParams['image.cmap'] = 'NPS_mag'
	cmap = copy.copy(mpl.cm.get_cmap("NPS_mag"))
	cmap.set_bad(color='black')
	
	#Fisheye plot setting
	fig0 = plt.figure('fisheye') 
	ax0 = fig0.add_subplot(111, projection='polar')
	ax0.set_rlim(0,90)
	ax0.set_yticklabels([])
	ax0.tick_params(colors='darkgray')
	ax0.set_theta_zero_location('N') 
	
	#Hammer plot setting
	fig1 = plt.figure('hammer',figsize=(15,5.2))
	ax1 = fig1.add_subplot(111, projection="upper_hammer")
	fig1.tight_layout(rect=(0.03,-0.6,0.98,0.97))
	
	#Suppressing a MatPlotLib benign warning about pcolormesh shading
	warnings.filterwarnings("ignore",category=UserWarning)

	#--------------------------------------------------------------------------#
	#				Plot the image in fisheye and Hammer projections		   #
	#--------------------------------------------------------------------------#
	for i, f in enumerate(glob(filepath.data_cal+'*sky*.fit')):
		print(f"projecting image {i}")
		img = fits.open(f,uint=False)[0].data[yc-r0:yc+r0,xc-r0:xc+r0] 
		img_hammer = rotate(img.astype('float32'),-90,cval=n.nan)[inds,:]
		
		#plot fisheye
		ax0.pcolormesh(theta_f,r_deg,img,shading='auto',vmin=14,vmax=24)
		ax0.grid(True, color='gray', linestyle='dotted', linewidth=.5)
		fig0.savefig(f[:-4]+'_fisheye.png', dpi=300)
		
		#plot hammer
		ax1.pcolormesh(theta_s,r_s,img_hammer,shading='auto',vmin=14,vmax=24)
		ax1.grid(True)
		fig1.savefig(f[:-4]+'_hammer.png')


if __name__ == '__main__':
	main()	