#-----------------------------------------------------------------------------#
#projection.py
#
#NPS Night Skies Program
#
#Last updated: 2021/02/22
#
#This script plots the fits images in fisheye and Hammer projections. 
#
#Input: 
#   (1) reading in the mask to get the x,y center and the fisheye view radius
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
from matplotlib.transforms import Transform
import numpy as n
import warnings

from astropy.io import fits
from astropy.time import Time, TimeDelta
from glob import glob
from matplotlib import pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from skimage.transform import rotate

# Local Source
import colormaps
import process_input as p 
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
	mask = fits.open(p.mask,uint=False)[0].header
	xc, yc, r0 = int(mask['CENTERX']), int(mask['CENTERY']), int(mask['RADIUS'])
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
	#fig0 = plt.figure('fisheye') 
	#ax0 = fig0.add_subplot(111, projection='polar')
	#ax0.set_rlim(0,90)
	#ax0.set_yticklabels([])
	#ax0.tick_params(colors='darkgray')
	#ax0.set_theta_zero_location('N') 
	
	#Hammer plot setting
	fig1 = plt.figure('hammer',figsize=(15,5.2))
	ax1 = fig1.add_subplot(111, projection="upper_hammer")
	fig1.tight_layout(rect=(0.03,-0.6,0.98,0.97))
	colorbar = plt.imread(p.calibration+'colorbar.jpg')
	imagebox = OffsetImage(colorbar, zoom=0.47)
	imagebox.image.axes = ax1
	ab = AnnotationBbox(imagebox,(0.19, 0.94),xycoords='figure fraction',
						frameon=False)
	ax1.add_artist(ab)
	for mag in range(14,25):
		ax1.text(-0.017+0.0348*(mag-14),1.06,mag,color='gray',fontsize=10,
				 ha='left',va='top', transform=ax1.transAxes)
	ax1.text(0.14,1.04,r'V mags arcsec$^{-2}$',color='gray',fontsize=10,
			 ha='center',va='top', transform=ax1.transAxes)

	#Suppressing a MatPlotLib benign warning about pcolormesh shading
	warnings.filterwarnings("ignore",category=UserWarning)

	#--------------------------------------------------------------------------#
	#				Plot the image in fisheye and Hammer projections		   #
	#--------------------------------------------------------------------------#
	for f in glob(p.data_cal+'img-00*-sky*.fit'):

		print('projecting ' + f[len(p.data_cal):])
		imgf = fits.open(f,uint=False)[0]
		hdr = imgf.header
		img = imgf.data[yc-r0:yc+r0,xc-r0:xc+r0]
		img_hammer = rotate(img.astype('float32'),-90,cval=n.nan)[inds,:]
		
		#plot fisheye
		#ax0.pcolormesh(theta_f,r_deg,img,shading='flat',vmin=14,vmax=24)
		#ax0.grid(True, color='gray', linestyle='dotted', linewidth=.5)
		#fig0.savefig(f[:-4]+'_fisheye.png', dpi=300)
		
		#plot hammer
		ax1.pcolormesh(theta_s,r_s,img_hammer,vmin=14,vmax=24)
		ax1.grid(True)
		ax1.text(0.5,1.07,hdr['LOCATION'], color='gray', fontsize=16, 
				 ha='center',va='center', transform=ax1.transAxes)
		t = Time(hdr['DATE-OBS'])+TimeDelta(p.UTCoffset*3600,format='sec')
		obstime = str(t.datetime.date())+'  '+ \
				  str(round(t.datetime.hour+t.datetime.minute/60,1))+' hours'
		ta = ax1.text(1.0,1.07,obstime, color='gray', fontsize=16, 
				 ha='right',va='center', transform=ax1.transAxes)
		fig1.savefig(f[:-4]+'_hammer.png')
		ta.remove()#ax1.clear()
		#plt.show(block=False)

if __name__ == '__main__':
	main()	