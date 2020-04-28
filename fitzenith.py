#-----------------------------------------------------------------------------#
#fitzenith.py
#
#NPS Night Skies Program
#
#Last updated: 2019/09/17
#
#This script...
#	(1) 
#	(2)  
#	(3) 
#   (4) 
#
#Note: 
#
#
#Input: 
#   (1) 
#   (2) 
#   (3) 
#
#Output:
#   (1) 
#   (2) 
#
#History:
#	Li-Wei Hung -- Created
#
#-----------------------------------------------------------------------------#
from astropy.io import fits
#from datetime import datetime as Dtime
from glob import glob, iglob
from scipy import ndimage
from scipy import optimize
#from scipy.misc import imsave 

#import pdb
import matplotlib.pyplot as plt
import numpy as n

# Local Source
from standard import StandardStarMap
from sphericalgeometry import DistanceAndBearing
#-----------------------------------------------------------------------------#
#for f in iglob(imgpath+'*.FIT'):
#	img = fits.open(f)[0]
#	meanADU = img.data.mean()#[840:1240,1348:1748].mean()
#	print(f, meanADU)

#-----------------------------------------------------------------------------#
#Load the test image
imgpath = 'C:/Users/lhung/Pictures/20190823_ROMO _astronomy/'
img = fits.open(imgpath+'ASICAP_2019-08-23_21_39_22_281.FIT')[0]
rotated_img = n.flip(ndimage.rotate(img.data, -90),0) #shape is (3096,2080)
stars_only = rotated_img.copy()[1048:2048,540:1540] #shape is (1000,1000)
stars_only[n.where(stars_only<65500)] = 0
stars_only[n.where(stars_only>65500)] = 1

fits.writeto(imgpath+'ASICAP_2019-08-23_21_39_22_281_center400.FIT', rotated_img.copy()[1348:1748,840:1240], overwrite=True)

#plt.figure()
#plt.imshow(stars_only, origin='bottom')
#plt.colorbar()

#variables: 1. plate scale 2. rotation angle 3. x offset 4. y offset 

#Glacier Basin Campground
lat = 40.325979
lon = -105.598273
elevation = 2591
time = img.header['DATE-OBS'] #UTC

#Standard stars above the horizon
#Hv.columns is ['HIP', 'Vmag', 'RA', 'DE', 'B-V', 'V-I', 'AZ', 'ALT'],
Hv = StandardStarMap(lat,lon,elevation,time)
ps = 0.052 #pixel scale [deg/pixel]
Hv['X'] = -(90-Hv.ALT) * n.sin(n.deg2rad(Hv.AZ)) / ps 
Hv['Y'] = (90-Hv.ALT) * n.cos(n.deg2rad(Hv.AZ)) / ps

def PolarPlot(az, alt):
	"""
	Plot the given azimuth and altitude points with the correct polar setting
	az = azimuth in degrees
	alt = altitude in degrees
	"""
	fig = plt.figure()
	ax = plt.subplot(111, projection='polar')
	c = ax.scatter(n.deg2rad(az), 90-alt, c=Hv.Vmag, cmap='jet', alpha=0.75)
	ax.set_rlim(0,90)
	ax.set_yticks(n.linspace(0, 90, 10))
	ax.set_yticklabels(ax.get_yticks()[::-1])
	ax.set_theta_zero_location('N') 
	ax.grid(True)
	cbar = fig.colorbar(c, ax=ax)
	cbar.set_clim(-1.5,3.5)

#plot the standard stars in polar coordinate
#PolarPlot(Hv.AZ, Hv.ALT)

def ErrorFunction(p):#l, a, ps
	"""
	l: new latitude origin
	a: new azimuth origin
	ps: plate scale
	"""
	print(p)
	ALT_new, AZ_new = DistanceAndBearing(p[0],p[1],Hv.ALT,Hv.AZ)
	Hv['Xnew'] = (-(90-ALT_new)*n.sin(n.deg2rad(AZ_new))/p[2]).astype('int64')
	Hv['Ynew'] = ((90-ALT_new)*n.cos(n.deg2rad(AZ_new))/p[2]).astype('int64')
	Hw = Hv[(abs(Hv.Xnew)<500) & (abs(Hv.Ynew)<500)]
	err = -1*stars_only[(Hw.Ynew.values+500,Hw.Xnew.values+500)]
	print(len(err))
	return err

#plot the standard stars with the new origin in polar coordinate
#PolarPlot(AZ_new, ALT_new) 

p = optimize.leastsq(ErrorFunction, (85, 260., ps), epsfcn=1e-4, full_output=1) 
print(p)

#plt.show(block=False)

#run client.py -k cllxijkpvxsibace --upload center400_b.FIT --newfits center400_wcs.FIT
