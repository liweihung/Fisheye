#-----------------------------------------------------------------------------#
#projection.py
#
#NPS Night Skies Program
#
#Last updated: 2020/07/22
#
#This script 
#Input: 
#   (1) 
#
#Output:
#   (1) 
#   (2) 
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import numpy as n

from astropy.io import fits
from matplotlib import pyplot as plt

import colormaps


#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#			  Input files													  #
#-----------------------------------------------------------------------------#
#original image 
forig = 'Test_Images/20200427/40_sec_V_light_4x4.fit'
hdu_orig = fits.open(forig, fix=False)[0] #open the original image file
img = hdu_orig.data
img_center = [773,1180]
y = n.arange(img.shape[0])- img_center[0]
x = n.arange(img.shape[1])- img_center[1]
X, Y = n.meshgrid(x,y)
image_radius = 743# [pix]
r = n.pi/2 - n.sqrt(X**2+Y**2) / image_radius *n.pi/2
theta = n.arctan2(Y,X)


y0 = int(img_center[0]-image_radius)
y1 = int(img_center[0]+image_radius)
x0 = int(img_center[1]-image_radius)
x1 = int(img_center[1]+image_radius)

w = n.where(r<0)
img[w] = 220#n.nan

r = r[y0:y1,x0:x1]
theta = theta[y0:y1,x0:x1]
img = img[y0:y1,x0:x1]

img2 = img-209 #arbitrary reduction; bias? dark?
inds = n.argsort(theta[:,0])
theta=theta[inds,:]
r=r[inds,:]
img=img[inds,:]-209 #arbitrary reduction; bias? dark?

#brightness calibration
psa = 2.5*n.log10(380**2) # platescale adjustment
mc = 8.88+psa-2.5*n.log10(img/40)
mc2 = 8.88+psa-2.5*n.log10(img2/40)
#-----------------------------------------------------------------------------#
#							Plot the fitting results						  #
#-----------------------------------------------------------------------------#
plt.figure()
plt.rcParams['image.cmap'] = 'NPS_mag'
current_cmap = plt.cm.get_cmap()
current_cmap.set_bad(color='black')

plt.subplot(111, projection="hammer")
plt.pcolormesh(theta,r,mc,vmin=14,vmax=24)#,norm=norm,cmap=cmap)
plt.title("Hammer")
plt.grid(True)

plt.figure()
plt.imshow(mc2,vmin=14,vmax=24)

plt.colorbar()


#general plot setting
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('signal-to-noise ratio')
#plt.legend(loc='upper right', frameon=False, numpoints=1)
#plt.xlabel('Airmass')
#plt.ylabel('M-m')
#plt.title('Zeropoint and Extinction Coefficient')
#imgout = 'Plots/zeropoint_and_extinction_fit.png'
#plt.savefig(imgout, dpi=300)
plt.show(block=False)


#----------------old code below------------------------------------------------#
#from astropy.coordinates import EarthLocation
#elev = 1580. #meters
#location = EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=elev*u.m)
#c0 = [time.sidereal_time('mean',longitude=long).degree, lat]#[RA,Dec] in degrees






