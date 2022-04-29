#------------------------------------------------------------------------------#
#mask.py
#
#NPS Night Skies Program
#
#Last updated: 2022/3/3
#
#This script finds a circle to describe where the fisheye view is located. 
#
#Input: 
#   (1) mask_input.py, containing a flat image or an image with fisheye viewing 
#		area well lit and the brightness cutoff.
#
#Output:
#   (1) Image of a circular mask
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#
import importlib
import numpy as n

from astropy.io import fits
from matplotlib import pyplot as plt

# Local Source
import mask_input as mi

importlib.reload(mi)
#-----------------------------------------------------------------------------#

#Read in the file and the light pixels
flat = fits.open(mi.filein,uint=False)[0].data
bright = n.where(flat>mi.t) #location of bright pixels

#Find the xy center and the radius with no missing pixel 
center_y , center_x = n.mean(bright,axis=1)
radius_t = n.min(n.percentile(bright,99.99,axis=1) - n.mean(bright,axis=1))
radius = min(radius_t, center_y)


#Creat the mask image
x, y = n.meshgrid(n.arange(flat.shape[1]),n.arange(flat.shape[0]))
r = n.sqrt((x-center_x)**2 + (y-center_y)**2)
mask = n.zeros_like(r)
mask[n.where(r<=radius)] = 1
mask[n.where(mask==0)] = n.nan

#save the bestfit model mask
hdu = fits.PrimaryHDU()
hdu.header['THRESHOL'] = mi.t
hdu.header['CENTERX'] = center_x
hdu.header['CENTERY'] = center_y
hdu.header['RADIUS'] = radius
hdu.data = mask
hdu.writeto(mi.fileout, overwrite=True)

#plot
fig = plt.figure(1)
maskplot = n.zeros_like(mask)
maskplot[n.where(mask==1)]=n.nan
plt.imshow(flat-maskplot)
plt.show(block=False)



