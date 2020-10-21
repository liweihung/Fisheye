#------------------------------------------------------------------------------#
#mask.py
#
#NPS Night Skies Program
#
#Last updated: 2020/10/07
#
#This script fit a circle to the image to find out where the fisheye view is 
#located. 
#
#Input: 
#   (1) mask_input.py, containing a flat image or an image with fisheye viewing 
#		area well lit and fitting parameters
#
#Output:
#   (1) The best fit model image of a circular mask
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#
import numpy as n

from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.optimize import leastsq
   
# Local Source
import mask_input as mi
#-----------------------------------------------------------------------------#
	
def circlemodel(xc, yc, rc, imgshape):
	'''
	This module returns the modeled image of a solid-filled circle given the 
	center coordinates, the radius of the circle, and the overall image shape.
	
	Parameters
	----------
	xc : a number
		X center of the circle.
	yc : a number
		Y center of the circle.
	rc : a number
		radius of the circle.
	imgshape : list or tuple of numbers
		shape of the overall image.		
		
	Returns
	-------
	model : 2D array
		Model image of a solid circle.   
	'''
	x, y = n.meshgrid(n.arange(imgshape[1]),n.arange(imgshape[0]))
	r = n.sqrt((x-xc)**2 + (y-yc)**2)
	model = n.zeros_like(r)
	model[n.where(r<=rc)] = 1
	return model
	

def errorfunction(xyr, img):
	'''
	This module computes the pixelwise difference a modeled circled imgage and 
	the input image. 
	
	Parameters
	----------
	xyr : list of numbers
		X center, y center, and radius in pixels of the circle.
	img : 2D array
		Input image for fitting; usually a flatfield image.  
		
	Returns
	-------
	delta : list of numbers
		Difference between model and image for each pixel.
	'''
	model = circlemodel(*xyr, img.shape)
	delta = abs(model-img).ravel()
	print(n.sum(delta, dtype=int), xyr)
	return delta
	
#-----------------------------------------------------------------------------#

#Readin the file and assign light and dark pixels
flat = fits.open(mi.filein,uint=False)[0].data
flat[n.where(flat>mi.t)] = 1
flat[n.where(flat<=mi.t)] = 0

#fit to a circle
popt = leastsq(errorfunction, mi.p0, args=(flat), epsfcn=mi.s)[0]
bestfitmodel = circlemodel(*popt, flat.shape)
print('Best fit x, y, r:', n.round(popt,2))

#save the bestfit model mask
mask = bestfitmodel.copy()
mask[mask==0] = n.nan
hdu = fits.PrimaryHDU()
hdu.header['CENTERX'] = n.int(popt[0])
hdu.header['CENTERY'] = n.int(popt[1])
hdu.header['RADIUS'] = n.int(popt[2])
hdu.data = mask
hdu.writeto(mi.fileout, overwrite=True)
#fits.writeto(mi.fileout, mask, header=hdr, overwrite=True)

#plot
plt.imshow(flat-bestfitmodel)
plt.show(block=False)



