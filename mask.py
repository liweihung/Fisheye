#-----------------------------------------------------------------------------#
#mask.py
#
#NPS Night Skies Program
#
#Last updated: 2020/08/25
#
#This script fit a circle to the image to find out where the fisheye view is 
#located. 
#
#Input: 
#   (1) A flat image or an image with fisheye viewing area well lit
#
#Output:
#   (1) The best fit model of a circular mask
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import numpy as n

from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.optimize import leastsq

# Local Source
import filepath     

#-----------------------------------------------------------------------------#
#--------------       Change parameters here        --------------------------#
#-----------------------------------------------------------------------------#
#Inputfile 
#flatfile = 'flat_2x2_20200623.fit'
flatfile = 'flat_4x4_rebinned_20200623.fit'

t = 0.1 #pixels greater than this threshold will be assummed it received light

#fitting parameters:
p0 = [1186,773,747]  #initial parameter: x_center, y_center, radius 
s = 1e-1 #stepsize

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

if __name__ == "__main__":
	#Readin the file and assign light and dark pixels
	flat = fits.open(filepath.calibration+flatfile,uint=False)[0].data
	flat[n.where(flat>t)] = 1
	flat[n.where(flat<=t)] = 0
	
	#fit to a circle
	popt = leastsq(errorfunction, p0, args=(flat), epsfcn=s)[0]
	bestfitmodel = circlemodel(*popt, flat.shape)
	print('Best fit x, y, r:', n.round(popt,2))
	
	#save the bestfit model mask
	maskfile = filepath.calibration+flatfile[:-4]+'_mask.fit'
	mask = bestfitmodel.copy()
	mask[mask==0] = n.nan
	fits.writeto(maskfile, mask, overwrite=True)
	
	#plot
	plt.imshow(flat-bestfitmodel)
	plt.show(block=False)
	
	

