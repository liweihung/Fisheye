#-----------------------------------------------------------------------------#
#cropimages.py
#
#NPS Night Skies Program
#
#Last updated: 2020/04/28
#
#This script crops the center of the images and copy the header with 
#modifications on the image size entries.
#
#Input: Test_Images/20200427/40_sec_V_light_4x4.fit
#
#Output:
#   (1) A set of cropped images in the Test_Images/Coyote_Ridge_Crop folder
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import numpy as n

from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
#-----------------------------------------------------------------------------#

# Load the FITS hdulist using astropy.io.fits
filename = "Test_Images/20200427/40_sec_V_light_4x4.fit"
hdulist = fits.open(filename)

outfolder = "Test_Images/20200427_Crop/40_sec_V_light_4x4"

def center_crop(hdulist, interval):
	"""
	Crop the center of the image
	hdulist: input image file
	interval: list of croped image sizes 
	"""
	hdr = hdulist[0].header
	img = hdulist[0].data
	(y,x) = img.shape
	for s in interval:
		crop = img[int(y/2-s/2):int(y/2+s/2),int(x/2-s/2):int(x/2+s/2)]
		hdr['NAXIS1'] = s
		hdr['NAXIS2'] = s
		outfile = outfolder+'_crop%i.fit'%s
		fits.writeto(outfile, crop, hdr, overwrite=True)


def general_crop(hdulist,side,nsquare):
	"""
	Crop the image into smaller squares
	"""
	hdr = hdulist[0].header
	img = hdulist[0].data
	(y,x) = img.shape
	hdr['NAXIS1'] = side
	hdr['NAXIS2'] = side
	yi, xi = map(int,0.5*(n.array([y,x])-side*nsquare))
	yf, xf = map(int,0.5*(n.array([y,x])+side*nsquare))
	for i in n.arange(xi,xf,side):
		for j in n.arange(yi,yf,side):
			crop = img[j:j+side,i:i+side]
			outfile = outfolder+'_crop%i_%i.fit'%(i,j)
			fits.writeto(outfile, crop, hdr, overwrite=True)

	
#center_crop(hdulist,n.arange(500,1100,100))
general_crop(hdulist,280,5)	

	
