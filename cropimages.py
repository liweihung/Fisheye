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
#Input: Test_Images/Coyote_Ridge/30sec.4_2x2_dark_subtracts.fit
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
filename = "Test_Images/Coyote_Ridge/30sec.4_2x2_dark_subtracts.fit"
hdulist = fits.open(filename)
hdr = hdulist[0].header
img = hdulist[0].data
(y,x) = img.shape

outfolder = "Test_Images/Coyote_Ridge_Crop/30sec.4_2x2_dark_subtracts"

for s in n.arange(200,1500,100):
	crop = img[int(y/2-s/2):int(y/2+s/2),int(x/2-s/2):int(x/2+s/2)]
	hdr['NAXIS1'] = s
	hdr['NAXIS2'] = s

	outfile = outfolder+'_crop%i.fit'%s
	fits.writeto(outfile, crop, hdr, overwrite=True)

