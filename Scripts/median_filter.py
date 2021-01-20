#-----------------------------------------------------------------------------#
#median_filter.py
#
#NPS Night Skies Program
#
#Last updated: 2020/10/16
#
#This script applies a median filter with 1 degree window to the fisheye images.
#
#Input: 
#   (1) platescale.csv 	#platescale associated with the binning factor
#	(2) fisheye images to be filtered
#
#Output:
#   (1) Median filtered (MF_*) fits images
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#
import numpy as n
import pandas as pd

from astropy.io import fits
from glob import glob
from scipy.ndimage import median_filter

# Local Source
import filepath     
#------------------------------------------------------------------------------#

def main():
	"""
	This function applies a median filter to the images. See the script 
	description for detail.
	"""
	#read in the pixel scale associated with the binning factor
	imgfile = glob(filepath.data_cal+'*sky*')[0]
	binning = fits.open(imgfile,uint=False)[0].header['XBINNING']
	P = pd.read_csv(filepath.calibration+'platescale.csv',index_col=0)
	pixscale = P['Scale'][binning] #[arcsec/pix]
	
	#median filter window size
	w = n.ceil(3600/pixscale).astype(n.int) #[pix]
	
	#Median filter all the fisheye images
	l = len(filepath.data_cal)
	for f in glob(filepath.data_cal+'[!MF]*sky*.fit'):
		image = fits.open(f,uint=False)
		hdr = image[0].header
		filtered_img = median_filter(image[0].data, size=w)
		hdr['history'] = f'Median filtered with 1 degree ({w} pixel) window'
		fits.writeto(f[:l]+'MF_'+f[l:],filtered_img,header=hdr,overwrite=1)