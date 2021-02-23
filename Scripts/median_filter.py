#-----------------------------------------------------------------------------#
#median_filter.py
#
#NPS Night Skies Program
#
#Last updated: 2021/02/22
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
import process_input as p 

#------------------------------------------------------------------------------#

def main():
	"""
	This function applies a median filter to the images. See the script 
	description for detail.
	"""
	#read in the pixel scale associated with the binning factor
	imgfile = glob(p.data_cal+'*sky*')[0]
	binning = fits.open(imgfile,uint=False)[0].header['XBINNING']
	F = pd.read_csv(p.calibration+'platescale.csv',index_col=0)
	pixscale = F['Scale'][binning] #[arcsec/pix]
	
	#median filter window size -- 1 degree
	w = n.ceil(3600/pixscale).astype(n.int) #[pix]
	
	#Median filter all the fisheye images
	l = len(p.data_cal)
	print('Applying 1-degree median filter to image:')
	for f in glob(p.data_cal+'[!MF]*sky*.fit'):
		print(f[l:])
		image = fits.open(f,uint=False)
		hdr = image[0].header
		filtered_img = median_filter(image[0].data, size=w)
		hdr['history'] = f'Median filtered with 1 degree ({w} pixel) window'
		fits.writeto(f[:l]+'MF_'+f[l:],filtered_img,header=hdr,overwrite=1)

		
if __name__ == '__main__':
	main()