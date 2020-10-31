#-----------------------------------------------------------------------------#
#photometric_calibration.py
#
#NPS Night Skies Program
#
#Last updated: 2020/10/15
#
#This script performs absolute photometric calibration using the image pixel 
#scale and the instrumental zeropoint. This script convers the image from 
#ADU/sec to magnitude per square arcsec. 
#  
#Input: 
#   (1) zeropoint.csv 	#zeropoint associated with the camera system
#	(2) platescale.csv 	#platescale associated with the binning factor
#	(3) fisheye images to be photometrically calibrated
#
#Output:
#   (1) photometrically calibrated images 
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import numpy as n
import pandas as pd

from astropy.io import fits
from glob import glob

# Local Source
import filepath

#-----------------------------------------------------------------------------#

def main():
	"""
	This script performs absolute photometric calibration. See the script 
	description for detail.
	"""
	
	#read in the zeropoint
	if filepath.use_default_zeropoint:	#use the default value
		Z = pd.read_csv(filepath.calibration+'zeropoint.csv',index_col=0)
		zp = Z['Zeropoint'][filepath.camera]
	else:								#use the measured value from the dataset
		Z = pd.read_csv(filepath.data_cal+'zeropoint.csv')
		zp = Z['Zeropoint'].mean()
		
	#read in the pixel scale associated with the binning factor
	imgfile = glob(filepath.data_cal+'*light*')[0]
	binning = fits.open(imgfile,uint=False)[0].header['XBINNING']
	P = pd.read_csv(filepath.calibration+'platescale.csv',index_col=0)
	pixscale = P['Scale'][binning]
	
	#platescale adjustment
	psa = 2.5*n.log10(pixscale**2) 
	
	#brightness calibration
	for f in glob(filepath.data_cal+'*light*.fit'):
		image = fits.open(f,uint=False,mode='update')
		hdr = image[0].header
		image[0].data = zp + psa - 2.5*n.log10(image[0].data/hdr['EXPTIME'])
		hdr['history'] = f'Zeropoint used for calibration is {zp}'
		hdr['history'] = f'Pixel scale used for calibration is {pixscale} "/pix'
		image.flush()
		image.close()
		
						 
if __name__ == '__main__':
	main()