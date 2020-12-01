#-----------------------------------------------------------------------------#
#centering.py
#
#NPS Night Skies Program
#
#Last updated: 2020/11/04
#
#This script corrects for fisheye lens distortion, puts the zenith to the center
#of images, and rotate the images to make the north pointing up. Correction or 
#the fisheye lens distortion and aligning the zenith to the image center have 
#not been implemented. The output overwrites the inmages in the calibrated data 
#folder. 
#
#Input: 
#   (1) filepath.data_cal+'center.txt'
#	(2) filepath.data_cal+filepath.reference
#	(3) filepath.mask
#   (4) filepath.data_cal+'*light*.fit'
#
#Output:
#   (1) filepath.data_cal+'*light*.fit'
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import ast
import numpy as n

from astropy.io import fits
from astropy.time import Time
from glob import glob
from skimage.transform import rotate

# Local Source
import filepath     
from sphericalgeometry import DistanceAndBearing

#-----------------------------------------------------------------------------#

def main():
	"""
	This script performs distortion correction and centering. See the script 
	description for detail.
	"""
	
	#Read in the coordinates of the image center
	file = open(filepath.data_cal+'center.txt', "r")
	center = ast.literal_eval(file.read())
	file.close()
	center_ra = center['ra']
	center_de = center['dec']
	ori = round(center['orientation'],1)
	
	#Compute zenith RA and Dec based on the observing location and time
	hdu = fits.open(filepath.data_cal+filepath.reference, fix=False)
	time = Time(hdu[0].header['DATE-OBS'])  #UTC observing date and time
	hdu.close()
	zenith_ra = time.sidereal_time('mean',longitude=filepath.long).degree
	zenith_de = filepath.lat
	
	#Compute the offset of the image centroid
	dab = DistanceAndBearing(center_de,center_ra,zenith_de,zenith_ra)	
	offset_distance = round(90-dab[0],2)
	offset_bearing = dab[1]
	print(f'Zenith is offset from the center by {offset_distance} degrees')
	
	#Mask - read in the fisheye mask center coordinates and radius
	mask = fits.open(filepath.mask,uint=False)[0].header
	xc, yc = mask['CENTERX'], mask['CENTERY']
	
	#Position calibration
	for f in glob(filepath.data_cal+'*sky*.fit'):
		image = fits.open(f,uint=False,mode='update')
		
		#correct for fisheye lens distorsion - need to be implemented
		pass	
		
		#align image centroid with the zenith - need to be implemented
		pass
		
		#correct for orientation; put north up
		image[0].data = rotate(image[0].data,ori,center=(xc,yc),mode='edge')
		image[0].header['history']=f'Image is rotated by {ori} degrees'	
		image[0].header['history']=f'Image is processed by {filepath.processor}'
		
		image.flush()
		image.close()
		
if __name__ == '__main__':
	main()

