#-----------------------------------------------------------------------------#
#positional_calibration.py
#
#NPS Night Skies Program
#
#Last updated: 2021/02/19
#
#This script corrects for fisheye lens distortion, puts the zenith to the center
#of images, and rotate the images to make the north pointing up. Correction or 
#the fisheye lens distortion and aligning the zenith to the image center have 
#not been implemented. The output overwrites the inmages in the calibrated data 
#folder. 
#
#Input: 
#   (1) center.txt
#	(2) reference image for observing time
#	(3) mask
#   (4) sky images
#
#Output:
#   (1) sky images corrected for rotation, (centering, and distortion) 
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import ast
import numpy as n

from astropy.coordinates import EarthLocation
from astropy.io import fits
from astropy.time import Time
from glob import glob
from skimage.transform import rotate

# Local Source
import process_input as p 
from sphericalgeometry import DistanceAndBearing

#-----------------------------------------------------------------------------#

def main():
	"""
	This script performs distortion correction and centering. See the script 
	description for detail.
	"""
	
	#Read in the coordinates of the image center
	try:
		file = open(p.data_cal+'center.txt', "r")
	except(FileNotFoundError):
		print('center.txt is not found. Position calibration is not performed.')
		return
	center = ast.literal_eval(file.read())
	file.close()
	center_ra = center['ra']
	center_de = center['dec']
	ori = round(center['orientation'],1)
	
	#Compute zenith RA and Dec based on the observing location and time
	hdu = fits.open(p.data_cal+p.reference, fix=False)
	hdr = hdu[0].header
	time = Time(hdr['DATE-OBS'])  #UTC observing date and time
	c = EarthLocation(lon=hdr['SITELONG'] , lat=hdr['SITELAT'])
	zenith_ra = time.sidereal_time('mean',longitude=c.lon.deg).degree
	zenith_de = c.lat.deg
	hdu.close()
	
	#Compute the offset of the image centroid
	dab = DistanceAndBearing(center_de,center_ra,zenith_de,zenith_ra)	
	offset_distance = round(90-dab[0],2)
	offset_bearing = dab[1]
	print(f'Zenith is offset from the center by {offset_distance} degrees')
	
	#Mask - read in the fisheye mask center coordinates and radius
	mask = fits.open(p.mask,uint=False)[0].header
	xc, yc = mask['CENTERX'], mask['CENTERY']
	
	#Position calibration
	for f in glob(p.data_cal+'*sky*.fit'):
		image = fits.open(f,uint=False,mode='update')
		
		#correct for fisheye lens distorsion - need to be implemented
		pass	
		
		#align image centroid with the zenith - need to be implemented
		pass
		
		#correct for orientation; put north up
		imgdata = image[0].data.astype('float32')
		image[0].data = rotate(imgdata,ori,center=(xc,yc),mode='edge')
		image[0].header['history']=f'Image is rotated by {ori} degrees'	
		
		image.flush()
		image.close()
		
if __name__ == '__main__':
	main()

