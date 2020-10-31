#-----------------------------------------------------------------------------#
#centering.py
#
#NPS Night Skies Program
#
#Last updated: 2020/09/22
#
#This script centers the zenith to the center of the fisheye image
#
#Input: 
#   (1) 
#
#Output:
#   (1) 
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
from matplotlib import pyplot as plt
from scipy.ndimage import shift
from skimage.transform import rotate
#from scipy.optimize import leastsq

# Local Source
import filepath     
from sphericalgeometry import DistanceAndBearing

#-----------------------------------------------------------------------------#

#Read in the coordinates of the image center
file = open(filepath.data_cal+'center.txt', "r")
center = ast.literal_eval(file.read())
file.close()
center_ra = center['ra']
center_de = center['dec']
orientation = round(center['orientation'],1)

#Compute zenith RA and Dec based on the observing location and time
hdu = fits.open(filepath.data_cal+filepath.reference, fix=False)
time = Time(hdu[0].header['DATE-OBS'])  #UTC observing date and time
hdu.close()
zenith_ra = time.sidereal_time('mean',longitude=filepath.long).degree #[degree]
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
for f in glob(filepath.data_cal+'*light*.fit'):
	image = fits.open(f,uint=False,mode='update')
	
	#correct for orientation; put north up
	image[0].data = rotate(image[0].data,orientation,center=(xc,yc))#,cval=n.nan)
	image[0].header['history'] = f'Image is rotated by {orientation} degrees'
	
	#correct for fisheye lens distorsion - need to be implemented
	pass 
	
	#align image centroid with the zenith - need to be implemented
	pass
	
	image[0].header['history'] = f'Image is processed by {filepath.processor}'
	image.flush()
	image.close()

plt.show()

#shift and then rotate
