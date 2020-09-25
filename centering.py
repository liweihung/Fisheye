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
#import numpy as n

from astropy.io import fits
from astropy.time import Time
#from matplotlib import pyplot as plt
#from scipy.optimize import leastsq

# Local Source
import filepath     
#-----------------------------------------------------------------------------#

#Read in the coordinates of the image center
file = open(filepath.data_cal+'center.txt', "r")
center = ast.literal_eval(file.read())
file.close()
center_ra = center['ra']
center_de = center['dec']
print(center_ra, center_de)


#Compute zenith RA and Dec based on the observing location and time
hdu = fits.open(filepath.data_cal+'40_sec_V_light_4x4.fit', fix=False)[0] 
time = Time(hdu.header['DATE-OBS'])  #UTC observing date and time
latitude, longitude = 40.696034, -104.600301 #degrees
zenith_ra = time.sidereal_time('mean',longitude=longitude).degree #[degree]
zenith_de = latitude
print (zenith_ra, zenith_de)											  #[degree]


#shift and then rotate
