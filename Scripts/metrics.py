#------------------------------------------------------------------------------#
#metrics.py
#
#NPS Night Skies Program
#
#Given the calibrated fisheye images, the script extracts various night sky 
#brightness metrics 
#
#Input: 
#   (1) 
#	(2) 
#
#Output:
#   (1) 
#	(2) 
#	(3) 
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#
import numpy as n
import os
import pandas as pd

from astropy.io import fits
from glob import glob

# Local Source
import process_input as p     

#------------------------------------------------------------------------------#

def zenith(img, r=4.75):
    """
	Calculate the median zenith brightness. 

	Parameters
	----------
	img : 2D array
        Fisheye image
        
	r : number, optional
        radius in pixel; 1 pix = 379 arcsec; aperture diameter should be set to
        about 1 degree, comparable to CCD zenith aperture 
		
	Returns
	-------
	zenith_mag : float
		median zenith brightness in mag per sequared arcsec
  
	"""
    zenith_mag = n.median(img[n.where(R<r)])
    return zenith_mag

def illuminance_horizontal(img):
    """
	Calculate the horizontal illuminance. 

	Parameters
	----------
	img : 2D array
        Fisheye image
        
	Returns
	-------
	Eh : float
		horizontal illuminance in mlx
  
	"""
    bi = 108.48*n.exp(20.7233-0.92104*img) # [ucd m-2], Dan's conversion formula
    Ei = bi*Omega/1000 # [mlx]
    Eh = n.nansum(Ei * n.cos(Theta)) #Horizontal illuminance [mlx]
    return Eh


def illuminance_vertical(img):
    """
	Calculate the illuminance_vertical illuminance. 

	Parameters
	----------
	img : 2D array
        Fisheye image
        
	Returns
	-------
	Eh : float
		illuminance_vertical illuminance in mlx
  
	"""
    pass

#------------------------------------------------------------------------------#

#Read in the image center coordinates
C = pd.read_csv('../Calibration/imagecenter.csv',index_col=0)
xc = C['Xcenter'][p.camera]
yc = C['Ycenter'][p.camera]
Radius = C['Radius'][p.camera]

#List of files to process
files = glob(p.data_cal+f'Light*.fit') + glob(p.data_cal+f'img*sky*.fit')

#Convert to polar coordinates
image_shape  = fits.open(files[0],uint=False)[0].data.shape
x, y = n.meshgrid(n.arange(image_shape[1]),n.arange(image_shape[0]))
R = n.sqrt((x-xc)**2 + (y-yc)**2)

#Read in or create the aperture mask 
R[n.where(R>Radius)] = n.nan 

#illuminance calculation parameters
Theta = R/Radius * n.pi/2
Omega = n.deg2rad(450/3600)**2 * n.sin(Theta)


for f in files:
    
    image  = fits.open(f,uint=False)[0].data
    
    #zenith brightness in mag per squared arcsec
    zenith_mag = zenith(image)
    print(f)
    print("Zenith:", zenith_mag)
    
    #horizontal illuminance in mlx
    horizontal_illum = illuminance_horizontal(image)
    print("Horizontal Illuminance:", horizontal_illum)
    
    
    