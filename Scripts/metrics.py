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

#Read in the image center coordinates
C = pd.read_csv('../Calibration/imagecenter.csv',index_col=0)
xc = C['Xcenter'][p.camera]
yc = C['Ycenter'][p.camera]

#Set aperture diameter to about 1 degree 
r = 4.75 #radius in pixel; 1 pix = 379 arcsec 

#Creat the aperture mask 
x, y = n.meshgrid(n.arange(2392),n.arange(1596)) # 4x4 binning
R = n.sqrt((x-xc)**2 + (y-yc)**2)
mask = n.ones_like(R)
mask[n.where(R>r)] = n.nan

for f in glob(p.data_cal+f'Light*.fit'):
    image  = fits.open(f,uint=False)[0].data
    zenith = n.median(image[n.where(R<r)])
    print(zenith)
    