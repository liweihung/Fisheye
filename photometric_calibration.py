#-----------------------------------------------------------------------------#
#photometric_calibration.py
#
#NPS Night Skies Program
#
#Last updated: 2020/10/14
#
#This script 
#Input: 
#   (1) 
#
#Output:
#   (1) 
#   (2) 
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import ast
import numpy as n
import pandas as pd

from astropy.io import fits
from glob import glob

from matplotlib import pyplot as plt


# Local Source

import colormaps

import filepath

#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#			  Input files													  #
#-----------------------------------------------------------------------------#
#original image 

#read in zeropoint
if filepath.use_default_zeropoint:
	Z = pd.read_csv(filepath.calibration+'zeropoint.csv',index_col=0)
	zp = Z['Zeropoint'][filepath.camera]
else:
	Z = pd.read_csv(filepath.data_cal+'zeropoint.csv')
	zp = Z['Zeropoint'].mean()
	
#Read in the pixel scale
file = open(filepath.data_cal+'center.txt', "r")
center = ast.literal_eval(file.read())
file.close()
pixscale = center['pixscale'] #["/pix]

#platescale adjustment
psa = 2.5*n.log10(pixscale**2) 

#brightness calibration
for f in glob(filepath.data_cal+'*light*'):
	image = fits.open(f,uint=False)[0]
	exptime = image.header['EXPTIME'] #[sec]
	data = image.data
	mc = zp+psa-2.5*n.log10(data/exptime)


plt.figure()
plt.rcParams['image.cmap'] = 'NPS_mag'
current_cmap = plt.cm.get_cmap()
current_cmap.set_bad(color='black')
plt.imshow(mc,vmin=14,vmax=24)

plt.colorbar()
plt.show(block=False)






