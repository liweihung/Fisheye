#------------------------------------------------------------------------------#
#mask.py
#
#NPS Night Skies Program
#
#Last updated: 2023/1/6
#
#This script generates a circular mask based on the given center and radius.
#
#Input: 
#   (1) mask_input.py
#   (2) ../Calibration/imagecenter.csv
#
#Output:
#   (1) A circular mask in fits format
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#
import numpy as n

from astropy.io import fits
from matplotlib import pyplot as plt
%matplotlib qt

# Local Source
import mask_input as mi

#-----------------------------------------------------------------------------#

#Read in the field od view information
C = pd.read_csv('../Calibration/imagecenter.csv',index_col=0)
xc = C['Xcenter'][mi.camera]
yc = C['Ycenter'][mi.camera]
r = C['Radius'][mi.camera]

#Creat the mask image
x, y = n.meshgrid(n.arange(2394),n.arange(1597)) # 4x4 binning
R = n.sqrt((x-xc)**2 + (y-yc)**2)
mask = n.ones_like(R)
mask[n.where(R>r)] = n.nan

#save the bestfit model mask
hdu = fits.PrimaryHDU()
hdu.header['XCENTER'] = xc
hdu.header['YCENTER'] = yc
hdu.header['RADIUS'] = r
hdu.data = mask
hdu.writeto('../Calibration/mask_%s_%i_%i_%i.fit'%(mi.camera, xc, yc, r), overwrite=True)

#plot
fig = plt.figure()
plt.imshow(mask)
plt.show(block=False)