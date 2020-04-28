#-----------------------------------------------------------------------------#
#wcsconversion.py
#
#NPS Night Skies Program
#
#Last updated: 2019/10/03
#
#This script is developed based on the tutorials on the astropy website
#https://docs.astropy.org/en/stable/wcs/
#
#Input: 
#   (1) 'Test_Images/Coyote_Ridge_Test_Images/Coyote_Ridge_Astrometry/wcs.fits'
#	(2) 'Test_Images/Coyote_Ridge_Astrometry/new-image.fits'
#	(3) 'Test_Images/Coyote_Ridge_Astrometry/rdls.fits'
#	(4) 'Test_Images/Coyote_Ridge_Astrometry/axy.fits'
#	(5) 'Test_Images/Coyote_Ridge_Astrometry/corr.fits'
#
#Output:
#   (1) 
#   (2) 
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import numpy as n

from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
#-----------------------------------------------------------------------------#
#wcs.fits file contains only the wcs information in the header
#Load the WCS information from a fits header, and use it to convert pixel 
#coordinates to world coordinates.

# Load the FITS hdulist using astropy.io.fits
filename = 'Test_Images/Coyote_Ridge_Astrometry/wcs.fits'
hdulist = fits.open(filename)

# Parse the WCS keywords in the primary HDU
w = WCS(hdulist[0].header)

# Three pixel coordinates of interest.
# Note we've silently assumed an NAXIS=2 image here.
# The pixel coordinates are pairs of [X, Y].
# The "origin" argument indicates whether the input coordinates
# are 0-based (as in Numpy arrays) or
# 1-based (as in the FITS convention, for example coordinates
# coming from DS9).
pixcrd = n.array([[0, 0], [24, 38], [45, 98]], dtype=n.float64)

# Convert pixel coordinates to world coordinates
# The second argument is "origin" -- in this case we're declaring we
# have 0-based (Numpy-like) coordinates.
world = w.wcs_pix2world(pixcrd, 0)
print(world)

# Convert the same coordinates back to pixel coordinates.
pixcrd2 = w.wcs_world2pix(world, 0)
print(pixcrd2)

# These should be the same as the original pixel coordinates, modulo
# some floating-point error.
assert n.max(n.abs(pixcrd - pixcrd2)) < 1e-6

#-----------------------------------------------------------------------------#
#new-image.fits is the image with the WCS header attached
#Matplotlib plots with correct WCS projection
#The WCSAxes framework, previously a standalone package, allows the WCS to be 
#used to define projections in Matplotlib. More information on using WCSAxes can be found here.

filename2 = 'Test_Images/Coyote_Ridge_Astrometry/new-image.fits'

hdu2 = fits.open(filename2)[0]
wcs = WCS(hdu2.header)

fig = plt.figure()
fig.add_subplot(111, projection=wcs)
plt.imshow(hdu2.data, origin='lower', cmap=plt.cm.viridis)
plt.xlabel('RA')
plt.ylabel('Dec')

plt.show(block=False)

#-----------------------------------------------------------------------------#
#rdls.fits contains reference stars nearby (RA,Dec table); The table is stored 
#under hdu[1].data. 
filename3 = 'Test_Images/Coyote_Ridge_Astrometry/rdls.fits'
hdu3 = fits.open(filename3)
print(hdu3.info())

#-----------------------------------------------------------------------------#
#axy.fits contains Stars detected in your images (x,y table); The table is 
#stored under hdu[1].data. 
filename4 = 'Test_Images/Coyote_Ridge_Astrometry/axy.fits'
hdu4 = fits.open(filename4)
print(hdu4.info())

#-----------------------------------------------------------------------------#
#corr.fits contains correspondences between image and reference stars (table); 
#The table is stored under hdu[1].data. "field_x", "field_y", "FLUX", and 
#"Background" are measured in the input image.
filename5 = 'Test_Images/Coyote_Ridge_Astrometry/corr.fits'
hdu5 = fits.open(filename5)
print(hdu5.info())


