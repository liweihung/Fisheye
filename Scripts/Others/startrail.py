#-----------------------------------------------------------------------------#
#startrail.py
#
#NPS Night Skies Program
#
#Last updated: 2022/04/11
#
#This script all the input images to make a star trail image. 
#
#Input: 
#   (1) Raw field data -- sky images
#
#Output:
#   (1) Star trail image
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#
import numpy as n

from astropy.io import fits
from glob import glob

# Local Source
import process_input as p  

#------------------------------------------------------------------------------#
#Raw sky images 
skyimg = [
	fits.open(i,uint=False)[0].data 
	for i in glob(p.data_raw+'*-s*')]
startrail = n.sum(skyimg,axis=0)
fits.writeto(p.data_raw+'startrail.fit',startrail,overwrite=1)