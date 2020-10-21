#-----------------------------------------------------------------------------#
#main.py
#
#NPS Night Skies Program
#
#Last updated: 2020/10/13
#
#This script 
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
#import numpy as n
#import os

#from astropy.io import fits
#from glob import glob

# Local Source
import astrometry
import filepath
import median_filter
import photometric_calibration
import photometry
import reduction

#------------------------------------------------------------------------------#
reduction.main()

if filepath.measure_reference:
	astrometry.main()
	photometry.main()
	
photometric_calibration.main()
median_filter.main()
	
	
