#-----------------------------------------------------------------------------#
#reduction.py
#
#NPS Night Skies Program
#
#Last updated: 2020/08/25
#
#This script performs basic image reduction and applies mask on the fisheye 
#image data collected by the NPS Night Skies Program. The script corrects for:
#	(1) Bias 
#	(2) Dark 
#	(3) Flat
#   (4) Linearity response of the detector
#
#Note: All the calibration files should be taken at the same temperature as the 
#science images. 
#
#Input: 
#   (1) Raw field data -- dark, bias, and science frames
#   (2) Flat
#   (3) Linearity curves
#	(4) Mask
#
#Output:
#   (1) Reduced science images
#
#History:
#	Li-Wei Hung -- Created in 2020
#
#------------------------------------------------------------------------------#
import numpy as n
import os

from astropy.io import fits
from glob import glob

# Local Source
import filepath     

#------------------------------------------------------------------------------#
#Linearity - load the linearity curve (ADU, multiplying factor)
xp, fp = n.loadtxt(filepath.linearity, unpack=True)


#Flat - load in flat
flat = fits.open(filepath.flatV,uint=False)[0].data


#Bias - generate averaged bias 
biaslist = [
	fits.open(i,uint=False)[0].data 
	for i in glob(filepath.data_raw+'*bias*.fit')]
bias = n.average(biaslist,axis=0)


#Dark - average dark that is bias subtracted and linearity response corrected
darklist = [
	fits.open(i,uint=False)[0].data 
	for i in glob(filepath.data_raw+'*dark*.fit')]
dark_bias_subtracted = n.average(darklist,axis=0)-bias     
dark = dark_bias_subtracted * n.interp(dark_bias_subtracted,xp,fp) 


#Mask - read in the fisheye mask
mask = fits.open(filepath.mask,uint=False)[0].data


#Output - create calibrated images folder
if not os.path.isdir(filepath.data_cal): os.makedirs(filepath.data_cal)


#Calibration - calibrate the science images
for i in glob(filepath.data_raw+'*light*.fit'):
	image  = fits.open(i,uint=False)[0]
	light  = image.data 						#science image
	light *= mask								#apply fisheye mask
	light -= bias								#subtract bias
	light *= n.interp(light,xp,fp)				#correct for linearity
	light -= dark								#subtract dark
	light /= flat								#divide by flat
	fits.writeto(filepath.data_cal+i[len(filepath.data_raw):],light,
				 header=image.header,overwrite=1)
