#-----------------------------------------------------------------------------#
#reduction.py
#
#NPS Night Skies Program
#
#Last updated: 2021/01/22
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
#   (3) Linearity curve
#	(4) Mask
#
#Output:
#   (1) Reduced science images
#
#History:
#	Li-Wei Hung -- Created in 2020
#	Li-Wei Hung -- Modified the code to support processing V and B images 
#
#------------------------------------------------------------------------------#
import numpy as n
import os

from astropy.io import fits
from glob import glob

# Local Source
import process_input as p     

#------------------------------------------------------------------------------#
def main():
	"""
	This script performs basic image reduction. See the description on the top 
	of this script. All the inputs are defined and passed through process_input.
	"""
	
	#Linearity - load the linearity curve (ADU, multiplying factor)
	xp, fp = n.loadtxt(p.linearity, unpack=True)
	
	
	#Flat - read in the flat filenames
	flatfile = {'V':p.flatv, 'B':p.flatb}

	
	#Bias - generate averaged bias 
	biaslist = [
		fits.open(i,uint=False)[0].data 
		for i in glob(p.data_raw+'Bias*.fit')]
	bias = n.average(biaslist,axis=0)
	
	
	#Dark - average dark is bias subtracted and linearity response corrected
	darklist = [
		fits.open(i,uint=False)[0].data 
		for i in glob(p.data_raw+'Dark*.fit')]
	dark_bias_subtracted = n.average(darklist,axis=0)-bias     
	dark = dark_bias_subtracted * n.interp(dark_bias_subtracted,xp,fp) 
	
	
	#Mask - read in the fisheye mask
	mask = fits.open(p.mask,uint=False)[0].data
	
	
	#Output - create calibrated images folder
	if not os.path.isdir(p.data_cal): os.makedirs(p.data_cal)
	
	
	#Calibration - calibrate the science images
	for k in flatfile.keys(): 							#loop through filters
		
		try: 
			flat = fits.open(flatfile[k],uint=0)[0].data
			flat /= n.nanmean(flat)
		except: continue
		
		for f in glob(p.data_raw+f'Light*{k}*.fit'):
		
			image  = fits.open(f,uint=False)[0]
			light  = image.data 					#science image
			light *= mask							#apply fisheye mask
			light -= bias							#subtract bias
			light *= n.interp(light,xp,fp)			#correct for linearity
			light -= dark							#subtract dark
			light  = n.clip(light,0.1,n.inf)		#replace negatives w/ 1
			light /= flat							#divide by flat
			
			hdr = image.header
			hdr['PROCESS'] = p.processor
			hdr['history'] = f'Linearity curve used is {p.fn_linearity}'
			hdr['history'] = f'Flat used is {flatfile[k][len(p.calibration):]}'
			hdr['history'] = f'Mask used is {p.fn_mask}'			
			fits.writeto(p.data_cal+f[len(p.data_raw):],light,
						 header=hdr,overwrite=1)

						 
if __name__ == '__main__':
	main()
