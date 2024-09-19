#-----------------------------------------------------------------------------#
#astrometry.py
#
#NPS Night Skies Program
#
#Last updated: 2021/01/26
#
#This script solves the image pointing and identifies the standard stars. First,
#the script crops the input into a smaller subset of images. Then, some selected
#and croped images are sent through client.py to astrometry.net to solve. The 
#detected lists of standard stars and the center coordinates are downloaded. The
#standard star lists are combined while only the central coordinates of the 
#originl image is saved. Interim processing files are deleted. 
#
#Input: 
#   (1) imagecenter.csv
#	(2) Reference fisheye images to be solved
#
#Output:
#   (1) detected_stars.csv
#   (2) center.txt
#
#History:
#	Li-Wei Hung -- Created 
#-----------------------------------------------------------------------------#
import numpy as n
import os
import pandas as pd
import re
import shutil

from astropy.io import fits
from datetime import datetime
from glob import glob
from matplotlib import pyplot as plt

# Local Source
import process_input as p     

#-----------------------------------------------------------------------------#

def main():
	"""
	This script solves the reference image and outputs a list of detected 
	standard stars and the center coordinates. See the script description for 
	detail. All the inputs are defined and passed through process_input.
	"""
	
	#Skip this script if default calibration constants will be used
	if p.measure_reference == False: 
		return
	 						
	#Read in the field of view parameters
	C = pd.read_csv(p.calibration+'imagecenter.csv',index_col=0)
	xc = C['Xcenter'][p.camera]
	yc = C['Ycenter'][p.camera]
	radius = C['Radius'][p.camera]
	scale  = 90*60*60/radius			#estimated scale ["/pix]

	#Define the cutting parameters
	nsquare = 5 #image will be cut into nsquare x nsquare pieces
	side = int(n.floor_divide(2*radius,nsquare)) #length of each side [pix]
	hs = (nsquare-2)*side/2 #half of total cropped side length
	
	#Crop and solve the reference image
	t = []
	xbound = n.arange(xc-hs,xc+hs,side, dtype=int)
	ybound = n.arange(yc-hs,yc+hs,side, dtype=int)

	ref_img = glob(p.data_cal+p.reference)[0]
	print("Reference image is: ", ref_img)
	img = fits.open(ref_img,uint=False)[0].data
	
	for i in xbound:
		for j in ybound:
		
			#skip the images at the 4 corners because the solving tends to fail 
			if i in [xbound[0],xbound[-1]] and j in [ybound[0],ybound[-1]]:
				continue
			print(i,j)
			
			#crop the image
			crop = img[j:j+side,i:i+side]
			fn = p.data_cal+'crop%i_%i.fit'%(i,j)
			fits.writeto(fn, crop, overwrite=True)
			
			#solve and astrometry; see "python client.py -help"
			cmd = f'python client.py -k {p.apikey} \
									--upload "{fn}" \
									--parity 1 \
									--scale-est {scale} \
									--corr "{fn[:-4]}_corr.fit" \
									--calibrate "{fn[:-4]}_calib.txt"'

			t1 = datetime.now()
			os.system(cmd)
			t2 = datetime.now()
			t.append(f'crop{i}_{j}.fit: {t2-t1}')
				
	print('Solving time:', t)
	
	#combine the detected reference stars in corr files from astrometry.net
	fcor = glob(p.data_cal+'*corr*')			
	D = pd.DataFrame()
	for i,f in enumerate(fcor):
		crop_shift = [int(s) for s in re.findall(r'\d+', f)[-2:]] 
		hdu = fits.open(f)
		F = pd.DataFrame.from_records(hdu[1].data).astype('float64').round(3)
		F['field_x'] += crop_shift[0] - 1 #Start counting from 0 instead of 1
		F['field_y'] += crop_shift[1] - 1 #Offset to the uncroped img position
		D = pd.concat([D,F],axis=0,join='outer')
		#copy the center RA and DEC coordinates to a new file
		if crop_shift == [xbound[int(len(xbound)/2)],ybound[int(len(xbound)/2)]]:
			shutil.copyfile(f[:-8]+'calib.txt',p.data_cal+'center.txt')
		hdu.close()
		
	D.rename(columns={'index_ra':'RA','index_dec':'DE'},inplace=True)
	D.rename(columns={'FLUX':'Flux','BACKGROUND':'Background'},inplace=True)
	D = D[['field_x','field_y','RA','DE','Flux','Background']]
	D.set_index(['RA','DE'], inplace=True)	
	D.to_csv(p.data_cal+'detected_stars.csv')
	
	#delete interim processing files
	[os.remove(f) for f in glob(p.data_cal+'*crop*')]


if __name__ == '__main__':
	main()
