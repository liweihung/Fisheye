#-----------------------------------------------------------------------------#
#astrometry.py
#
#NPS Night Skies Program
#
#Last updated: 2020/09/22
#
#This script solves the image pointing and identifies the standard stars in the 
#images. First, the script crops the input into a smaller subset of images. 
#Then, some selected croped images are sent through client.py to astrometry.net 
#to solve. The detected standard stars list and the center coordinates are 
#downloaded. The standard star lists re combined while only the central 
#coordinates of the originl imageis save. Interim processing files are deleted. 
#
#Input: 
#   (1) Mask for the fisheye view
#	(2) Raw fisheye images to be solved
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
import filepath     

#-----------------------------------------------------------------------------#

#Mask - read in the fisheye mask to find the center of the view
mask   = fits.open(filepath.mask,uint=False)[0].data
yc, xc = n.round(n.mean(n.where(mask==1),axis=1))	#center of fisheye view
xmin   = n.min(n.where(mask==1),axis=1)[1] 			#x edge of fisheye view
scale  = 90*60*60/(xc-xmin) 						#estimated scale ["/pix]

#Define the cutting parameters
nsquare = 5 #image will be cut into nsquare x nsquare pieces
side = n.int(n.floor_divide(2*(xc-xmin),nsquare)) #length of each side in pixel
hs = (nsquare-2)*side/2 #half of total cropped side length

#Loop through all the fisheye images to crop and solve
t = []
xbound = n.arange(xc-hs,xc+hs,side, dtype=n.int)
ybound = n.arange(yc-hs,yc+hs,side, dtype=n.int)
for f in glob(filepath.data_raw+'*light*.fit'):
	img = fits.open(f,uint=False)[0].data 	#science image
	for i in xbound:
		for j in ybound:
		
			#skip the images at the 4 corners because the solving tends to fail 
			if i in [xbound[0],xbound[-1]] and j in [ybound[0],ybound[-1]]:
				continue
			print(i,j)
			
			#crop the image
			crop = img[j:j+side,i:i+side]
			fn = filepath.data_cal+f[len(filepath.data_raw):-4]+'_crop%i_%i.fit'%(i,j)
			fits.writeto(fn, crop, overwrite=True)
			
			#solve and astrometry; see "python client.py -help"
			cmd = f"python client.py -k {filepath.apikey} \
									 --upload {fn} \
									 --parity 1 \
									 --scale-est {scale} \
									 --corr {fn[:-4]}_corr.fit \
									 --calibrate {fn[:-4]}_calib.txt"
			t1 = datetime.now()
			os.system(cmd)
			t2 = datetime.now()
			t.append([fn,f'Total time:{t2-t1}'])
			
print(t)

#combine the detected reference stars in corr files from astrometry.net
fcor = glob(filepath.data_cal+'*corr*')			
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
		shutil.copyfile(f[:-8]+'calib.txt',filepath.data_cal+'center.txt')
	hdu.close()
	
D.rename(columns={'index_ra':'RA','index_dec':'DE'},inplace=True)
D.rename(columns={'FLUX':'Flux','BACKGROUND':'Background'},inplace=True)
D = D[['field_x','field_y','RA','DE','Flux','Background']]
D.set_index(['RA','DE'], inplace=True)	
D.to_csv(filepath.data_cal+'detected_stars.csv')

#delete interim processing files
[os.remove(f) for f in glob(filepath.data_cal+'*crop*')]
