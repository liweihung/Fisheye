#-----------------------------------------------------------------------------#
#solve.py
#
#NPS Night Skies Program
#
#Last updated: 2020/08/19
#
#This script uses the API (client.py) to upload images to astrometry.net for 
#solving the images. 
#
#Input: 
#   (1) Fits images
#
#Output:
#   (1) Fits images with wcs attached
#   (2) WCS files
#	(3) Corr files - correspondences between image and reference stars (table)
#
#History:
#	Li-Wei Hung -- Created 
#-----------------------------------------------------------------------------#
import os
from datetime import datetime
from glob import glob, iglob

#-----------------------------------------------------------------------------#

filedir = "Test_Images/20200427_Astrometry_small_set/failed/"

for i, fn in enumerate(glob(filedir+"*")):
	cmd = f"python client.py -k your_API_key --upload {fn} --newfits {fn}_new.FIT --wcs {fn}_wcs.FIT --corr {fn}_corr.FIT"
	t1 = datetime.now()
	os.system(cmd)
	t2 = datetime.now()
	print(i,fn,f'Total time:{t2-t1}')

#cmd2 = "python client.py -help"


