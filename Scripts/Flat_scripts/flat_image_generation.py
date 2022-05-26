#------------------------------------------------------------------------------#
#flat_image_generation.py
#
#NPS Night Skies Program
#
#Last updated: 2021/12/20
#
#Given the input of a set of images with randomly orientated and evenly 
#illumiated strips, the script automaticall detect the strip orientation and 
#select the pixels that are colse to the center line. These pixels are then used 
#to determine the flat profile as the distance to the center of the lens. The 
#flat profile is a parabola in the center + two lines + a constant.  
#
#Input: 
#   (1) a set of flat strip image in flat_image_generation_input.flatstrips
#	(2) flat_image_generation_input with three break points
#
#Output:
#   (1) radial profile.png
#	(2) 2D flat model.png
#	(3) 2D flat model in fits file
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#

import numpy as n

from astropy.io import fits
from glob import glob
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm

# Local Source
import flat_image_generation_input as fi

#------------------------------------------------------------------------------#
#read in all the images
img_folder = fi.calibration + fi.flatstrips
imglist = glob(img_folder+'*.fit')
hdr = fits.open(imglist[0])[0].header
imgcube = n.empty((len(imglist),hdr['NAXIS2'],hdr['NAXIS1']))

print('Reading in %i images in %s' %(len(imglist),img_folder))
for i,f in tqdm(enumerate(imglist)):
	imgcube[i] = fits.open(f,uint=False)[0].data
	

#find the center of the fisheye and define the radial distance grid R
sum = n.sum(imgcube,axis=0)
center = n.unravel_index(n.argmax(sum, axis=None), sum.shape)
print('The fisheye image center:', center) #(778, 1160)
X,Y = n.meshgrid(n.arange(hdr['NAXIS1']),n.arange(hdr['NAXIS2'])) 
R = n.sqrt((X-center[1])**2+(Y-center[0])**2)
	
#set the detection rows and columns to find two ends of the stript
row1,col1,row2,col2 = 400, 800, 1200, 1600

#initialize r and l
r = n.array([]) # pixel distance from the fisheye center 
l = n.array([]) # brightness

#identify the pixels near the center line of the strips
print('Detecting %i flatstrip orientations...' %len(imgcube))
for i,img in tqdm(enumerate(imgcube)):

	#finding two points (x1,y1) and (x2,y2) on the ends of a flat stript. 
	h = max(img[row1])
	v = max(img[:,col1])
	h2 = max(img[row2])
	v2 = max(img[:,col2])
	if h>v:
		y1 = row1
		x1 = n.median(n.where(img[row1]>0.9*h)[0]).astype('int')
		y2 = row2
		x2 = n.median(n.where(img[row2]>0.9*h2)[0]).astype('int')
	else:
		y1 = n.median(n.where(img[:,col1]>0.9*v)[0]).astype('int')
		x1 = col1
		y2 = n.median(n.where(img[:,col2]>0.9*v2)[0]).astype('int')
		x2 = col2
		
	#calculating the distance d of each pixel to the line Ax+By+C = 0 
	A = (y1-y2)/(x1-x2)
	B = -1
	C = y1-A*x1
	d = abs(A*X+B*Y+C)/n.sqrt(A**2+B**2)		
	
	#recording the pixels near the center line
	w = n.where(d<3)
	r = n.hstack((r,R[w]))
	l = n.hstack((l,img[w]))
	
print('Fitting the flat profile...')
#flat profile: r<r1 -- parabola; r1<r<r2 -- linear; r>r2 -- nan
[r1,r2,r3] = fi.linear_region

#fit a 2nd degree polynomial + two linear lines + a constant
def func3(x, a, b, c, d, e, f, g):
	crosspt1 = a*r1**2+b*r1+c
	crosspt2 = d*r2+e
	crosspt3 = f*r3+g

	m12 = (crosspt1-crosspt2)/(r1-r2)
	m23 = (crosspt2-crosspt3)/(r2-r3)

	y = n.zeros(len(x))
	y += (a*x**2+b*x+c) * (x<=r1)  
	y += (m12*x-m12*r2+crosspt2) * (x>r1) * (x<=r2)
	y += (m23*x-m23*r3+crosspt3) * (x>r2) * (x<=r3)
	y +=  crosspt3* (x>r3)
	return y

popt, pcov = curve_fit(func3, r, l)
print(popt)
xp = n.linspace(r.min(), r.max(), 200)


#plot the radial profile of the flat and the best fits
fig = plt.figure(1)
plt.plot(r,l,'.')
plt.plot(xp, func3(xp, *popt), 'y-', label='break points=%5i, %5i, and %5i' % (r1,r2,r3))
plt.title('Radial profile of the flat and the best fits')
plt.legend()
plt.xlabel('Radial distance (pix)')
plt.ylabel('Brightness (ADU)')
plt.savefig(fi.calibration+fi.flatstrips+'radial profile.png')

print('Generating the flat model image')	
#generate a flat model image
model = n.empty_like(R)
model[n.where(R<=r3)] = func3(R[n.where(R<=r3)], *popt)
model[n.where(R>r3)] = n.nan
model = model/func3(n.array([0.,]), *popt)	#normalize the center to 1

#save the bestfit flat model as a fits file
hdu = fits.PrimaryHDU()
hdu.header['BREAKPT1'] = r1
hdu.header['BREAKPT2'] = r2
hdu.header['BREAKPT3'] = r3
hdu.data = model
hdu.writeto(fi.calibration+fi.flatstrips[:-1]+'.fit', overwrite=True)

#plot the model image 
fig = plt.figure(2)
plt.imshow(model)
plt.colorbar()
plt.title('2D flat model')
plt.savefig(fi.calibration+fi.flatstrips+'2D flat model.png')
plt.show(block=False)
