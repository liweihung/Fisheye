#------------------------------------------------------------------------------#
#flat.py
#
#NPS Night Skies Program
#
#Last updated: 2022/11/30
#
#Given a set of images with randomly orientated and evenly illumiated strips, 
#the script automatically detects the strip orientation and select the pixels 
#that are colse to the center line. These pixels are then used 
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
import pandas as pd

from astropy.io import fits
from glob import glob
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from tqdm import tqdm

# Local Source
import flat_input as fi

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                      Reading in flat strips images                           #
#------------------------------------------------------------------------------#

#read in all the images
img_folder = fi.calibration + fi.flatstrips
imglist = glob(img_folder+'*.fit')
hdr = fits.open(imglist[0])[0].header
imgcube = n.empty((len(imglist),hdr['NAXIS2'],hdr['NAXIS1']))

print('Reading in %i images in %s' %(len(imglist),img_folder))
for i,f in tqdm(enumerate(imglist)):
	imgcube[i] = fits.open(f,uint=False)[0].data
	

#find the center of the strips
sum = n.sum(imgcube,axis=0)
center = n.unravel_index(n.argmax(sum, axis=None), sum.shape)
dark = n.median(n.mean(imgcube,axis=0)[:,0:500])
print('The fisheye image center:', center)

#define the radial distance grid R
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

l = l-dark #dark subtraction
#------------------------------------------------------------------------------#
#                             Flat model fitting                               #
#------------------------------------------------------------------------------#
print('Fitting the flat profile...')

#flat profile: r<r1 -- spline; r1<r<r2 -- linear; r>r2 -- constant
[r1,r2] = fi.linear_region

#fit a spline function in r < r1
df = pd.DataFrame({'r': r[n.where(r<r1)], 'l': l[n.where(r<r1)]})
bins = df.groupby(df.r.round(0)).mean()
spl = UnivariateSpline(bins.index, bins.l, s=2.5e7)

#fit linear lines in r >= r1
def func(x, a, b):
	crosspt1 = spl(r1)
	crosspt2 = a*r2+b
	m = (crosspt1-crosspt2)/(r1-r2)
	y = n.zeros(len(x))
	y += (m*x-m*r2+crosspt2) * (x>=r1) * (x<r2)
	y +=  crosspt2* (x>=r2)
	return y

popt, pcov = curve_fit(func, r[n.where(r<r2+200)], l[n.where(r<r2+200)])

#combined model curve
def modelf(x):
	y = func(x, *popt) + spl(x)*(x<r1)
	return y 

#------------------------------------------------------------------------------#
#                             Generate Flat model                              #
#------------------------------------------------------------------------------#
print('Generating the flat model image')	
#generate a flat model image
R2 = n.sqrt((X-center[1])**2+(Y-center[0])**2)
model = n.empty_like(R2)
model[n.where(R2<r2)] = modelf(R2[n.where(R2<r2)]) 
model[n.where(R2>=r2)] = n.nan
model = model/modelf(n.array([0.,]))	#normalize the center to 1

#save the bestfit flat model as a fits file
hdu = fits.PrimaryHDU()
hdu.header['XCENTER'] = center[1]
hdu.header['YCENTER'] = center[0]
hdu.header['BREAKPT1'] = r1
hdu.header['BREAKPT2'] = r2
hdu.data = model
hdu.writeto(fi.calibration+fi.flatstrips[:-1]+'.fit', overwrite=True)


#------------------------------------------------------------------------------#
#                                  Plotting                                    #
#------------------------------------------------------------------------------#
#xp = n.linspace(r.min(), r2+100, 1000)
xp = n.linspace(r.min(), r1, 100)

#plot the radial profile of the flat and the best fits
fig = plt.figure(1)
#plt.plot(r,l,'.', label='strip center at (%i,%i)' %(center))
#plt.plot(xp, modelf(xp), 'y-', lw=2, label='break points at %3i and %3i' %(r1,r2))
plt.plot(r,l,'.', label='illuminated pixels')
plt.plot(xp, modelf(xp), 'y-', lw=2, label='1-D smoothing spline fit')
plt.xlim(0,900)
plt.tick_params(axis='both', which='major', labelsize=13)
#plt.title('Radial profile of the flat and the best fits')
plt.legend(fontsize=13)
plt.xlabel('Radial Distance (pixel)', size=13)
plt.ylabel('Brightness (ADU)', size=13)
plt.savefig(fi.calibration+fi.flatstrips+'radial profile.png', dpi=300)


#plot the model image 
fig = plt.figure(2)
ax = plt.gca()
im = ax.imshow(model)
ax.tick_params(axis='both', which='major', labelsize=13)
plt.xlabel('X (pixel)', size=13)
plt.ylabel('Y (pixel)', size=13)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)

plt.colorbar(im, cax=cax)
#plt.title('2D flat model')

plt.savefig(fi.calibration+fi.flatstrips+'2D flat model.png', dpi=300)

plt.show(block=False)
