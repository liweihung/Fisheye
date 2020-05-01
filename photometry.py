#-----------------------------------------------------------------------------#
#photometry.py
#
#NPS Night Skies Program
#
#Last updated: 2020/04/15
#
#This script performs absolute photometry on fisheye images  
#
#
#Input: 
#   (1) 'hipparcos_bright_standards.txt' -- Standard star catalog
#	(2) 'Test_Images/20200427_Astrometry/corr.fits' -- matched star list
#	(3) 'Test_Images/20200427_Astrometry/wcs.fits' -- wcs coordinates
#	(4) 'Test_Images/20200427_Astrometry/30_sec_V_light_2x2_crop1300.fit'
#
#	(2) 'Test_Images/20200427_Astrometry/new-image.fits'
#	(3) 'Test_Images/20200427_Astrometry/rdls.fits'	
#	(5) 
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
import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS 
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

# Local Source
#from gaussian import Gaussian_2d
from sphericalgeometry import DistanceAndBearing

#-----------------------------------------------------------------------------#
def Gaussian_2d(xy, x0, y0, sigma, V):
	'''
	This module returns the (x,y) value of the 2D gaussian function with the 
	given parameters. V is the volume under the curve. 
	'''
	x, y = xy
	g = V/(2*n.pi*sigma**2)*n.exp(-((x-x0)**2+(y-y0)**2)/(2*sigma**2))
	return g.ravel()
	
#-----------------------------------------------------------------------------#
#			  Merge the standard stars and astrometry star files			  #
#-----------------------------------------------------------------------------#
#Standard stars
hips = n.loadtxt('hipparcos_bright_standards.txt')
H = pd.DataFrame(hips, columns=['HIP','Vmag','RA','DE','B-V','V-I'])
H = H.round({'RA':2,'DE':2})
H.set_index(['RA','DE'], inplace=True)

#Astrometry file contains correspondences between image and reference stars 
hdu = fits.open('Test_Images/20200427_Astrometry/corr.fits')
C = pd.DataFrame.from_records(hdu[1].data).astype('float64').round(2)
C.rename(columns={'index_ra':'RA','index_dec':'DE'},inplace=True)
C.rename(columns={'FLUX':'Flux','BACKGROUND':'Background'},inplace=True)
C = C[['field_x','field_y','RA','DE','Flux','Background']]
C.set_index(['RA','DE'], inplace=True)

#Merge the data frames to get the overlapped stars
T = pd.concat([H,C],axis=1,join='inner').sort_values('Vmag')
T.reset_index(inplace=True)

#Start counting from 0 instead of 1
T.field_x -= 1
T.field_y -= 1


#-----------------------------------------------------------------------------#
#							Zenith angle and airmass 						  #
#-----------------------------------------------------------------------------#
#Read in the WCS coordinates
hdu2 = fits.open('Test_Images/20200427_Astrometry/wcs.fits', fix=False)
w = WCS(hdu2[0].header)

#Convert the central pixel to world coordinates w/ 0-based coordinates
imgw = hdu2[0].header['IMAGEW']
imgh = hdu2[0].header['IMAGEH']
c0 = w.wcs_pix2world([[imgw/2-0.5,imgh/2-0.5],],0)[0] #[RA,Dec] in degrees

#Zenith angle: compute and add to the dataframe
AL = DistanceAndBearing(c0[1],c0[0],T.DE,T.RA)[0]  #star altitude
T['ZA'] = n.round(90-AL, 2)  #star zenith angle

#Airmass: compute and add to the dataframe
T['Airmass'] = 1/n.cos(n.deg2rad(T['ZA'])) #airmass


#-----------------------------------------------------------------------------#
#	 Absolute photometry with Gaussian PSF: calculate flux and background	  #
#-----------------------------------------------------------------------------#
#open the original image
imgname = '30_sec_V_light_2x2_crop1300.fit'
hdu3 = fits.open('Test_Images/20200427_Astrometry/'+imgname, fix=False)
data = hdu3[0].data
exptime = hdu3[0].header['EXPTIME']

#set aperture radius
sa = 5   	#[pix] source fitting aperture radius
bai = 6  	#[pix] bacuground aperture ring inner radius
bao = 10 	#[pix] bacuground aperture ring outer radius

#define window size for fitting
s = n.arange(-bao, bao+1) 
x, y = n.meshgrid(s, s)
r = n.sqrt(x**2+y**2)
sw = n.where(r<=sa)               #source window
bw = n.where((r>bai) & (r<bao))   #background window

T.loc[:, 'fit'] = False
for i in range(len(T)):	
	#star location on the image
	x_center = int(T.field_x[i])
	y_center = int(T.field_y[i])
	
	#cropped image centered on the star for fitting 
	star = data[y_center+s[:,n.newaxis],x_center+s] / exptime #counts per second
	
	#measuring background
	bg = n.median(star[bw])
	
	#fitting background-subtracted flux
	f = star[sw].ravel() - bg 				  #source pixels
	p0 = (T.field_x[i], T.field_y[i], 3, 3000)#initial parameters (x,y,std,flux)
	x, y = n.meshgrid(x_center+s,y_center+s)
	popt = curve_fit(Gaussian_2d, [x[sw],y[sw]], f, p0=p0)[0]
		
	#set the acceptance threshold to record the measurement
	delta_position = n.sum(((popt-p0)[0:2]**2))**0.5   #position diff [pix]
	sigma = abs(popt[2])          				       #sigma of the 2d gaussian
	npix = n.pi*(3*sigma)**2 	   #numberof pixels in the approximated aperture
	source_noise = popt[3]*exptime #noise from the source
	sky_noise = npix*bg*exptime	   #noise from the sky background
	dark_and_read_noise = 0		   #dark and read noise http://www2.lowell.edu/rsch/LMI/ETCMethod.pdf
								   #FIXME http://www.astro.wisc.edu/~sheinis/ast500/AY500_lect5.ppt.pdf	
	totalnoise = n.sqrt(source_noise+sky_noise+dark_and_read_noise)
	signal_to_noise = popt[3]*exptime/totalnoise       						
	
	if delta_position<1 and sigma<2 and signal_to_noise>60:		#FIXME Adjust the SN ratio once the noise reduction is implimented!
		T.loc[i,['field_x','field_y']] = popt[0], popt[1]
		T.loc[i,['Flux','Background']] = popt[3], bg   #[counts/sec]
		T.loc[i,'fit'] = True
		T.loc[i,'SN'] = signal_to_noise
		T.loc[i,'Aperture'] = npix #background aperture radius [pix]		

#delete the rows with bad photometry measurements
print('%s stars are used in the photometric calibration.' % sum(T.fit==True))
print('%s stars rejected due to bad photometry.' % sum(T.fit==False))
T.drop(T[T.fit==False].index, inplace=True)
T.drop(columns='fit', inplace=True)
T = T.round({'field_x':2,'field_y':2,'Flux':2,'Background':2,'SN':2,'Sigma':2})


#-----------------------------------------------------------------------------#
#						Zeropoint and extinction fitting					  #
#-----------------------------------------------------------------------------#
#apparent magnitude of the background-subtraced stellar flux [counts/sec]
T['m'] = -2.5*n.log10(T.Flux)
T.dropna(inplace=True)

#fitting for zeropoint and extinction
param, cov = n.polyfit(T.Airmass, T.Vmag-T.m, 1, cov=True)
c, z = param                         #bestfit coefficient and zeropoint
c_err, z_err = n.sqrt(cov.diagonal())#uncertainties

#plotting
fig = plt.figure('zeropoint')
plt.plot(T.Airmass, T.Vmag-T.m, 'o', label='Hipparcos standard stars')
a = n.arange(n.ceil(max(T.Airmass))+1)
plt.plot(a,c*a+z,'-',lw=2,label='Best fit: %.2fx+%.3f' %(c,z))
plt.errorbar(0,z,z_err,fmt='o',label='zeropoint: %.3f+-%.3f'%(z,z_err))
plt.legend(loc=0, numpoints=1)
plt.xlabel('Airmass',fontsize=14)
plt.ylabel('M-m',fontsize=14)
plt.title('Zeropoint and Extinction Coefficient',fontsize=14)
#imgout = filepath.calibdata+dnight+'/extinction_fit_%s_%s.png' %(filter,s[0])
#plt.savefig(imgout)
#plt.close('zeropoint')
plt.show(block=False)


#-----------------------------------------------------------------------------#
#import warnings
#warnings.filterwarnings("ignore")
#hdu3 = fits.open('Test_Images/20200427_Astrometry/wcs.fits')








