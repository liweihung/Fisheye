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
#   (1) 'hipparcos_bright_standards_vmag6.txt' -- Standard star catalog
#	(2) 'Test_Images/20200427_Astrometry/corr.fits' -- matched star list
#	(3) 'Test_Images/20200427_Astrometry/wcs.fits' -- wcs coordinates
#	(4) 'Test_Images/20200427_Astrometry/40_sec_V_light_4x4_crop800.fit'
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
import re

from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS 
from glob import glob
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy import stats
from sklearn.linear_model import (
    LinearRegression, TheilSenRegressor, RANSACRegressor, HuberRegressor)
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline



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


def match_stars(file_std, files_detected):
	"""
	Merge the standard stars and stars detected in the images.
	
	Parameters
	----------
	file_std : string
		File name of the standard stars. Columns in the fie should be 'HIP', 
		'Vmag', 'RA', 'DE', 'B-V', 'V-I'.
	files_detected : list of strings
		A list of file names containing stars detected in the image. Each file
		name should have the starting x and y pixel positions as the last two 
		numbers in the file name. Each file	should have columns 'field_x',
		'field_y','index_ra','index_dec','FLUX','BACKGROUND'.
		
	Returns
	-------
	H : Pandas dataframe
		Dataframe containing the input standard stars.
	D : Pandas dataframe
		Dataframe containing the reference stars detected in the image. 
	T : Pandas dataframe
		Dataframe containing all the standard stars info in the image.
	"""
	#Standard stars
	f = n.loadtxt(file_std)
	H = pd.DataFrame(f, columns=['HIP','Vmag','RA','DE','B-V','V-I'])
	H = H.round({'RA':3,'DE':3})
	H.set_index(['RA','DE'], inplace=True)
	
	#Astrometry files containing star position and brightness in the image 
	D = pd.DataFrame()
	for i,f in enumerate(files_detected):
		crop_shift = [int(s) for s in re.findall(r'\d+', f)[-2:]] 
		hdu = fits.open(f)
		F = pd.DataFrame.from_records(hdu[1].data).astype('float64').round(3)
		F['field_x'] += crop_shift[0] - 1 #Start counting from 0 instead of 1
		F['field_y'] += crop_shift[1] - 1 #Offset to the uncroped img position
		D = pd.concat([D,F],axis=0,join='outer')
	D.rename(columns={'index_ra':'RA','index_dec':'DE'},inplace=True)
	D.rename(columns={'FLUX':'Flux','BACKGROUND':'Background'},inplace=True)
	D = D[['field_x','field_y','RA','DE','Flux','Background']]
	D.set_index(['RA','DE'], inplace=True)

	#Merge the data frames to get the overlapped stars
	T = pd.concat([H,D],axis=1,join='inner').sort_values('Vmag')
	T.reset_index(inplace=True)
	
	return H, D, T

#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
#			  Input files													  #
#-----------------------------------------------------------------------------#

#standard stars
fstd = 'hipparcos_bright_standards_vmag6.txt' 

#detected reference stars from astrometry.net
fcor = glob('Test_Images/20200427_Astrometry_small_set/*corr*') 
#fcor = glob('Test_Images/20200427_Astrometry/*800*_corr*')

#original image 
forig = 'Test_Images/20200427/40_sec_V_light_4x4.fit'

#observing location
lat, long = 40.696034, -104.600301 #degrees
elev = 1580. #meters

#-----------------------------------------------------------------------------#
#			  Merge the standard stars and astrometry star files			  #
#-----------------------------------------------------------------------------#
H, C, T = match_stars(fstd, fcor)

#-----------------------------------------------------------------------------#
#							Zenith angle and airmass 						  #
#-----------------------------------------------------------------------------#
#Compute RA and Dec at the zenith based on the observing location and time
hdu3 = fits.open(forig, fix=False)
time = Time(hdu3[0].header['DATE-OBS']) #UTC date and time
#location = EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=elev*u.m)
c0 = [time.sidereal_time('mean',longitude=long).degree, lat]#[RA,Dec] in degrees

#Zenith angle: compute and add to the dataframe
AL = DistanceAndBearing(c0[1],c0[0],T.DE,T.RA)[0]  #star altitude
T['ZA'] = n.round(90-AL, 2)  #star zenith angle

#Airmass: compute and add to the dataframe
T['Airmass'] = 1/n.cos(n.deg2rad(T['ZA'])) #airmass


#-----------------------------------------------------------------------------#
#	 Absolute photometry with Gaussian PSF: calculate flux and background	  #
#-----------------------------------------------------------------------------#
#open the original image
imgname = '40_sec_V_light_4x4.fit'
hdu3 = fits.open('Test_Images/20200427/'+imgname, fix=False)

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
	
	if delta_position<1 and sigma<2 and signal_to_noise>5:		#FIXME Adjust the SN ratio once the noise reduction is implimented!
		T.loc[i,['field_x','field_y']] = popt[0], popt[1]
		T.loc[i,['Flux','Background']] = popt[3], bg   #[counts/sec]
		T.loc[i,'fit'] = True
		T.loc[i,'SN'] = signal_to_noise
		T.loc[i,'Aperture'] = npix #background aperture radius [pix]		

#delete the rows with bad photometry measurements
print('%s stars are used in the photometric calibration.' % sum(T.fit==True))
print('Additional %s stars rejected due to bad photometry.' % sum(T.fit==False))
T.drop(T[T.fit==False].index, inplace=True)
T.drop(columns='fit', inplace=True)
T = T.round({'field_x':2,'field_y':2,'Flux':2,'Background':2,'SN':2,'Sigma':2})


#-----------------------------------------------------------------------------#
#						Zeropoint and extinction fitting					  #
#-----------------------------------------------------------------------------#
#apparent magnitude of the background-subtraced stellar flux [counts/sec]
T['m'] = -2.5*n.log10(T.Flux)
T.dropna(inplace=True)

T_orig = T.copy()

#T.drop(T[T.Vmag-T.m > 8.3].index, inplace=True)
#T.drop(T[T.Vmag-T.m < 7.8].index, inplace=True)
T.drop(T[n.abs(stats.zscore(T.Vmag-T.m)) > 1].index, inplace=True)

T_dropped = T_orig.drop(T.index)

#fitting for zeropoint and extinction
param, cov = n.polyfit(T.Airmass, T.Vmag-T.m, 1, cov=True)
c, z = param                         #bestfit coefficient and zeropoint
c_err, z_err = n.sqrt(cov.diagonal())#uncertainties

#---------------------------

#plotting
fig = plt.figure('zeropoint3')
#---------------------------
estimators = [('OLS', LinearRegression()),
              ('Theil-Sen', TheilSenRegressor(random_state=42)),
              ('RANSAC', RANSACRegressor(random_state=42)),
              ('HuberRegressor', HuberRegressor())]
colors = {'OLS': 'turquoise', 'Theil-Sen': 'gold', 'RANSAC': 'lightgreen', 'HuberRegressor': 'black'}
linestyle = {'OLS': '-', 'Theil-Sen': '-.', 'RANSAC': '--', 'HuberRegressor': '--'}
lw = 3

x_plot = n.linspace(T_orig.Airmass.min(), T_orig.Airmass.max())
for name, estimator in estimators:
	print(name)
	model = make_pipeline(PolynomialFeatures(1), estimator)
	model.fit(T_orig.Airmass[:, n.newaxis], T_orig.Vmag-T_orig.m)
	#mse = mean_squared_error(model.predict(X_test), y_test)
	y_plot = model.predict(x_plot[:, n.newaxis])
	plt.plot(x_plot, y_plot, color=colors[name], linestyle=linestyle[name],
			linewidth=lw, label='%s: error = %.3f' % (name, 999))

legend_title = 'Error of Mean\nAbsolute Deviation\nto Non-corrupt Data'
legend = plt.legend(loc='upper right', frameon=False, title=legend_title,
                    prop=dict(size='x-small'))

#---------------------------
#plt.plot(T.Airmass, T.Vmag-T.m, 'o', alpha=T.SN/T.SN.max(), label='Hipparcos standard stars')
plt.plot(T_dropped.Airmass, T_dropped.Vmag-T_dropped.m, 'o', c='k', 
		 fillstyle='none', label='reference stars not used')
plt.scatter(T.Airmass, T.Vmag-T.m, marker='o' ,c=T.SN, cmap='nipy_spectral_r')
#plt.plot(T_orig.Airmass, T_orig.Vmag-T_orig.m, 'ko', alpha=0.1, label='reference stars not used')

plt.colorbar()

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

#star zenith angle

J = pd.concat([H,C],axis=1).dropna(subset=['Flux']).sort_values('Flux')
J.reset_index(inplace=True)
J['ZA'] = n.round(90-DistanceAndBearing(c0[1],c0[0],J.DE,J.RA)[0], 2)


#import warnings
#warnings.filterwarnings("ignore")
#hdu3 = fits.open('Test_Images/20200427_Astrometry/wcs.fits')








