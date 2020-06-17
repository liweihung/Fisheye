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
import operator
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
	file_std : str
		File name of the standard stars. Columns in the file should be 'HIP', 
		'Vmag', 'RA', 'DE', 'B-V', 'V-I'.
	files_detected : list of str
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
		Dataframe containing all the matched standard stars info in the image.
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


def compute_za_airmass(time, latitude, longitude, ra, dec):
	"""
	Compute zenith angles and airmass of the input list of stars given their 
	RA and Dec from the specified time and observing location.
	
	Parameters
	----------
	time : Time object from astropy.time
		UTC observing date and time.
	latitude : number
		Latitude in degrees of the observing site.
	longitude : number
		Longitude in degrees of the observing site.
	ra : list of numbers
		List of the right ascension in degrees.
	dec : list of numbers
		List of the declination in degrees.
		
	Returns
	-------
	za : list of numbers
		A list of zenith angle in degrees.
	airmass : list of numbers
		A list of airmass.
	"""
	
	#Zenith RA and Dec: compute based on the observing location and time
	zenith_ra = time.sidereal_time('mean',longitude=longitude).degree #[degree]
	zenith_de = latitude											  #[degree]

	#Zenith angle: compute and add to the dataframe
	al = DistanceAndBearing(zenith_de,zenith_ra,dec,ra)[0]  #star altitude
	za = n.round(90-al, 2)  #star zenith angle [degree]

	#Airmass: compute and add to the dataframe
	airmass = 1/n.cos(n.deg2rad(za)) #airmass
	
	return za, airmass


def photometry(x, y, img, exptime, sa=5, bai=6, bao=10):
	"""
	Photometry of the stars with Gaussian PSF
	
	Parameters
	----------
	x : list of numbers
		X coordinates of star location on the image.
	y : list of numbers
		Y coordinates of star location on the image. 
	img : 2d array
		Input image of stars to have photometry performed on.
	exptime : number
		Exposure time in second.
	sa : int, optional 
		Source window radius in pixels.
	bai : int, optional
		Bacuground ring inner radius in pixels.
	bao : int, optional
		Bacuground ring outer radius in pixels.
		
	Returns
	-------
	bestfit_x : list of numbers
		Best fit x coordinates of the stars.
	bestfit_y : list of numbers
		Best fit y coordinates of the stars.
	flux : list of numbers
		Best fit flux of the background-subtracted stars in counts per second, 
		assuming Gaussian PSF.
	background : list of numbers
		Best fit background around the stars in counts per second.
	delta_position : list of numbers
		Difference in pixels between the bestfit position and initial guess.
	sigma : list of numbers
		Sigma of the 2d gaussian in pixels.
	SN : list of numbers
		Estimated signal-to-noise ration of the stars.
	"""
	#define fitting window size 
	s = n.arange(-bao, bao+1) 
	xi, yi = n.meshgrid(s, s)
	r = n.sqrt(xi**2+yi**2)
	sw = n.where(r<=sa)               #source fitting window
	bw = n.where((r>bai) & (r<bao))   #background window
	
	bestfit_x, bestfit_y, flux, background = [], [], [], [] 
	delta_position, sigma, SN = [], [], []
	for xc, yc in zip(x,y):	
		#cropp image centered on the star for fitting 
		star = img[int(yc)+s[:,n.newaxis],int(xc)+s] / exptime #counts/sec
	
		#measure background
		bg = n.median(star[bw])
		
		#fit background-subtracted flux
		f = star[sw].ravel() - bg 				  #source pixels
		p0 = (xc, yc, 3, 3000)#initial parameters (x,y,std,flux)
		x, y = n.meshgrid(int(xc)+s,int(yc)+s)
		popt = curve_fit(Gaussian_2d, [x[sw],y[sw]], f, p0=p0)[0]

		#set the acceptance threshold to record the measurement

		npix = n.pi*(3*popt[2])**2 	   #numberof pixels in the approximated aperture
		source_noise = popt[3]*exptime #noise from the source
		sky_noise = npix*bg*exptime	   #noise from the sky background
		dark_and_read_noise = 0		   #dark and read noise http://www2.lowell.edu/rsch/LMI/ETCMethod.pdf
									#FIXME http://www.astro.wisc.edu/~sheinis/ast500/AY500_lect5.ppt.pdf	
		totalnoise = n.sqrt(source_noise+sky_noise+dark_and_read_noise)
		signal_to_noise = popt[3]*exptime/totalnoise       						
		
		#record the fitting results
		bestfit_x.append(round(popt[0],2))
		bestfit_y.append(round(popt[1],2))
		flux.append(round(popt[3],2))   #[counts/sec]
		background.append(round(bg,2))  #[counts/sec]
		delta_position.append(n.sum(((popt-p0)[0:2]**2))**0.5)  #[pix]
		sigma.append(abs(popt[2]))          				    #[pix]
		SN.append(round(signal_to_noise,2))
		
	return bestfit_x, bestfit_y, flux, background, delta_position, sigma, SN
	

def fit_zeropoint_and_extinction(df, selection=True, dp=1, sig=2, snr=5, z=1):
	"""
	Find the bestfit of zeropoint and extinction.

	Parameters
	----------
	df : Pandas dataframe
		Pandas dataframe containing the airmass and photometry of star. Must 
		have df['Airmass'] and df['Flux'] columns.
	selection : bool
		Use the specified criteria (dp, sig, and snr) to select photometry 
		points for fitting. If True, df must have 'deltap', 'sigma', and 'SN'
		columns.
	dp : number, optional
		Position shift threshold in pixels for accepting the photometry point.
	sig : number, optional
		Gaussian sigma threshold in pixels for accepting the photometry point.
	snr : number, optional
		Signal-to-noise threshold for accepting the photometry point.
	z : number, optional
		Zscore threshold for accepting the photometry point.
		
	Returns
	-------
	bestfit : dict
		Bestfit parameters formated as {'modelname':[intercept,slope]}
	df_drop : Pandas dataframe
		Contain rows of rejected stars.
	"""
	
	#apparent magnitude of the background-subtraced stellar flux [counts/sec]
	df['m'] = -2.5*n.log10(df.Flux)
	df.dropna(inplace=True)
	ntotal = len(df)

	#use only selected points for fitting
	if selection:
		df_drop = df[(df['deltap']>dp) | (df['sigma']>sig) | (df['SN']<snr) |\
					 (n.abs(stats.zscore(df.Vmag-df.m)) > z)]
		df = df.drop(df_drop.index)
		
	nselect = len(df)
	nreject = ntotal-len(df)

	print('Total: %s stars have photometric measurements.' %ntotal)
	print('Used: %s (%s%%) for fitting.'%(nselect, round(100*nselect/ntotal)))
	print('Rejected: %s (%s%%).' %(nreject, round(100*nreject/ntotal)))
	
	##fitting for zeropoint and extinction using n.polyfit
	#param, cov = n.polyfit(df.Airmass, df.Vmag-df.m, 1, cov=True)
	#c, z = param                         #bestfit coefficient and zeropoint
	#c_err, z_err = n.sqrt(cov.diagonal())#uncertainties
	
	estimators = [('OLS', LinearRegression()),
				  ('Theil-Sen', TheilSenRegressor(random_state=42)),
				  ('RANSAC', RANSACRegressor(random_state=42)),
				  ('HuberRegressor', HuberRegressor())]

	bestfit = {}
	for name, estimator in estimators:
		estimator.fit(df.Airmass[:, n.newaxis], df.Vmag-df.m)
		if name == 'RANSAC':
			bestfit[name] = [estimator.estimator_.intercept_, 
							  estimator.estimator_.coef_]
		else:
			bestfit[name] = [estimator.intercept_, estimator.coef_]
	
	return bestfit, df_drop
	

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
#						Zenith angle and airmass 							  #
#-----------------------------------------------------------------------------#
#Open the original image and get the observing time
hdu_orig = fits.open(forig, fix=False)[0] #open the original image file
time = Time(hdu_orig.header['DATE-OBS'])  #UTC observing date and time

#Zenith RA and Dec: compute based on the observing location and time
T['ZA'], T['Airmass'] = compute_za_airmass(time, lat, long, T.RA, T.DE)


#-----------------------------------------------------------------------------#
#	   	 Photometry with Gaussian PSF: calculate flux and background	 	  #
#-----------------------------------------------------------------------------#
#photometry with Gaussian PSF
T.field_x, T.field_y, T.Flux, T.Background, T['deltap'], T['sigma'], T['SN'] =\
photometry(T.field_x, T.field_y, hdu_orig.data, hdu_orig.header['EXPTIME'])

#-----------------------------------------------------------------------------#
#						Zeropoint and extinction fitting					  #
#-----------------------------------------------------------------------------#
#fit the zeropoint and extinction
bestfit, T_drop = fit_zeropoint_and_extinction(T, selection=True, 
											   dp=1, sig=1, snr=5, z=1.5)
T_use = T.drop(T_drop.index)
	
#-----------------------------------------------------------------------------#
#							Plot the fitting results						  #
#-----------------------------------------------------------------------------#
fig = plt.figure('zeropoint_and_extinction_fit')

#data points
plt.plot(T_drop.Airmass, T_drop.Vmag-T_drop.m, 'o', c='k', fillstyle='none', 
		 label='reference stars not used')
plt.scatter(T_use.Airmass, T_use.Vmag-T_use.m, marker='o', c=T_use.SN, 
			cmap='nipy_spectral_r')

#bestfit models
colors = {'OLS': 'turquoise', 'Theil-Sen': 'gold', 'RANSAC': 'lightgreen', 
		  'HuberRegressor': 'black'}
linestyle = {'OLS':'-', 'Theil-Sen':'-.', 'RANSAC':'--', 'HuberRegressor':'--'}
a = n.array([1,max(T.Airmass)])
for name in bestfit:
	intercept, slope = bestfit[name]
	plt.plot(a, intercept+a*slope, color=colors[name], linestyle=linestyle[name],
		linewidth=3, label='%s: %.2fx+%.3f ' % (name,slope,intercept))

#general plot setting
plt.legend(loc='upper right', frameon=False, numpoints=1)
plt.colorbar()
plt.xlabel('Airmass',fontsize=14)
plt.ylabel('M-m',fontsize=14)
plt.title('Zeropoint and Extinction Coefficient',fontsize=14)
#imgout = filepath.calibdata+dnight+'/extinction_fit_%s_%s.png' %(filter,s[0])
#plt.savefig(imgout)
#plt.close('zeropoint')
plt.show(block=False)


#----------------old code below------------------------------------------------#

#T_orig = T.copy()
#T.drop(T[T.Vmag-T.m > 8.3].index, inplace=True)
#T.drop(T[T.Vmag-T.m < 7.8].index, inplace=True)
#T_dropped = T_orig.drop(T.index)

#plt.plot(T_orig.Airmass, T_orig.Vmag-T_orig.m, 'ko', alpha=0.1, label='reference stars not used')
#plt.errorbar(0,z,z_err,fmt='o',label='zeropoint: %.3f+-%.3f'%(z,z_err))

#location = EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=elev*u.m)
#c0 = [time.sidereal_time('mean',longitude=long).degree, lat]#[RA,Dec] in degrees






