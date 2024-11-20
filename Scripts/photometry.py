#-----------------------------------------------------------------------------#
#photometry.py
#
#NPS Night Skies Program
#
#Last updated: 2021/01/27
#
#This script finds the bestfit zeropoint and extinction coefficient. The script
#first matches the detected stars in the image to a list of standard stars. The
#merged list contains the stardard stars' V magnitudes and the approximated x 
#and y positions in the image. Next, the zenith angle and airmass of the stars 
#are calculated based on the observing time and location. Then Gaussian PSF 
#photometry is performed to measure the detected brightness of each star. 
#Finally, the script uses different robust leanear fitting algorithems to fit 
#for the zeropoint and extintion coefficient on the M-m vs. airmass plot. 
#
#Input: 
#   (1) 'hipparcos_bright_standards_vmag6.txt' -- Standard star catalog
#	(2) 'detected_stars.csv' -- Detected stars from astrometry.net
#	(3) p.data_cal+p.reference -- reference fisheye image
#	(4) observing latitude and longitude.
#
#Output:
#   (1) p.data_cal+'zeropoint.png'
#   (2) p.data_cal+'zeropoint.csv'
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import numpy as n
import os
import pandas as pd

import astropy
# `astropy` relies on IERS (International Earth Rotation and Reference Systems Service)
# by default it will throw an error if it encounters download issues
# we override to accomodate SSL encryption
# see https://docs.astropy.org/en/stable/api/astropy.utils.iers.IERSDegradedAccuracyWarning.html
astropy.utils.iers.conf.iers_degraded_accuracy = "warn" 
from astropy.coordinates import EarthLocation
from astropy.io import fits
from astropy.time import Time
from glob import glob
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
from sklearn.linear_model import (
     TheilSenRegressor, RANSACRegressor, HuberRegressor)

# Local Source
import process_input as p     
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


def match_stars(file_std, file_detected):
	"""
	Merge the standard stars and stars detected in the images.
	
	Parameters
	----------
	file_std : str
		File name of the standard stars. Columns in the file should be 'HIP', 
		'Vmag', 'RA', 'DE', 'B-V', 'V-I'.
	file_detected : str
		File name of the stars detected in the image. Columns in the file are 
		'RA', 'DE', 'field_x', 'field_y', 'Flux', and 'Background'.
		
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
	D = pd.read_csv(file_detected,index_col=['RA','DE'])

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
		#crop image centered on the star for fitting 
		star = img[int(yc)+s[:,n.newaxis],int(xc)+s] / exptime #counts/sec
	
		#measure background
		bg = n.median(star[bw])
		
		#fit background-subtracted flux
		f = star[sw].ravel() - bg 				  #source pixels
		p0 = (xc, yc, 3, 3000)#initial parameters (x,y,std,flux)
		x, y = n.meshgrid(int(xc)+s,int(yc)+s)
		try: 
			popt = curve_fit(Gaussian_2d, [x[sw],y[sw]], f, p0=p0)[0]
		except:
			popt = n.array([xc, yc, 3, n.nan])

		#set the acceptance threshold to record the measurement
		npix = n.pi*(3*popt[2])**2 	   #numberof pixels in the aperture
		source_noise = popt[3]*exptime #noise from the source
		sky_noise = npix*bg*exptime	   #noise from the sky background
		dark_and_read_noise = 200	   #dark and read noise, approximated dummy
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
	z_err: float
		Uncertainty of the best fit intercept from OLS
	c_err: float
		Uncertainty of the best fit slope from OLS
	df : Pandas dataframe
		Dataframe containing stars used in the fitting process.
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
	
	estimators = [('Theil-Sen', TheilSenRegressor(random_state=42)),
				  ('HuberRegressor', HuberRegressor())]

	#fitting with OLS using n.polyfit
	param, cov = n.polyfit(df.Airmass, df.Vmag-df.m, 1, cov=True)
	c_err, z_err = n.sqrt(cov.diagonal())#uncertainties
	bestfit = {'OLS':[param[1], param[0]]}

	for name, estimator in estimators:
		estimator.fit(df.Airmass.values[:, n.newaxis], df.Vmag-df.m)
		if name == 'RANSAC':
			bestfit[name] = [estimator.estimator_.intercept_, 
							  estimator.estimator_.coef_[0]]
		else:
			bestfit[name] = [estimator.intercept_, estimator.coef_[0]]
	
	return bestfit, z_err, c_err, df


def main():
	"""
	This script performs Gaussian PSF photometry and fit for the photometric 
	zeropoint and extinction coefficient. See the script description for detail.
	All the inputs are defined and passed through process_input.
	"""
	
	#Skip this script if default calibration constants will be used
	if p.measure_reference == False: 
		return
	elif not os.path.exists(p.data_cal+'detected_stars.csv'):
		print('Astrometry solutions do not exist. Default zeropoint is used.')
		return

	#--------------------------------------------------------------------------#
	#		   Merge the standard star and detected star lists				   #
	#--------------------------------------------------------------------------#
	fstd = p.calibration+'standards_hipparcos_vmag6.txt' #standard 
	fcor = p.data_cal+'detected_stars.csv'  #detected from astrometry.net

	H, C, T = match_stars(fstd, fcor)

	#--------------------------------------------------------------------------#
	#						Zenith angle and airmass 						   #
	#--------------------------------------------------------------------------#
	#Open the original image and get the observing time and location
	hdu_orig = fits.open(glob(p.data_cal+p.reference)[0], fix=False)[0] 
	hdr = hdu_orig.header
	time = Time(hdr['DATE-OBS'])  		#UTC observing date and time
	c = EarthLocation(lon=hdr['SITELONG'] , lat=hdr['SITELAT'])
	lat, long = c.lat.deg, c.lon.deg 	#degrees
	
	#Zenith RA and Dec: compute based on the observing location and time
	T['ZA'], T['Airmass'] = compute_za_airmass(time, lat, long, T.RA, T.DE)
	
	#--------------------------------------------------------------------------#
	# 	 	Photometry with Gaussian PSF: calculate flux and background	 	   #
	#--------------------------------------------------------------------------#
	#photometry with Gaussian PSF
	T.field_x, T.field_y, T.Flux, T.Background, T['deltap'],T['sigma'],T['SN']=\
	photometry(T.field_x, T.field_y, hdu_orig.data, hdr['EXPTIME'])
	
	#--------------------------------------------------------------------------#
	#						Zeropoint and extinction fitting				   #
	#--------------------------------------------------------------------------#
	#fit the zeropoint and extinction
	bestfit, z_err, e_err, T_use = fit_zeropoint_and_extinction(T, 
								   selection=True, dp=1, sig=2, snr=5, z=1.5)
	T_drop = T.drop(T_use.index)
	print('Bestfit zeropoints:', bestfit)
	
	#save the bestfit to a file
	F = pd.DataFrame.from_dict(bestfit, orient='index')
	F.columns=['Zeropoint','Extinction']
	F.to_csv(p.data_cal+'zeropoint.csv')
		
	#--------------------------------------------------------------------------#
	#						Plot the fitting results						   #
	#--------------------------------------------------------------------------#
	fig = plt.figure('zeropoint_and_extinction_fit')
	
	#data points
	plt.plot(T_drop.Airmass, T_drop.Vmag-T_drop.m, 'o', c='k', fillstyle='none', 
			label='reference stars not used')
	plt.scatter(T_use.Airmass, T_use.Vmag-T_use.m, marker='o', c=T_use.SN, 
				cmap='nipy_spectral_r')
	
	#bestfit models
	colors = {'OLS':'turquoise', 'Theil-Sen':'gold', 'HuberRegressor':'black'}
	linestyle = {'OLS':'-','Theil-Sen':'-.','HuberRegressor':'--'}
	a = n.array([1,max(T.Airmass)])
	for name in bestfit:
		intercept, slope = bestfit[name]
		if name == 'OLS':
			label = r'%s: (%.2f$\pm$%.2f)x+(%.2f$\pm$%.2f)' \
					% (name,slope,e_err,intercept,z_err)
		else: 
			label = '%s: %.2fx+%.2f ' % (name,slope,intercept)
		plt.plot(a, intercept+a*slope, color=colors[name], 
				linestyle=linestyle[name], linewidth=3, label=label)
	
	#general plot setting
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('signal-to-noise ratio')
	plt.ylim(6.8,10.1)
	plt.legend(loc='upper right', fontsize=8, frameon=True, edgecolor='lightgray', numpoints=1, facecolor='white', framealpha=1)
	plt.xlabel('Airmass')
	plt.ylabel('M-m')
	plt.title('Zeropoint and Extinction Coefficient')
	imgout = p.data_cal+'zeropoint.png'
	plt.savefig(imgout, dpi=300)
	plt.show(block=False)



if __name__ == '__main__':
	main()
