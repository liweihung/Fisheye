#-----------------------------------------------------------------------------#
#standard.py
#
#NPS Night Skies Program
#
#Last updated: 2019/09/17
#
#This script has a funtion to gerenerate the map of the brightest (Vmag < 3.5) 
#Hipparcos standard stars at the give time and location.
#
#Input: 
#   (1) 'hipparcos_bright_standards.txt'
#
#Output:
#   (1) Standard stars maps standard_north.png, standard_south.png
#   (2) Example location-specific standard star map standard_chelan.png
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as n
import pandas as pd

#-----------------------------------------------------------------------------#
#read in the standard star catalog
hips = n.loadtxt('hipparcos_bright_standards.txt')
H = pd.DataFrame(hips, columns=['HIP','Vmag','RA','DE','B-V','V-I'])
H.HIP = H.HIP.astype('int64')
print(H.info())

#-----------------------------------------------------------------------------#
utcoffset = {'PT':-8*u.hour, 'MT':-7*u.hour, 'ET':-6*u.hour}

def StandardStarMap(lat,lon,ele,time,timetype='UTC',timezone=None,dls=True):
	"""
	Gerenerate the map of the standard stars at the give time and location.
	lat = latitude in degress
	lon = longitude in degrees, west is negative
	ele = elevation in meters
	time = UTC date and time; ie. '2019-6-01 2:00:00'
	timetype = 'UTC' or 'local'
	timezone = 'PT','MT',or 'ET'; Pacific Time, Mountain Time, Eatern Time
	dls = True or False; daylightsavings
	"""
	site = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=ele*u.m)
	if timetype != 'UTC':
		if timetype == 'local':
			time = Time(time)-utcoffset[timezone]-dls*u.hour #UTC time
		else:
			print('timetype input error: timetype must be "UTC" or "local"')
	
	for i in H.index:
		star = SkyCoord(ra=H.loc[i,'RA']*u.degree, dec=H.loc[i,'DE']*u.degree, frame='icrs')
		staraltaz = star.transform_to(AltAz(obstime=time,location=site))
		H.loc[i,'AZ'] = staraltaz.az.value
		H.loc[i,'ALT'] = staraltaz.alt.value
		
	Hv = H[H.ALT>=0].copy() # Stars above the horizon
	
	return Hv
	

#-----------------------------------------------------------------------------#

if __name__ == "__main__":

	#------Plot the standard stars in the northern and southern hemisphere----#
	
	Hn = H[H.DE>=0] # Stars with positive Declination
	Hs = H[H.DE<0]  # Stars with negative Declination
	
	#plot the standard stars in polar coordinate
	fig = plt.figure() # Stars with positive Declination
	ax = plt.subplot(111, projection='polar')
	c = ax.scatter(n.deg2rad(Hn.RA), 90-Hn.DE, c=Hn.Vmag, cmap='jet', alpha=0.75)
	ax.set_rlim(0,90)
	ax.set_yticks(n.arange(0, 91, 10))
	ax.set_yticklabels(ax.get_yticks()[::-1])
	#ax.set_theta_zero_location('S') 
	ax.set_theta_direction(-1)  # theta increasing clockwise
	ax.grid(True)
	cbar = fig.colorbar(c, ax=ax)
	cbar.set_clim(-1.5,3.5)
	plt.savefig('standard_north.png', dpi=200)
	
	fig2 = plt.figure() # Stars with negative Declination
	ax2 = plt.subplot(111, projection='polar')
	c2 = ax2.scatter(n.deg2rad(Hs.RA), Hs.DE, c=Hs.Vmag, cmap='jet', alpha=0.75)
	ax2.set_rlim(-90,0)
	ax2.set_theta_direction(-1)  # theta increasing clockwise
	ax2.grid(True)
	cbar2 = fig2.colorbar(c2, ax=ax2)
	cbar2.set_clim(-1.5,3.5)
	plt.savefig('standard_south.png', dpi=200)
	
	#------Case Study---------------------------------------------------------#

	#Example: Standard stars on the sky in Chelan simulated for the CHLN190601 data 
	Hv = StandardStarMap(47.53991,-120.35551,1251,'2019-6-01 2:00:00','local','PT')
	
	#plot the standard stars in polar coordinate
	fig = plt.figure() # Stars above the horizon
	ax = plt.subplot(111, projection='polar')
	c = ax.scatter(n.deg2rad(Hv.AZ), 90-Hv.ALT, c=Hv.Vmag, cmap='jet', alpha=0.75)
	ax.set_rlim(0,90)
	ax.set_yticks(n.arange(0, 91, 10))
	ax.set_yticklabels(ax.get_yticks()[::-1])
	ax.set_theta_zero_location('N') 
	ax.grid(True)
	cbar = fig.colorbar(c, ax=ax)
	cbar.set_clim(-1.5,3.5)
	plt.savefig('standard_chelan.png', dpi=200)
	
	plt.show(block=False)
