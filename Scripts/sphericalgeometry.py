#-----------------------------------------------------------------------------#
#sphericalgeometry.py
#
#NPS Night Skies Program
#
#Last updated: 2019/09/17
#
#This script has a funtion to calculate the great circle angular distances 
#between two points on a sphere.
#
#History:
#	Li-Wei Hung -- Created
#
#-----------------------------------------------------------------------------#
import numpy as n

#-----------------------------------------------------------------------------#

def DistanceAndBearing(alt, az, alt_array, az_array):
	"""
	This function calculates the great circle angular distances between two 
	points on a sphere. The formulas are adopted from Bullock, R. 
	"Great Circle Distances and Bearings Between Two Locations." 
	https://dtcenter.org/met/users/docs/write_ups/gc_simple.pdf
	
	Input:
	alt = altitude in degrees of the new origin
	az = azimuth in degrees of the new origin
	alt_array = array of altitude values to be transformed to the new coordinate
	az_array = array of azimuth values to be transformed to the new coordinate
	
	Output:
	a = array of new angular distance (altitude)
	b = array of new bearing (azimuth)
	"""
	
	#convert the input angles from degrees to radians
	p = n.deg2rad(alt)
	l = n.deg2rad(az)
	pa = n.deg2rad(alt_array)
	la = n.deg2rad(az_array)
	
	#calculate the great circle angular distances
	d = n.rad2deg(n.arccos(n.sin(p)*n.sin(pa)+n.cos(p)*n.cos(pa)*n.cos(l-la)))
	a = 90 - d
	
	#calculate the bearing
	S = n.cos(pa)*n.sin(l-la)
	C = n.cos(p)*n.sin(pa)-n.sin(p)*n.cos(pa)*n.cos(l-la)
	t = n.rad2deg(n.arctan2(S,C))
	b = t+az-180 #correction for the azimuthal angle offset
	
	return a, b
	
	
if __name__ == '__main__':

	#Glacier Basin Campground
	from astropy.io import fits
	imgpath = 'C:/Users/lhung/Pictures/20190823_ROMO _astronomy/'
	img = fits.open(imgpath+'ASICAP_2019-08-23_21_39_22_281.FIT')[0]
	time = img.header['DATE-OBS'] #UTC
	lat = 40.325979
	lon = -105.598273
	elevation = 2591
	
	#Find standard stars above the horizon
	from standard import StandardStarMap
	Hv = StandardStarMap(lat,lon,elevation,time) #HIP,Vmag,RA,DE,B-V,V-I,AZ,ALT
	
	import matplotlib.pyplot as plt
	def PolarPlot(az, alt):
		"""
		Plot the given azimuth and altitude points with the correct polar setting
		az = azimuth in degrees
		alt = altitude in degrees
		"""
		fig = plt.figure()
		ax = plt.subplot(111, projection='polar')
		c = ax.scatter(n.deg2rad(az), 90-alt, c=Hv.Vmag, cmap='jet', alpha=0.75)
		ax.set_rlim(0,90)
		ax.set_yticks(n.linspace(0, 90, 10))
		ax.set_yticklabels(ax.get_yticks()[::-1])
		ax.set_theta_zero_location('N') 
		ax.grid(True)
		cbar = fig.colorbar(c, ax=ax)
		cbar.set_clim(-1.5,3.5)

	#Plot the standard stars in polar coordinate
	PolarPlot(Hv.AZ, Hv.ALT)
	plt.savefig('standard_glacierbasin.png', dpi=200)
	
	#Test the sperical distance and bearing calculation
	lat_neworigin = 30.
	lon_neworigin = 350.
	
	ALT_new, AZ_new = DistanceAndBearing(lat_neworigin,lon_neworigin,Hv.ALT,Hv.AZ)
	
	#plot the standard stars with the new origin in polar coordinate
	PolarPlot(AZ_new, ALT_new) 
	plt.savefig('standard_glacierbasin_offset.png', dpi=200)
	
	plt.show(block=False)


