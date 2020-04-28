#-----------------------------------------------------------------------------#
#astrometry_results.py
#
#NPS Night Skies Program
#
#Last updated: 2020/04/07
#
#This script plots the astrometry outputs for a series of test images
#
#Input: 
#   (1) 'Test_Images/Coyote_Ridge_Astrometry/solved_results.xlsx'
#
#Output:
#   (1) 'Test_Images/Coyote_Ridge_Astrometry/solved_results.png'
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import numpy as n
import pandas as pd

#from astropy.io import fits
#from astropy.wcs import WCS
from matplotlib import pyplot as plt
#-----------------------------------------------------------------------------#

column_names = ['width','time','RA','Dec','scale','orientation']
df = pd.read_excel('Test_Images/Coyote_Ridge_Astrometry/solved_results.xlsx', names=column_names)
df['orientation'] = (df['orientation']+360)%360 #gets rid of negative values

dfc = df.subtract(df.mean()) #subtract the column mean

#plt.plot(df['width'], dfc['time'], 'o-', label='processing time (s)')
plt.plot(df['width'], dfc['RA'], 'o-', label='RA (degree)')
plt.plot(df['width'], dfc['Dec'], 'o-', label='Dec (degree)')
plt.plot(df['width'], dfc['scale'], 'o-', label='pixel scale (arcsec/pix)')
plt.plot(df['width'], dfc['orientation'], 'o-', label='orientation (up, degrees E of N)')

plt.legend(loc=0)
plt.xlabel('Image Width (pix)')
plt.ylabel('Deviation from Mean')
plt.savefig('Test_Images/Coyote_Ridge_Astrometry/solved_results.png')
plt.show(block=False)