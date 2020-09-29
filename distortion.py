#-----------------------------------------------------------------------------#
#distortion.py
#
#NPS Night Skies Program
#
#Last updated: 2020/09/28
#
#This script charaterizes the radial and angular distortions of the 
#fisheye image.  
#
#
#Input: 
#   (1) original fisheye image
#	(2) detected_Stars.csv
#	(3) mask
#	(4) center.txt
#
#Output:
#   (1) 
#
#History:
#	Li-Wei Hung -- Created in 2020
#
#------------------------------------------------------------------------------#
import ast

import numpy as n
import pandas as pd
import seaborn as sns

from astropy.io import fits
from astropy.time import Time
from matplotlib import pyplot as plt
from scipy import stats
from sklearn.linear_model import LinearRegression

# Local Source
import filepath
from sphericalgeometry import DistanceAndBearing

#------------------------------------------------------------------------------#

S = pd.read_csv(filepath.data_cal+'detected_stars.csv') #RA, DEC, X, Y of stars

#observing location
#lat, long = 40.696034, -104.600301 #degrees

#-----------------------------------------------------------------------------#
#						Zenith angle and azimuth 							  #
#-----------------------------------------------------------------------------#
#Open the original image and get the observing time
forig = filepath.data_raw+'40_sec_V_light_4x4.fit' #original image 
hdu_orig = fits.open(forig, fix=False)[0] #open the original image file

#Read in the coordinates of the image center
file = open(filepath.data_cal+'center.txt', "r")
center = ast.literal_eval(file.read())
file.close()
center_ra = center['ra']
center_de = center['dec']

al, S['Azimuth'] = DistanceAndBearing(center_de,center_ra,S.DE,S.RA)
S['ZA'] = 90-al

#Mask - read in the fisheye mask to find the center of the view
mask   = fits.open(filepath.mask,uint=False)[0].data
yc, xc = n.round(n.mean(n.where(mask==1),axis=1))	#center of fisheye view

#Radial distance
S['RD'] = n.sqrt((S['field_x']-xc)**2+(S['field_y']-yc)**2)

#Get coeffs of linear fit
lr_fi_false = LinearRegression(fit_intercept=False)
X = n.array(S['RD']).reshape(-1,1)
lr_fi_false.fit(X,S['ZA'])
lr_fi_false_yhat = n.dot(X, lr_fi_false.coef_) + lr_fi_false.intercept_

#------------------- First plot -------------------------------#

#plot
#plt.figure()
ax1 = plt.subplot2grid((5, 1), (0, 0), rowspan=4)
plt.scatter(S['RD'],S['ZA'], label='Reference stars')
plt.plot(S['RD'], lr_fi_false_yhat, 'r-', label='fit_intercept=False')

plt.legend()
plt.ylabel('Zenith Angle [degree]')


#------------------- Plot the residuals -------------------------------#

plt.subplot2grid((25, 1), (21, 0), rowspan=5, sharex=ax1)
plt.plot([S['RD'].min(),S['RD'].max()],n.zeros(2), color= 'black',linestyle ='-', lw=1)
plt.plot(S['RD'],S['ZA']-lr_fi_false_yhat,'.')
plt.ylabel('Residual')
plt.xlabel('Distance from the center [pix]')


#------------------- Second plot -------------------------------#
#Angle 
S['Angle'] = n.rad2deg(n.arctan2(S['field_y']-yc,S['field_x']-xc))+180
S['Azimuth'] = (S['Azimuth']+180-12.3)%360
slope, intercept, r, p, std_err = stats.linregress(S['Angle'],S['Azimuth'])

plt.figure()
ax2 = plt.subplot2grid((5, 1), (0, 0), rowspan=4)

k = 'y = {0:.2f}x - {1:.2f} \nr$^2$ = {2:.2f}; p = {3:.2e}'\
	.format(slope,abs(intercept),r**2,p)
sns.regplot(S['Angle'],S['Azimuth'], ci=68, marker='o', line_kws={'label':k})

plt.legend()
plt.ylabel('True Azimuth [degree]')

#------------------- Plot the residuals -------------------------------#

plt.subplot2grid((25, 1), (21, 0), rowspan=5, sharex=ax2)
plt.plot([S['Angle'].min(),S['Angle'].max()],n.zeros(2), color= 'black',linestyle ='-', lw=1)
y_residual = S['Azimuth']-(S['Angle']*slope+intercept)
plt.plot(S['Angle'],y_residual,'.')
plt.ylabel('Residual')
plt.xlabel('Angle in the image [degree]')



plt.show(block=False)



