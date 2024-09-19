#------------------------------------------------------------------------------#
#linearity.py
#
#NPS Night Skies Program
#
#Last updated: 2023/2/6
#
#Given a set or sets of images taken with a range of exposure time, this script 
#characterizes the detector's linearity response. If multiple brightness 
#settings were used, only images taken with the same brightness setting should 
#be save in the same folder. 
#
#Input: 
#   (1) linearity_input.py
#   (2) Paths to all folders of images taken with different brightness settings.
#       Each folder contains a set of images with different exposure time.
#   (3) Output filename
#
#Output:
#   (1) linearity multiplier txt file
#	(2) linearity multiplier plot
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
from scipy.interpolate import UnivariateSpline

# Local Source
import linearity_input as L

#------------------------------------------------------------------------------#

#Read in the center of the field of view
C = pd.read_csv('../Calibration/imagecenter.csv',index_col=0)
xc = C['Xcenter'][L.camera]
yc = C['Ycenter'][L.camera]
r = 30 

#Define the image region for characterizing linearity
x, y = n.meshgrid(n.arange(2392),n.arange(1596)) # 4x4 binning
R = n.sqrt((x-xc)**2 + (y-yc)**2)

#Determine bias
biasimg = fits.open(glob(L.calibration + L.imagefolder[0]+'*.fit')[0])[0].data
bias = n.median(biasimg[n.where(R>800)])
print('Midian bias is ', bias)

#------------------------------------------------------------------------------#
#                           Reading in images                                  #
#------------------------------------------------------------------------------#
dflist = []

for i in range(len(L.imagefolder)): 
    
    #read in all the images in a folder
    img_folder = L.calibration + L.imagefolder[i]
    linearitylist = []
    for f in glob(img_folder+'*.fit'):	
        image   = fits.open(f,uint=False)[0]
        light   = image.data[n.where(R<r)]		#science 
        median  = n.median(light) - bias        #bias-subtracted
        exptime = image.header['EXPTIME']
        if median < (65535-bias):               #exclude saturated cases
            linearitylist.append([exptime, median])

    #compute multiplier for linearity
    df = pd.DataFrame(linearitylist, columns=['Exptime','Count'])
    df['Rate'] = df['Count'] / df['Exptime']
    df.sort_values('Count', inplace=True)
    df.reset_index(drop=True, inplace=True)
    dflist.append(df)

#Adjust for different brightness settings used 
for i in range(1,len(dflist)):
    b = dflist[i].iloc[0].Rate / dflist[i-1].iloc[-1].Rate
    dflist[i].Rate /= b

#Combine the dataframes 
df_combined = pd.concat(dflist)
df_combined.sort_values('Count',inplace=True)
df_combined.reset_index(drop=True, inplace=True)

#Groupby count bins
c = n.array([0,10,20,40,50,75,100,125,175,200,350,500,1000,1500,2000,5000,9500,
             10000,15000,20000,30000,40000,50000,55000,60000,62500,64500,65000,
             65535])
c_middle = n.array([(c[i]+c[i+1])/2 for i in range(len(c)-1)])
df2 = df_combined.groupby(pd.cut(df_combined['Count'],c_middle)).mean()
df = pd.DataFrame({'Count':df2.Count,'Rate':df2.Rate})

#Add the first and last points
df = df.append({'Count':0,'Rate':df_combined.Rate[0]},ignore_index=True)
df = df.append({'Count':65535,'Rate':df_combined.Rate.iat[-1]},ignore_index=True)
df.sort_values('Count',inplace=True)
df.reset_index(drop=True, inplace=True)

#Compute the linearity multiplier
df['Multiplier'] = 1 / df['Rate'] 
df['Multiplier'] *= df.loc[df['Count']>=10000,"Rate"].values[0] # normalize
df = df.round({'Multiplier':3})

#------------------------------------------------------------------------------#
#                            Linearity Model                                   #
#------------------------------------------------------------------------------#

#fit a spline function for counts between 10 and 25000
spl = UnivariateSpline(df.Count[(df.Count>10)&(df.Count<25000)], 
                       df.Multiplier[(df.Count>10)&(df.Count<25000)], 
                       bbox=[10, 25000], s=8e-5)

def linearity_model(x):
    y = n.interp(x,df.Count, df.Multiplier)
    w = n.where((x>10)&(x<25000))
    y[w] = spl(x[w])
    return y.round(3)


#Save the linearity multiplier
linearity_multiplier = linearity_model(c)
df3 = pd.DataFrame({'Count':c, 'Multiplier':linearity_multiplier})
df3.to_csv(L.calibration+L.outfilename, sep='\t', header=False, index=False)


#plot
plt.plot(df.Count, df.Multiplier, 'o')
plt.plot(c, linearity_multiplier, '-')
plt.xlabel('Pixel Value [counts]')
plt.ylabel('Linearity Multiplier')
plt.savefig(L.calibration+L.outfilename[:-4]+'.png', dpi=300)
plt.show(block='False')

