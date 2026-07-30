#------------------------------------------------------------------------------#
#metrics.py
#
#NPS Night Skies Program
#
#Given the calibrated fisheye images, the script extracts various night sky 
#brightness metrics 
#
#Input: 
#   (1) 
#	(2) 
#
#Output:
#   (1) 
#	(2) 
#	(3) 
#
#History:
#	Li-Wei Hung -- Created 
#
#------------------------------------------------------------------------------#
import importlib
import numpy as n
import os
import pandas as pd
import shutil

from astropy.io import fits
from datetime import datetime
from glob import glob

# Local Source
#import process_input as p     
import projection

#------------------------------------------------------------------------------#

def zenith(img, a=4.75):
    """
	Calculate the median zenith brightness. 

	Parameters
	----------
	img : 2D array
        Fisheye image
        
	a : number, optional
        radius in pixel; 1 pix = 379 arcsec; aperture diameter should be set to
        about 1 degree, comparable to CCD zenith aperture 
		
	Returns
	-------
	zenith_mag : float
		median zenith brightness in mag per sequared arcsec
  
	"""
    zenith_mag = round(n.median(img[n.where(r<a)]), 2)
    return zenith_mag


def build_illuminance_geometry(xc, yc, radius, image_shape):
	"""
	Precompute the zenith angle and per-pixel solid angle grids needed for
	illuminance integration, assuming an equidistant fisheye projection
	(zenith angle proportional to pixel radius, as used in projection.py:
	r_deg = 90 * r_normalized).

	Parameters
	----------
	xc, yc : float
		Fisheye image center coordinates [pix].
	radius : float
		Fisheye field-of-view radius in pixels (90 degrees zenith angle).
	image_shape : tuple
		(height, width) of the image.

	Returns
	-------
	Theta : 2D array
		Zenith angle [radians] at each pixel. NaN beyond the FoV.
	Omega : 2D array
		Solid angle subtended by each pixel [steradians]. NaN beyond the FoV.
	"""
	ny, nx = image_shape
	x, y = n.meshgrid(n.arange(nx), n.arange(ny))

	#pixel radius and zenith angle (equidistant fisheye: theta linear in r)
	Rpix = n.sqrt((x-xc)**2 + (y-yc)**2)
	Theta = (n.pi/2) * (Rpix/radius) 	#[radians], 0 at zenith

	#mask out everything beyond the 90-degree FoV
	Theta[Rpix > radius] = n.nan

	#angular scale: radians of zenith angle per pixel of image radius
	dtheta_dr = (n.pi/2) / radius 		#constant for equidistant fisheye

	#Per-pixel solid angle for an equidistant projection:
	#   dOmega = sin(theta) * dtheta * dphi
	#For a pixel at radius Rpix, azimuthal extent dphi corresponds to an arc
	#length of 1 pixel: dphi = 1/Rpix (radians), and radial extent dtheta =
	#dtheta_dr (radians per pixel). So the pixel's solid angle is obtained
	#via the Jacobian: dOmega/dA_pix = sin(theta) * dtheta_dr / Rpix
	#(valid everywhere except Rpix -> 0, handled separately below)
	with n.errstate(divide='ignore', invalid='ignore'):
		Omega = n.sin(Theta) * dtheta_dr / Rpix

	#special-case the zenith pixel(s) where Rpix -> 0 (avoid singularity);
	#the solid angle density -> dtheta_dr^2 in the limit (finite), since
	#sin(theta)/Rpix = sin(dtheta_dr*Rpix)/Rpix -> dtheta_dr as Rpix -> 0
	near_zenith = Rpix < 1.0
	Omega[near_zenith] = dtheta_dr**2

	Omega[Rpix > radius] = n.nan

	return Theta, Omega



def illuminance_horizontal(img):
    """
	Calculate the horizontal illuminance. 

	Parameters
	----------
	img : 2D array
        Fisheye image
        
	Returns
	-------
	E_h : float
		horizontal illuminance in mlx
  
	"""

    L = 108.48*n.exp(20.7233-0.92104*img) # [ucd m-2], Dan's conversion formula
    dE = L*n.cos(Theta)*dOmega 
    E_h = round(n.nansum(dE)/1000, 2) #Horizontal illuminance [mlx]

    return E_h


def illuminance_vertical(img):
    """
	Calculate the illuminance_vertical illuminance. 

	Parameters
	----------
	img : 2D array
        Fisheye image
        
	Returns
	-------
	Eh : float
		illuminance_vertical illuminance in mlx
  
	"""
    pass

#------------------------------------------------------------------------------#



def parse_lat_lon(coord):
    if isinstance(coord, float) or isinstance(coord, int):
        # Already a number, return as-is
        return float(coord)
    
    # Otherwise, assume it's a string
    parts = coord.strip().split()
    if len(parts) == 3:
        d, m, s = map(float, parts)
        return d + m / 60 + s / 3600
    else:
        return float(coord)


#------------------------------------------------------------------------------#


#List of files to process
datasets = glob('../Data_processed/WHSA_20220227/process_input.py')
#    '../Data_processed/GUIS_20231210B/process_input.py') #+ \
#			  glob('../Data_processed/EGLI*/process_input.py') + \
#             glob('../Data_processed/BOSE*/process_input.py')

# List to store each row of data
data_rows = []


for pi in datasets:
    shutil.copy(pi, 'process_input_metric.py')
    import process_input_metric as p
    importlib.reload(p)
    
    #Read in the image center coordinates
    C = pd.read_csv('../Calibration/imagecenter.csv',index_col=0)
    xc = C['Xcenter'][p.camera]
    yc = C['Ycenter'][p.camera]
    R  = C['Radius'][p.camera]
    
    #Read in the terrain mask
    with fits.open(p.data_cal+f'mask.fit', uint=False, memmap=False) as hdul:
        mask = hdul[0].data.astype(float, copy=True)

    #Define the image geometry in polar coordinates
    image_shape = mask.shape
    x, y = n.meshgrid(n.arange(image_shape[1]),n.arange(image_shape[0]))
    r = n.sqrt((x-xc)**2 + (y-yc)**2)

    #Define horizontal illuminance parameters
    dTheta = n.pi/2/R
    dPhi = 1/r
    Theta = n.pi/2*r/R
    dOmega = dTheta * dPhi * n.sin(Theta)
    dOmega[r==0] = dTheta**2 #special case for zenith pixel(s)

    #List of images to process        
    images = glob(p.data_cal+f'Light*.fit') + glob(p.data_cal+f'img*sky*.fit')

    for file in images:

        print(f"Processing {file} ...")
        #Read in the image and apply the terrian mask
        with fits.open(file, uint=False, memmap=False) as hdul:
            img = hdul[0].data.astype(float, copy=True) * mask
        
        #Zenith brightness [mag per squared arcsec]
        zenith_mag = zenith(img)
        print("Zenith:", zenith_mag)

        #Horizontal illuminance [mlx]
        horizontal_illum = illuminance_horizontal(img)
        print("Horizontal Illuminance:", horizontal_illum)
    
    
    
'''    
    for f in files:
                
        hdr = fits.open(f,uint=False)[0].header

        try:
            utc = datetime.strptime(hdr['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
        except ValueError:
            try:
                utc = datetime.strptime(hdr['DATE-OBS'], '%Y-%m-%dT%H:%M:%S')
            except ValueError:
                print(f"Unrecognized DATE-OBS format in file: {f}")

        
        lat = parse_lat_lon(hdr['SITELAT'])
        lon = -abs(parse_lat_lon(hdr['SITELONG']))
        
       # t = projection.get_local_time_from_utc(lat, lon, utc)
       # date = str(t.date())
       # hour = t.strftime("%H:%M")
       # print(hdr['PARKNAME'])
       # print(hdr['LOCATION'])        
       # print(date, hour)
        
        
       # # Append row to list
       # data_rows.append({
       #     'Parkname': hdr['PARKNAME'],
       #     'Location': hdr['LOCATION'],
       #     'Date': date,
       #     'Time (LMT)': hour,
       #     'Zenith': zenith_mag
       # })

        
        #horizontal illuminance in mlx
        #print(f)
        #horizontal_illum = illuminance_horizontal(image)
        #print("Horizontal Illuminance:", horizontal_illum)
        

# Create DataFrame and save to Excel
#df = pd.DataFrame(data_rows)
#df.to_csv('../zenith_brightness_report_DWH.csv', index=False)
'''    
    
    