#------------------------------------------------------------------------------#
#metrics.py
#
#NPS Night Skies Program
#
#Given the calibrated fisheye images, the script extracts various night sky 
#brightness metrics. 
#
#
#Input: 
#   (1) theta_r.xlsx -- empirically measured theta(r) calibration table
#	(2) 
#
#Output:
#   (1) p.data_cal+'metrics.xlsx' -- Date, Time (LMT), File name, Zenith,
#       Horizontal Illuminance, Max Vertical Illuminance, and its azimuth,
#       for every processed image in the dataset
#	(2) p.data_cal+'vertical_illuminance.xlsx' -- full azimuth sweep of
#	    vertical illuminance for the dataset's reference image only
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
from openpyxl.styles import Alignment
from openpyxl.utils import get_column_letter
from scipy.interpolate import interp1d

# Local Source
#import process_input as p     
import projection

#------------------------------------------------------------------------------#

def get_obs_date_time_lmt(hdr, file):
    """
	Extract the observing date and local mean time from a FITS header's
	DATE-OBS (UTC) and site coordinates.

	Parameters
	----------
	hdr : astropy FITS header
		Header of the image being processed, expected to contain DATE-OBS,
		SITELAT, and SITELONG.
	file : str
		File name of the image, used only for the warning message if
		DATE-OBS cannot be parsed.

	Returns
	-------
	date : str or None
		Observing date (YYYY-MM-DD), or None if DATE-OBS could not be
		parsed or the local time could not be determined.
	hour : str or None
		Observing time (HH:MM) in local mean time, or None under the same
		conditions as date.
	"""
    try:
        utc = datetime.strptime(hdr['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f')
    except ValueError:
        try:
            utc = datetime.strptime(hdr['DATE-OBS'], '%Y-%m-%dT%H:%M:%S')
        except ValueError:
            print(f"Unrecognized DATE-OBS format in file: {file}")
            utc = None

    date, hour = None, None
    if utc is not None:
        lat = parse_lat_lon(hdr['SITELAT'])
        lon = -abs(parse_lat_lon(hdr['SITELONG']))
        t = projection.get_local_time_from_utc(lat, lon, utc)
        if t is not None:
            date = str(t.date())
            hour = t.strftime("%H:%M")

    return date, hour

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


def load_theta_r_calibration(path):
	"""
	Load the experimentally measured theta(r) calibration table and build
	interpolation functions for theta(r) and dtheta/dr(r). This replaces
	the equidistant fisheye assumption (theta = pi/2 * r/R, constant
	dtheta/dr) with the actual, non-equidistant relationship measured for
	this lens.

	theta_r.xlsx contains:
		Zenith_Angle (deg) | Delta_Pixel | Cumulative_Pixel
	sampled every 3 degrees from 3 to 90. Cumulative_Pixel is r(theta) at
	that zenith angle; Delta_Pixel is the pixel distance covered by each
	3-degree step, i.e. Delta_r for that band.

	Parameters
	----------
	path : str, optional
		Path to the theta_r.xlsx calibration file.

	Returns
	-------
	theta_of_r : callable
		Function mapping pixel radius [pix] to zenith angle [radians].
		Valid for r in [0, R], where R is the maximum calibrated radius
		(corresponding to theta=90 degrees).
	dtheta_dr_of_r : callable
		Function mapping pixel radius [pix] to the local derivative
		d(theta)/dr [radians per pixel] at that radius.
	R : float
		Maximum calibrated pixel radius (at theta=90 degrees).
	"""
	cal = pd.read_excel(path)

	theta_deg = cal['Zenith_Angle (deg)'].values
	r_cumulative = cal['Cumulative_Pixel'].values
	delta_pixel = cal['Delta_Pixel'].values
	step_deg = theta_deg[0] #assumes a uniform step, e.g. 3 degrees

	#prepend the origin (r=0 at theta=0), since the table starts at 3 deg
	theta_deg_full = n.concatenate(([0], theta_deg))
	r_full = n.concatenate(([0], r_cumulative))

	#theta(r): interpolate zenith angle [radians] as a function of pixel
	#radius
	theta_of_r = interp1d(r_full, n.deg2rad(theta_deg_full), 
						  bounds_error=False,
						  fill_value=(0, n.nan),
						  kind='linear')
						  

	#dtheta/dr(r): local angular scale [radians per pixel] within each
	#3-degree band, computed directly from the table (step_deg / Delta_Pixel
	#for each band), then interpolated as a function of r, associated with
	#the midpoint of each band's r-range
	dtheta_dr_per_band = n.deg2rad(step_deg) / delta_pixel
	r_prev = n.concatenate(([0], r_cumulative[:-1]))
	r_band_mid = (r_prev + r_cumulative) / 2

	dtheta_dr_of_r = interp1d(r_band_mid, dtheta_dr_per_band,
							  kind='linear', bounds_error=False,
							  fill_value=(dtheta_dr_per_band[0],
										  dtheta_dr_per_band[-1]))

	R = r_cumulative[-1]

	return theta_of_r, dtheta_dr_of_r, R


def build_calibrated_geometry(r, theta_of_r, dtheta_dr_of_r, R):
	"""
	Build the per-pixel r, Theta, dTheta, dPhi, and dOmega arrays using the
	calibrated theta(r) relationship, replacing the equidistant assumption.

	Parameters
	----------
	r : 2D array
		Pixel radius from the image center.
	theta_of_r : callable
		Interpolated theta(r) function, from load_theta_r_calibration().
	dtheta_dr_of_r : callable
		Interpolated dtheta/dr(r) function, from load_theta_r_calibration().
	R : float
		Maximum calibrated pixel radius (at theta=90 degrees). Pixels
		beyond this radius are outside the calibrated field of view and are
		masked to NaN.

	Returns
	-------
	Theta : 2D array
		Zenith angle [radians] at each pixel, from the calibrated theta(r).
		NaN beyond the calibrated field of view.
	dTheta : 2D array
		Local angular step [radians per pixel] at each pixel's radius.
	dPhi : 2D array
		Azimuthal step [radians] per pixel, unaffected by the radial
		projection law: dPhi = 1/r.
	dOmega : 2D array
		Solid angle [steradians] subtended by each pixel. NaN beyond the
		calibrated field of view.
	"""

	Theta = theta_of_r(r)
	dTheta = dtheta_dr_of_r(r)

	with n.errstate(divide='ignore', invalid='ignore'):
		dPhi = 1/r

	dOmega = dTheta * dPhi * n.sin(Theta)

	#special case for the zenith pixel
	dOmega[r==0] = dtheta_dr_of_r(0)**2

	#mask out pixels beyond the calibrated field of view, matching the
	#equidistant geometry's behavior of excluding pixels beyond the FoV
	beyond_fov = r > R
	Theta[beyond_fov] = n.nan
	dOmega[beyond_fov] = n.nan

	return Theta, dTheta, dPhi, dOmega


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
    #mag to [ucd m-2] conversion Duriscoe 2016 equation 1
    L = 108.48*n.exp(20.7233-0.92104*img) # [ucd m-2], Dan's conversion formula
    dE = L*n.cos(Theta)*dOmega 
    E_h = round(n.nansum(dE)/1000, 2) #Horizontal illuminance [mlx]

    return E_h


def precompute_vertical_geometry(Phi, step_deg=5):
    """
	Precompute the azimuth sweep grid and cos(incidence) weighting array
	used by illuminance_vertical_max(), so they can be computed once per
	dataset (since Phi is fixed by camera geometry and does not change
	between images) rather than being rebuilt for every image processed.

	Parameters
	----------
	Phi : 2D array
		Azimuth [radians] at each pixel, compass convention (0 = North,
		clockwise through E/S/W), same shape as the images to be processed.
	step_deg : number, optional
		Azimuth sweep interval in degrees. Defaults to 5.

	Returns
	-------
	azimuths_deg : 1D array
		All facing azimuths swept, in degrees, from 0 up to (but not
		including) 360.
	cos_incidence : 3D array
		cos(Phi-facing_azimuth) for each swept azimuth, clipped to exclude
		light arriving from behind the surface. Shape: (n_az, ny, nx).

	"""

    azimuths_deg = n.arange(0, 360, step_deg)
    facing_azimuths = n.deg2rad(azimuths_deg)	  #shape: (n_az,)

    #broadcast Phi (ny, nx) against facing_azimuths (n_az,) by adding a new
    #axis: result shape (n_az, ny, nx)
    cos_incidence = n.cos(Phi[n.newaxis,:,:] - facing_azimuths[:,n.newaxis,n.newaxis])
    cos_incidence = n.clip(cos_incidence, 0, None) #exclude light from behind

    return azimuths_deg, cos_incidence


def illuminance_vertical_max(img, azimuths_deg, cos_incidence):
    """
	Calculate the vertical illuminance at every swept facing direction, and
	identify the maximum, using precomputed sweep geometry from
	precompute_vertical_geometry(). This matches the MAXVERT_MLX metric in
	the CCD database (NPMaps.csv).

	Parameters
	----------
	img : 2D array
		Fisheye image
	azimuths_deg : 1D array
		Facing azimuths swept, in degrees, from precompute_vertical_geometry().
	cos_incidence : 3D array
		Precomputed cos(incidence) weighting array, from
		precompute_vertical_geometry(). Shape: (n_az, ny, nx).

	Returns
	-------
	E_v_sweep : 1D array
		Vertical illuminance in mlx at each corresponding azimuth in
		azimuths_deg.
	E_v_max : float
		Maximum vertical illuminance in mlx, over all sweep directions.
	best_azimuth_deg : float
		Facing azimuth (degrees, compass convention, 0 = North) at which
		the maximum vertical illuminance occurs.

	"""
    #mag to [ucd m-2] conversion Duriscoe 2016 equation 1
    L = 108.48*n.exp(20.7233-0.92104*img) #[ucd m-2]

    #per-pixel quantities that don't depend on facing azimuth
    base = L*n.sin(Theta)*dOmega	 #shape: (ny, nx)

    dE = base[n.newaxis,:,:] * cos_incidence		 #shape: (n_az, ny, nx)
    E_v_sweep = n.nansum(dE, axis=(1,2)) / 1000	 #shape: (n_az,); [mlx]

    E_v_max = round(E_v_sweep.max(),2)
    best_azimuth_deg = azimuths_deg[n.argmax(E_v_sweep)]

    return E_v_sweep, E_v_max, best_azimuth_deg


def ALR(img, natural_reference=250):
    """
	Calculate the All-sky Light Pollution Ratio (ALR), as defined in 
	Duriscoe (2016). ALR is the total skyglow brightness divided by the 
	natural dark sky set at 250 ucd/m^2. 

	Parameters
	----------
	img : 2D array
		Fisheye image, calibrated in mag per square arcsec, already
		multiplied by the terrain mask, so that	terrain/obstructed pixels 
		are NaN and excluded from the average.
	natural_reference : number, optional
		Natural reference luminance value [ucd/m^2] used as the denominator
		of the ratio. Defaults to 250 ucd/m^2, the median natural all-sky
		brightness (sky+stars) from Duriscoe (2016).

	Returns
	-------
	alr : float
		All-sky Light Pollution Ratio (dimensionless).

	"""
    #mag to [ucd m-2] conversion Duriscoe 2016 equation 1
    L = 108.48*n.exp(20.7233-0.92104*img) #[ucd m-2]

    #average all-sky brightness, artificial light only
    b_average = n.nanmean(L) - natural_reference #[ucd m-2]

    #light pollution ratio
    alr = round(b_average/natural_reference, 2)

    return alr


def save_with_autofit_columns(df, outpath, min_width=10, decimal_cols=None):
    """
	Write a DataFrame to an .xlsx file with each column's width set to fit
	its content (based on the longest value in the column, not the header),
	using openpyxl. All columns are center-justified. Columns listed in
	decimal_cols are formatted to always display 2 decimal places, even for
	whole numbers.

	Parameters
	----------
	df : DataFrame
		Table to write.
	outpath : str
		Destination .xlsx file path.
	min_width : number, optional
		Minimum column width, regardless of content length. Defaults to 10
		(comparable to a 'YYYY-MM-DD' date column). 
	decimal_cols : list of str, optional
		Column names to format with a fixed 2 decimal places.
	"""
    if decimal_cols is None:
        decimal_cols = []

    with pd.ExcelWriter(outpath, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Sheet1')
        ws = writer.sheets['Sheet1']

        center = Alignment(horizontal='center', vertical='center')

        for i, col in enumerate(df.columns, start=1):
            col_letter = get_column_letter(i)

            values = df[col].tolist()
            max_len = max([len(str(v)) for v in values]) if values else 0
            ws.column_dimensions[col_letter].width = max(max_len+2, min_width)

            #center-justify every cell in this column, including the header
            for row in range(1, len(values)+2):
                cell = ws[f'{col_letter}{row}']
                cell.alignment = center

            #force 2 decimal places for designated numeric columns
            if col in decimal_cols:
                for row in range(2, len(values)+2):
                    ws[f'{col_letter}{row}'].number_format = '0.00'

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
#		   glob('../Data_processed/ROMO_20241003B_testcopy/process_input.py')
#    	   glob'../Data_processed/GUIS_20231210B/process_input.py') #+ \
#		   glob('../Data_processed/EGLI*/process_input.py') + \
#          glob('../Data_processed/BOSE*/process_input.py')

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

    #Image geometry in polar coordinates
    image_shape = mask.shape
    x, y = n.meshgrid(n.arange(image_shape[1]),n.arange(image_shape[0]))
    r = n.sqrt((x-xc)**2 + (y-yc)**2)
    Phi = -n.arctan2(y-yc, x-xc) + n.pi/2

    #Illuminance parameters, using the calibrated, non-equidistant
    #theta(r) relationship (theta_r.xlsx) rather than the equidistant
    #assumption theta = pi/2 * r/R
    theta_r_file = p.calibration+'theta_r.xlsx'
    theta_of_r, dtheta_dr_of_r, R_cal = load_theta_r_calibration(path=theta_r_file)
    Theta, dTheta, dPhi, dOmega = build_calibrated_geometry(
        r, theta_of_r, dtheta_dr_of_r, R_cal)
    azimuths_deg, cos_incidence = precompute_vertical_geometry(Phi, step_deg=5)

    #List of images to process        
    images = glob(p.data_cal+f'Light*.fit') + glob(p.data_cal+f'img*sky*.fit')
	
    reference_files = [os.path.normpath(f) for f in glob(p.data_cal+p.reference)]

    #reset per-dataset row collection for the metrics.xlsx summary table
    data_rows = []
    
    for file in images:

        print(f"Processing {file} ...")
        #Read in the image and apply the terrian mask
        with fits.open(file, uint=False, memmap=False) as hdul:
            hdr = hdul[0].header
            img = hdul[0].data.astype(float, copy=True) * mask

        #Observing date and time (LMT), from the FITS header
        date, hour = get_obs_date_time_lmt(hdr, file)
        
        #Zenith brightness [mag per squared arcsec]
        zenith_mag = zenith(img)
        print("Zenith:", zenith_mag)

        #Horizontal illuminance [mlx]
        horizontal_illum = illuminance_horizontal(img)
        print("Horizontal Illuminance:", horizontal_illum)

        #Vertical illuminance [mlx] and best facing azimuth [deg]
        E_v_sweep, E_v_max, best_azimuth_deg = illuminance_vertical_max(
        img, azimuths_deg, cos_incidence)
        print(f"Max vertical illuminance is at azimuth {best_azimuth_deg}°: {E_v_max:.2f} mlx")
		
        #All-sky average brightness Light Pollution Ratio (Duriscoe 2016)
        alr = ALR(img)
        print("ALR:", alr)

        #Append this image's metrics to the summary table
        data_rows.append({
            'Date': date,
            'Time': hour, #Local Mean Time
            'File': os.path.basename(file),
            'Zenith': zenith_mag, #[mag/arcsec2]
            'Horizontal': horizontal_illum, #Horizontal Illuminance [mlx]
            'Vertical': E_v_max, #Max Vertical Illuminance [mlx]
            'Azimuth': best_azimuth_deg, #Azimuth of Max Vertical Illuminance (degree)
            'ALR': alr, #All-sky Light Pollution Ratio (dimensionless)
        })

        #Save the full azimuth sweep of vertical illuminance for the
        #dataset's reference image only
        if os.path.normpath(file) in reference_files:
            df_vert = pd.DataFrame({
                'Azimuth': azimuths_deg,
                'Vertical Illuminance (mlx)': n.round(E_v_sweep, 2),
            })
            save_with_autofit_columns(
                df_vert, p.data_cal+'vertical_illuminance.xlsx',
                decimal_cols=['Vertical Illuminance (mlx)'])
            print('Saved vertical illuminance sweep to',
                  p.data_cal+'vertical_illuminance.xlsx')

    #Save the per-image metrics summary table for this dataset
    df_metrics = pd.DataFrame(data_rows)
    save_with_autofit_columns(
        df_metrics, p.data_cal+'metrics.xlsx',
        decimal_cols=['Zenith', 'Horizontal', 'Vertical', 'ALR'])
    print('Saved metrics summary to', p.data_cal+'metrics.xlsx')
			
    
    