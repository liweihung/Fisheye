#-----------------------------------------------------------------------------#
#concurrent_metric_comparison_AI.py
#
#NPS Night Skies Program
#
#Last updated: 2026/08/27
#
#This script computes the median zenith brightness (mag per square arcsec),
#horizontal illuminance (mlx), maximum vertical illuminance (mlx), and
#all-sky Light Pollution Ratio (ALR) of each fisheye image listed in
#concurrent_observations.xlsx, using the same math as metrics.py (including
#the calibrated, non-equidistant theta(r) relationship from theta_r.xlsx),
#and pairs each with the corresponding CCD zenith brightness
#(ZENITH_LUM_MSA), CCD horizontal illuminance (HORIZ_MLX), CCD maximum
#vertical illuminance (MAXVERT_MLX), and CCD ALR (ALR_POS) from the CCD
#database. All columns from concurrent_observations.xlsx are carried
#through, with eight new columns added: fisheye/CCD zenith brightness,
#fisheye/CCD horizontal illuminance, fisheye/CCD maximum vertical
#illuminance, and fisheye/CCD ALR.
#
#Matching metrics.py, zenith brightness and ALR use the terrain mask
#(mask.fit), so that they reflect only the actually-visible sky at this
#site. Horizontal and vertical illuminance use a simple 90-degree horizon
#mask instead, so that they reflect the full hemisphere, including light
#from terrain foreground. Per-pixel solid angle is computed on a Cartesian
#pixel grid as dOmega = dTheta^2 (no azimuthal/sin(theta) polar-grid
#correction), matching metrics.py.
#
#This script lives in Scripts/Others/, one level deeper than the other
#pipeline scripts in Scripts/, so its relative paths climb one extra level.
#
#CCD values are matched using both the CCD Dataset (DNIGHT) and CCD Dset
#columns from concurrent_observations.xlsx, since a single DNIGHT can
#contain multiple dsets/sites in the CCD database. DSET is normalized to
#string on both sides before matching, since the source spreadsheet may
#store it as text.
#
#Camera name per fisheye dataset (e.g. 'Fish5') is looked up from the
#'Camera' column of the fisheye database's Log sheet, which records the
#camera number under the 'Fish<N>' convention used in imagecenter.csv,
#rather than from a hardcoded list in this script.
#
#In the output .xlsx, fisheye-derived value columns are highlighted dark
#orange and CCD-derived value columns are highlighted dark gray, to make it
#easy to visually distinguish the two instruments' measurements.
#
#Note: this script requires the openpyxl package to read and write .xlsx
#files (not otherwise used elsewhere in this pipeline). Install with:
#   pip install openpyxl
#
#Input:
#   (1) ../../Performance/Concurrent_observations/
#       concurrent_observations.xlsx -- output of
#       list_bestpick_observing_times.py
#	(2) ../../Performance/Concurrent_observations/
#	    Data_summary_fisheye_20260827.xlsx -- fisheye database (Log
#	    sheet), for the Camera column
#	(3) ../../Calibration/imagecenter.csv -- fisheye center coordinates and
#	    FOV radius, keyed by camera
#	(4) ../../Calibration/theta_r.xlsx -- empirically measured theta(r)
#	    calibration table
#	(5) Calibrated fisheye images referenced in
#	    concurrent_observations.xlsx (data_cal+'Fisheye File' for each
#	    Fisheye Dataset)
#	(6) data_cal+'mask.fit' -- terrain mask for each fisheye dataset
#	(7) ../../Performance/Concurrent_observations/
#	    Data_summary_CCD_modified.xlsx -- CCD database, for
#	    ZENITH_LUM_MSA, HORIZ_MLX, MAXVERT_MLX, and ALR_POS
#
#Output:
#   (1) ../../Performance/Concurrent_observations/
#       concurrent_metric_comparison.xlsx -- all columns from
#       concurrent_observations.xlsx, plus fisheye/CCD zenith brightness,
#       fisheye/CCD horizontal illuminance, fisheye/CCD maximum vertical
#       illuminance, and fisheye/CCD ALR, with column widths auto-sized to
#       fit their content and color-coded by instrument
#
#History:
#	Li-Wei Hung -- Created
#
#-----------------------------------------------------------------------------#
import numpy as n
import pandas as pd

from astropy.io import fits
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter
from scipy.interpolate import interp1d

#-----------------------------------------------------------------------------#
#                              File locations                                 #
#-----------------------------------------------------------------------------#
#Root folder containing each dataset's processed data, e.g.
#DATA_ROOT+'NIOB_20240830C/'. This script lives in Scripts/Others/, so an
#extra '../' is needed to reach Data_processed/ relative to the Scripts/
#folder. Adjust to match your directory layout.
DATA_ROOT = '../../Data_processed/'

#Calibration folder, for imagecenter.csv and theta_r.xlsx
CALIBRATION = '../../Calibration/'

#Input table produced by list_bestpick_observing_times.py
CONCURRENT_OBS_PATH = '../../Performance/Concurrent_observations/concurrent_observations.xlsx'

#Fisheye database (Log sheet), for the Camera column
DATA_SUMMARY_PATH = '../../Performance/Concurrent_observations/Data_summary_fisheye_20260827.xlsx'

#CCD database location, for ZENITH_LUM_MSA, HORIZ_MLX, MAXVERT_MLX, and
#ALR_POS
NPMAPS_PATH = '../../Performance/Concurrent_observations/Data_summary_CCD_modified.xlsx'

#Output file
OUTPUT_PATH = '../../Performance/Concurrent_observations/concurrent_metric_comparison.xlsx'

#Azimuth sweep interval for maximum vertical illuminance, matching
#metrics.py
VERTICAL_STEP_DEG = 5

#Text color (ARGB hex, no leading '#') for fisheye vs. CCD derived value
#columns in the output .xlsx
FISHEYE_TEXT_COLOR = 'FFCC5500'  #dark orange
CCD_TEXT_COLOR = 'FF595959'	  #dark gray

#Columns to center-align in the output .xlsx
CENTER_ALIGN_COLS = ['CCD Dset', 'Date', 'CCD Mid Obs Time (LMT)',
					 'Fisheye Obs Time (LMT)', 'Time Difference (min)',
					 'F Zenith', 'C Zenith', 'F Horizontal', 'C Horizontal',
					 'F Vertical Max', 'C Vertical Max', 'F ALR', 'C ALR']

#Columns to force to display 2 decimal places, even for whole numbers
DECIMAL_COLS = ['F Zenith', 'C Zenith', 'F Horizontal', 'C Horizontal',
				'F Vertical Max', 'C Vertical Max', 'F ALR', 'C ALR']

#-----------------------------------------------------------------------------#

def zenith(img, r, a=4.75):
	"""
	Calculate the median zenith brightness, identical to metrics.py's
	zenith() function.

	Parameters
	----------
	img : 2D array
		Fisheye image, masked with the terrain mask, calibrated in mag
		per square arcsec.
	r : 2D array
		Pixel radius from the image center, same shape as img.
	a : number, optional
		Radius in pixels; 1 pix = 379 arcsec; aperture diameter should be
		set to about 1 degree, comparable to CCD zenith aperture.

	Returns
	-------
	zenith_mag : float
		Median zenith brightness in mag per square arcsec.
	"""
	zenith_mag = round(n.median(img[n.where(r<a)]), 2)
	return zenith_mag


def load_theta_r_calibration(path):
	"""
	Load the experimentally measured theta(r) calibration table and build
	interpolation functions for theta(r) and dtheta/dr(r), identical to
	metrics.py's load_theta_r_calibration() function.

	Parameters
	----------
	path : str
		Path to the theta_r.xlsx calibration file.

	Returns
	-------
	theta_of_r : callable
		Function mapping pixel radius [pix] to zenith angle [radians].
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
	step_deg = theta_deg[0]

	theta_deg_full = n.concatenate(([0], theta_deg))
	r_full = n.concatenate(([0], r_cumulative))

	theta_of_r = interp1d(r_full, n.deg2rad(theta_deg_full),
						  bounds_error=False,
						  fill_value=(0, n.nan),
						  kind='linear')

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
	Build the per-pixel Theta and dOmega arrays using the calibrated
	theta(r) relationship, identical to metrics.py's
	build_calibrated_geometry() function. Since the image is sampled on a
	Cartesian (not polar) pixel grid, dOmega is simply dTheta^2, with no
	separate azimuthal or sin(theta) weighting.

	Parameters
	----------
	r : 2D array
		Pixel radius from the image center.
	theta_of_r : callable
		Interpolated theta(r) function, from load_theta_r_calibration().
	dtheta_dr_of_r : callable
		Interpolated dtheta/dr(r) function, from load_theta_r_calibration().
	R : float
		Maximum calibrated pixel radius (at theta=90 degrees).

	Returns
	-------
	Theta : 2D array
		Zenith angle [radians] at each pixel. NaN beyond the calibrated FoV.
	dOmega : 2D array
		Solid angle [steradians] subtended by each pixel, dOmega = dTheta^2.
		NaN beyond the calibrated FoV.
	"""
	Theta = theta_of_r(r)
	dTheta = dtheta_dr_of_r(r)

	dOmega = dTheta**2

	beyond_fov = r > R
	Theta[beyond_fov] = n.nan
	dOmega[beyond_fov] = n.nan

	return Theta, dOmega


def illuminance_horizontal(img, Theta, dOmega):
	"""
	Calculate the horizontal illuminance, identical to metrics.py's
	illuminance_horizontal() function.

	Parameters
	----------
	img : 2D array
		Fisheye image, masked with the 90-degree horizon mask, calibrated
		in mag per square arcsec.
	Theta : 2D array
		Zenith angle [radians] at each pixel, same shape as img.
	dOmega : 2D array
		Solid angle [steradians] subtended by each pixel, same shape as img.

	Returns
	-------
	E_h : float
		Horizontal illuminance in mlx.
	"""
	L = 108.48*n.exp(20.7233-0.92104*img) # [ucd m-2], Duriscoe 2016 conversion
	dE = L*n.cos(Theta)*dOmega
	E_h = n.nansum(dE)/1000 #Horizontal illuminance [mlx]
	return E_h


def precompute_vertical_geometry(Phi, step_deg=5):
	"""
	Precompute the azimuth sweep grid and cos(incidence) weighting array
	used by illuminance_vertical_max(), identical to metrics.py's
	precompute_vertical_geometry() function.

	Parameters
	----------
	Phi : 2D array
		Azimuth [radians] at each pixel, compass convention (0 = North,
		clockwise through E/S/W).
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
	facing_azimuths = n.deg2rad(azimuths_deg)

	cos_incidence = n.cos(Phi[n.newaxis,:,:] - facing_azimuths[:,n.newaxis,n.newaxis])
	cos_incidence = n.clip(cos_incidence, 0, None)

	return azimuths_deg, cos_incidence


def illuminance_vertical_max(img, Theta, dOmega, azimuths_deg, cos_incidence):
	"""
	Calculate the maximum vertical illuminance over all swept facing
	directions, identical to metrics.py's illuminance_vertical_max()
	function. This matches the MAXVERT_MLX metric in the CCD database.

	Parameters
	----------
	img : 2D array
		Fisheye image, masked with the 90-degree horizon mask, calibrated
		in mag per square arcsec.
	Theta : 2D array
		Zenith angle [radians] at each pixel, same shape as img.
	dOmega : 2D array
		Solid angle [steradians] subtended by each pixel, same shape as img.
	azimuths_deg : 1D array
		Facing azimuths swept, in degrees, from precompute_vertical_geometry().
	cos_incidence : 3D array
		Precomputed cos(incidence) weighting array, from
		precompute_vertical_geometry().

	Returns
	-------
	E_v_max : float
		Maximum vertical illuminance in mlx, over all sweep directions.
	best_azimuth_deg : float
		Facing azimuth (degrees, compass convention, 0 = North) at which
		the maximum vertical illuminance occurs.
	"""
	L = 108.48*n.exp(20.7233-0.92104*img) # [ucd m-2], Duriscoe 2016 conversion

	base = L*n.sin(Theta)*dOmega
	dE = base[n.newaxis,:,:] * cos_incidence
	E_v_sweep = n.nansum(dE, axis=(1,2)) / 1000

	E_v_max = round(float(E_v_sweep.max()), 2)
	best_azimuth_deg = azimuths_deg[n.argmax(E_v_sweep)]

	return E_v_max, best_azimuth_deg


def ALR(img, dOmega, natural_reference=250):
	"""
	Calculate the All-sky Light Pollution Ratio (ALR), identical to
	metrics.py's ALR() function. ALR is the total skyglow brightness
	divided by the natural dark sky reference value. Each pixel's
	luminance is weighted by its solid angle dOmega, so that pixels
	representing more sky area contribute proportionally more to the
	average -- necessary because this fisheye image is not in an
	equal-area projection.

	Parameters
	----------
	img : 2D array
		Fisheye image, calibrated in mag per square arcsec, already
		multiplied by the terrain mask, so that terrain/obstructed pixels
		are NaN and excluded from the average.
	dOmega : 2D array
		Solid angle [steradians] subtended by each pixel, same shape as
		img, from build_calibrated_geometry().
	natural_reference : number, optional
		Natural reference luminance value [ucd/m^2] used as the denominator
		of the ratio. Defaults to 250 ucd/m^2, the median natural all-sky
		brightness (sky+stars) from Duriscoe (2016).

	Returns
	-------
	alr : float
		All-sky Light Pollution Ratio (dimensionless).
	"""
	L = 108.48*n.exp(20.7233-0.92104*img) #[ucd m-2]

	#solid-angle-weighted average all-sky brightness, artificial light only
	b_average = n.nansum(L*dOmega)/n.nansum(dOmega) - natural_reference #[ucd m-2]

	#light pollution ratio
	alr = round(float(b_average/natural_reference), 2)

	return alr


def get_fisheye_metrics(dataset, fname, camera, C, theta_of_r, dtheta_dr_of_r, R_cal):
	"""
	Compute the zenith brightness, horizontal illuminance, maximum
	vertical illuminance, and ALR for a single fisheye image, using the
	same geometry and math as metrics.py. Zenith brightness and ALR use
	the dataset's terrain mask (mask.fit), so that they reflect only the
	actually-visible sky at this site. Horizontal and vertical illuminance
	use a simple 90-degree horizon mask instead, so that they reflect the
	full hemisphere. The calibrated theta(r) relationship (theta_r.xlsx)
	is used instead of the equidistant assumption.

	Parameters
	----------
	dataset : str
		Fisheye dataset folder name, e.g. 'NIOB_20240830C'.
	fname : str
		File name of the image within the dataset's data_cal folder.
	camera : str
		Camera name matching an entry in imagecenter.csv (e.g. 'Fish5').
	C : DataFrame
		imagecenter.csv contents, indexed by camera name.
	theta_of_r : callable
		Interpolated theta(r) function, from load_theta_r_calibration().
	dtheta_dr_of_r : callable
		Interpolated dtheta/dr(r) function, from load_theta_r_calibration().
	R_cal : float
		Maximum calibrated pixel radius (at theta=90 degrees).

	Returns
	-------
	zenith_mag : float or None
		Median zenith brightness in mag per square arcsec, or None if the
		image could not be read.
	horiz_illum : float or None
		Horizontal illuminance in mlx, or None if the image could not be
		read.
	vert_illum_max : float or None
		Maximum vertical illuminance in mlx, or None if the image could
		not be read.
	alr : float or None
		All-sky Light Pollution Ratio (dimensionless), or None if the
		image could not be read.
	"""
	data_cal = DATA_ROOT+dataset+'/'
	fpath = data_cal+fname

	try:
		with fits.open(fpath, uint=False, memmap=False) as hdul:
			raw_img = hdul[0].data.astype(float, copy=True)
	except Exception as e:
		print('  Could not read %s: %s' %(fpath, e))
		return None, None, None, None

	xc = C['Xcenter'][camera]
	yc = C['Ycenter'][camera]
	R = C['Radius'][camera]

	ny, nx = raw_img.shape
	x, y = n.meshgrid(n.arange(nx), n.arange(ny))
	r = n.sqrt((x-xc)**2 + (y-yc)**2)
	Phi = -n.arctan2(y-yc, x-xc) + n.pi/2

	#terrain mask, for zenith and ALR
	try:
		with fits.open(data_cal+'mask.fit', uint=False, memmap=False) as hdul:
			terrain_mask = hdul[0].data.astype(float, copy=True)
	except Exception as e:
		print('  Could not read mask.fit for %s: %s' %(dataset, e))
		terrain_mask = n.ones_like(raw_img)

	img_terrain = raw_img * terrain_mask

	#simple 90-degree horizon mask, for horizontal/vertical illuminance
	horizon_mask = n.ones_like(raw_img)
	horizon_mask[r > R] = n.nan
	img_horizon = raw_img * horizon_mask

	#zenith brightness -- terrain mask
	zenith_mag = zenith(img_terrain, r)

	#calibrated horizontal/vertical illuminance geometry, identical to
	#metrics.py
	Theta, dOmega = build_calibrated_geometry(r, theta_of_r, dtheta_dr_of_r, R_cal)
	azimuths_deg, cos_incidence = precompute_vertical_geometry(Phi, step_deg=VERTICAL_STEP_DEG)

	#horizontal illuminance -- horizon mask
	horiz_illum = illuminance_horizontal(img_horizon, Theta, dOmega)
	if horiz_illum is not None:
		horiz_illum = round(float(horiz_illum), 2)

	#max vertical illuminance -- horizon mask
	vert_illum_max, _ = illuminance_vertical_max(
		img_horizon, Theta, dOmega, azimuths_deg, cos_incidence)

	#ALR -- terrain mask, solid-angle weighted
	alr = ALR(img_terrain, dOmega)

	return zenith_mag, horiz_illum, vert_illum_max, alr


def save_with_autofit_columns(df, outpath, fisheye_cols=None, ccd_cols=None):
	"""
	Write a DataFrame to an .xlsx file with each column's width set to fit
	its content (based on the longest value in the column, not the
	header), using openpyxl, with a minimum width equal to the width of
	the 'Date' column. Columns listed in fisheye_cols have their text
	colored dark orange; columns listed in ccd_cols have their text
	colored dark gray.

	Parameters
	----------
	df : DataFrame
		Table to write.
	outpath : str
		Destination .xlsx file path.
	fisheye_cols : list of str, optional
		Column names to color as fisheye-derived values (dark orange text).
	ccd_cols : list of str, optional
		Column names to color as CCD-derived values (dark gray text).
	"""
	if fisheye_cols is None:
		fisheye_cols = []
	if ccd_cols is None:
		ccd_cols = []

	fisheye_font = Font(color=FISHEYE_TEXT_COLOR)
	ccd_font = Font(color=CCD_TEXT_COLOR)
	center_align = Alignment(horizontal='center', vertical='center')

	#minimum column width, taken from the content length of the 'Date'
	#column, so every column is at least as wide as Date
	if 'Date' in df.columns:
		date_values = df['Date'].tolist()
		min_width = max([len('Date')] + [len(str(v)) for v in date_values]) + 2
	else:
		min_width = 10

	with pd.ExcelWriter(outpath, engine='openpyxl') as writer:
		df.to_excel(writer, index=False, sheet_name='Sheet1')
		ws = writer.sheets['Sheet1']

		for i, col in enumerate(df.columns, start=1):
			col_letter = get_column_letter(i)

			values = df[col].tolist()
			max_len = max([len(str(v)) for v in values]) if values else 0
			ws.column_dimensions[col_letter].width = max(max_len+2, min_width)

			if col in fisheye_cols:
				font = fisheye_font
			elif col in ccd_cols:
				font = ccd_font
			else:
				font = None

			if font is not None:
				for row in range(1, len(values)+2):
					ws[f'{col_letter}{row}'].font = font

			if col in CENTER_ALIGN_COLS:
				for row in range(1, len(values)+2):
					ws[f'{col_letter}{row}'].alignment = center_align

			if col in DECIMAL_COLS:
				for row in range(2, len(values)+2):
					ws[f'{col_letter}{row}'].number_format = '0.00'


def main():
	"""
	Computes the fisheye zenith brightness, horizontal illuminance,
	maximum vertical illuminance, and ALR for every row in
	concurrent_observations.xlsx, looks up the matching CCD metrics from
	the CCD database, and writes the combined table. See the script
	description for detail.
	"""

	#--------------------------------------------------------------------------#
	#							  Load input tables							   #
	#--------------------------------------------------------------------------#
	print('Loading concurrent observations from %s' %CONCURRENT_OBS_PATH)
	obs = pd.read_excel(CONCURRENT_OBS_PATH)

	#Date may be read in as a full Timestamp (with a 00:00:00 time
	#component); keep only the date portion for display
	if 'Date' in obs.columns:
		obs['Date'] = pd.to_datetime(obs['Date']).dt.date

	print('Loading fisheye database from %s' %DATA_SUMMARY_PATH)
	log = pd.read_excel(DATA_SUMMARY_PATH, sheet_name='Log')

	#build a Dataset -> camera name lookup (e.g. 'Fish5') from the Log
	#sheet's Camera column, which records the camera number under the
	#'Fish<N>' convention used in imagecenter.csv
	cameras = {row['Dataset']: 'Fish'+str(row['Camera']).strip()
			   for _, row in log.iterrows()}

	print('Loading CCD database from %s' %NPMAPS_PATH)
	npmaps = pd.read_excel(NPMAPS_PATH)

	#DSET is sometimes read as text (e.g. "1" instead of 1) depending on how
	#the source spreadsheet stores it; normalize both sides to string so
	#the concurrent_observations.xlsx CCD Dset column still matches
	npmaps['DSET'] = npmaps['DSET'].astype(str).str.strip()
	obs['CCD Dset'] = obs['CCD Dset'].astype(str).str.strip()

	print('Loading fisheye camera geometry from %simagecenter.csv' %CALIBRATION)
	C = pd.read_csv(CALIBRATION+'imagecenter.csv', index_col=0)

	print('Loading theta(r) calibration from %stheta_r.xlsx' %CALIBRATION)
	theta_of_r, dtheta_dr_of_r, R_cal = load_theta_r_calibration(CALIBRATION+'theta_r.xlsx')

	#--------------------------------------------------------------------------#
	#		  Compute fisheye metrics and look up matching CCD metrics	   #
	#--------------------------------------------------------------------------#
	fisheye_zeniths = []
	ccd_zeniths = []
	fisheye_horiz_illums = []
	ccd_horiz_illums = []
	fisheye_vert_illums = []
	ccd_vert_illums = []
	fisheye_alrs = []
	ccd_alrs = []

	for _, row in obs.iterrows():

		dataset = row['Fisheye Dataset']
		fname = row['Fisheye File']
		ccd_dnight = row['CCD Dataset']
		ccd_dset = row['CCD Dset']

		print('Processing %s / %s ...' %(dataset, fname))

		#fisheye zenith brightness, horizontal illuminance, max vertical
		#illuminance, and ALR
		camera = cameras.get(dataset)
		if camera is None:
			print('  No camera specified for %s, skipping.' %dataset)
			fisheye_zeniths.append(None)
			fisheye_horiz_illums.append(None)
			fisheye_vert_illums.append(None)
			fisheye_alrs.append(None)
		else:
			zenith_mag, horiz_illum, vert_illum_max, alr = get_fisheye_metrics(
				dataset, fname, camera, C, theta_of_r, dtheta_dr_of_r, R_cal)
			fisheye_zeniths.append(zenith_mag)
			fisheye_horiz_illums.append(horiz_illum)
			fisheye_vert_illums.append(vert_illum_max)
			fisheye_alrs.append(alr)

		#CCD zenith brightness, horizontal illuminance, max vertical
		#illuminance, and ALR -- match on both DNIGHT and DSET, since a
		#single DNIGHT can contain multiple dsets/sites
		ccd_match = npmaps[(npmaps['DNIGHT'] == ccd_dnight) &
							(npmaps['DSET'] == ccd_dset)]
		if ccd_match.empty:
			print('  No CCD match found for %s dset %s' %(ccd_dnight, ccd_dset))
			ccd_zeniths.append(None)
			ccd_horiz_illums.append(None)
			ccd_vert_illums.append(None)
			ccd_alrs.append(None)
		else:
			rec = ccd_match.iloc[0]
			ccd_zeniths.append(rec.get('ZENITH_LUM_MSA', None))
			ccd_horiz_illums.append(rec.get('HORIZ_MLX', None))
			ccd_vert_illums.append(rec.get('MAXVERT_MLX', None))
			ccd_alrs.append(rec.get('ALR_POS', None))

	obs['F Zenith'] = fisheye_zeniths
	obs['C Zenith'] = ccd_zeniths
	obs['F Horizontal'] = fisheye_horiz_illums
	obs['C Horizontal'] = ccd_horiz_illums
	obs['F Vertical Max'] = fisheye_vert_illums
	obs['C Vertical Max'] = ccd_vert_illums
	obs['F ALR'] = fisheye_alrs
	obs['C ALR'] = ccd_alrs

	#ensure fisheye values are rounded to 2 decimal places in the output,
	#even if float precision drifted during DataFrame assignment
	obs['F Horizontal'] = obs['F Horizontal'].round(2)
	obs['F Vertical Max'] = obs['F Vertical Max'].round(2)
	obs['F ALR'] = obs['F ALR'].round(2)

	#round the time difference to the nearest whole minute
	if 'Time Difference (min)' in obs.columns:
		obs['Time Difference (min)'] = obs['Time Difference (min)'].round(0).astype('Int64')

	#--------------------------------------------------------------------------#
	#								Save results							   #
	#--------------------------------------------------------------------------#
	save_with_autofit_columns(
		obs, OUTPUT_PATH,
		fisheye_cols=['F Zenith', 'F Horizontal', 'F Vertical Max', 'F ALR'],
		ccd_cols=['C Zenith', 'C Horizontal', 'C Vertical Max', 'C ALR'])
	print('\nSaved metric comparison table to', OUTPUT_PATH)
	print(obs)


if __name__ == '__main__':
	main()