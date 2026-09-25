#-----------------------------------------------------------------------------#
#concurrent_metric_comparison_AI.py
#
#NPS Night Skies Program
#
#Last updated: 2026/09/24
#
#This script computes the median zenith brightness, horizontal illuminance,
#maximum vertical illuminance, all-sky Light Pollution Ratio (ALR), Mean
#all-sky brightness, sky brightness percentiles (P1/P50/P99), percentage of
#naked-eye-visible stars (Star), and a synthetic Sky Quality Meter reading
#(SQM_syn) for each fisheye image listed in concurrent_observations.xlsx,
#using the same math as metrics.py (including the calibrated,
#non-equidistant theta(r) relationship from theta_r.xlsx), and pairs each
#with the corresponding CCD metric from the CCD database:
#   fisheye          CCD column
#   ---------------  ----------------
#   Zenith           ZENITH_LUM_MSA
#   Horizontal       HORIZ_MLX
#   Vertical Max     MAXVERT_MLX
#   ALR              ALR_POS
#   Mean             MEANLUM_ART_MAGS
#   P1               P01_ALL_MAGS
#   P50              P50_ALL_MAGS
#   P99              P99_ALL_MAGS
#   Star             VISSTARS_PCT
#   SQM_syn          SYN_SQM
#
#All columns from concurrent_observations.xlsx are carried through, with a
#fisheye/CCD pair of columns added for each metric above. If the CCD
#database also has a real (measured, not synthetic) SQM column, it is
#appended as an extra CCD-only column at the end, since there is no
#fisheye equivalent to pair it with.
#
#Matching metrics.py: Zenith, ALR, Mean, and the percentiles use the
#terrain mask (mask.fit), so that they reflect only the actually-visible
#sky at this site. Horizontal and vertical illuminance use a simple
#90-degree horizon mask instead, so that they reflect the full hemisphere,
#including light from terrain foreground. Per-pixel solid angle is
#dOmega = sin(Theta)*dTheta/r (with the r=0 limit dTheta^2), matching
#metrics.py's build_calibrated_geometry(). Star is predicted from ALR via
#star_visibility_model_AI.predict_visibility(). SQM_syn is computed via
#metrics.py's synthetic_sqm(), using the SQM's angular response
#(sqm_model_AI.sqm_response()) and the same SQM-V offset (0.07).
#
#This script lives in Scripts/Others/, one level deeper than the other
#pipeline scripts in Scripts/, so its relative paths climb one extra level.
#star_visibility_model_AI.py and sqm_model_AI.py live in Scripts/, so '..'
#is added to sys.path before importing them.
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
#	    Data_summary_CCD_modified.xlsx -- CCD database, for the CCD
#	    columns listed above, plus SQM if present
#
#Output:
#   (1) ../../Performance/Concurrent_observations/
#       concurrent_metric_comparison.xlsx -- all columns from
#       concurrent_observations.xlsx, plus fisheye/CCD pairs for every
#       metric listed above (and a CCD-only SQM column if available), with
#       column widths auto-sized to fit their content and color-coded by
#       instrument
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

import sys
sys.path.append('..')
from star_visibility_model_AI import predict_visibility
from sqm_model_AI import sqm_response

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

#CCD database location
NPMAPS_PATH = '../../Performance/Concurrent_observations/Data_summary_CCD_modified.xlsx'

#Output file
OUTPUT_PATH = '../../Performance/Concurrent_observations/concurrent_metric_comparison.xlsx'

#Azimuth sweep interval for maximum vertical illuminance, matching
#metrics.py
VERTICAL_STEP_DEG = 5

#Natural reference luminance [ucd/m^2], for ALR, matching metrics.py
NATURAL_REFERENCE = 250

#SQM-V conversion factor [mag/arcsec^2], matching metrics.py's
#synthetic_sqm() default
SQM_V_OFFSET = 0.07

#Fisheye metric name -> corresponding CCD column name in the CCD database
CCD_COLUMNS = {
	'Zenith':       'ZENITH_LUM_MSA',
	'Horizontal':   'HORIZ_MLX',
	'Vertical Max': 'MAXVERT_MLX',
	'ALR':          'ALR_POS',
	'Mean':         'MEANLUM_ART_MAGS',
	'P1':           'P01_ALL_MAGS',
	'P50':          'P50_ALL_MAGS',
	'P99':          'P99_ALL_MAGS',
	'Star':         'VISSTARS_PCT',
	'SQM_syn':      'SYN_SQM',
}

#CCD-only column, with no fisheye equivalent: a real (measured, not
#synthetic) SQM reading, listed as an extra column at the end if present
#in the CCD database
CCD_ONLY_COLUMN = 'SQM'

#Text color (ARGB hex, no leading '#') for fisheye vs. CCD derived value
#columns in the output .xlsx
FISHEYE_TEXT_COLOR = 'FFCC5500'  #dark orange
CCD_TEXT_COLOR = 'FF595959'	  #dark gray

#Columns to center-align in the output .xlsx
CENTER_ALIGN_COLS = ['CCD Dset', 'Date', 'CCD Mid Obs Time (LMT)',
					 'Fisheye Obs Time (LMT)', 'Time Difference (min)'] + \
					['F '+m for m in CCD_COLUMNS] + ['C '+m for m in CCD_COLUMNS] + \
					['C '+CCD_ONLY_COLUMN]

#Columns to force to display 2 decimal places, even for whole numbers
DECIMAL_COLS = ['F '+m for m in CCD_COLUMNS if m != 'Star'] + \
			   ['C '+m for m in CCD_COLUMNS if m != 'Star'] + \
			   ['C '+CCD_ONLY_COLUMN]

#All metric columns (including Star, which is excluded from DECIMAL_COLS)
#that should share a single uniform column width
METRIC_COLS = ['F '+m for m in CCD_COLUMNS] + ['C '+m for m in CCD_COLUMNS] + \
			  ['C '+CCD_ONLY_COLUMN]

#-----------------------------------------------------------------------------#

def zenith(img, r, a=4.75):
	"""
	Median zenith brightness, identical to metrics.py's zenith() function.

	Parameters
	----------
	img : 2D array
		Fisheye image, masked with the terrain mask, calibrated in
		mag/arcsec^2.
	r : 2D array
		Pixel radius from the image center, same shape as img.
	a : number, optional
		Aperture radius in pixels (~1 degree). Defaults to 4.75.

	Returns
	-------
	zenith_mag : float
		Median zenith brightness [mag/arcsec^2].
	"""
	zenith_mag = round(n.median(img[n.where(r<a)]), 2)
	return zenith_mag


def load_theta_r_calibration(path):
	"""
	Build theta(r) and dtheta/dr(r) interpolators from an empirically
	measured calibration table, identical to metrics.py's
	load_theta_r_calibration() function.

	Parameters
	----------
	path : str
		Path to the theta_r.xlsx calibration file.

	Returns
	-------
	theta_of_r : callable
		Pixel radius [pix] -> zenith angle [radians].
	dtheta_dr_of_r : callable
		Pixel radius [pix] -> local d(theta)/dr [radians/pix].
	R : float
		Maximum calibrated pixel radius (theta=90 degrees).
	"""
	cal = pd.read_excel(path)

	theta_deg = cal['Zenith_Angle (deg)'].values
	r_cumulative = cal['Cumulative_Pixel'].values
	delta_pixel = cal['Delta_Pixel'].values
	step_deg = theta_deg[0]

	theta_deg_full = n.concatenate(([0], theta_deg))
	r_full = n.concatenate(([0], r_cumulative))

	theta_of_r = interp1d(r_full, n.deg2rad(theta_deg_full),
						  bounds_error=False, fill_value=(0, n.nan),
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
	Build per-pixel Theta, dTheta, and dOmega using the calibrated
	theta(r) relationship, identical to metrics.py's
	build_calibrated_geometry() function.

	Parameters
	----------
	r : 2D array
		Pixel radius from the image center.
	theta_of_r, dtheta_dr_of_r : callable
		From load_theta_r_calibration().
	R : float
		Maximum calibrated pixel radius (theta=90 degrees); pixels beyond
		this are outside the field of view and masked to NaN.

	Returns
	-------
	Theta : 2D array
		Zenith angle [radians]. NaN beyond the field of view.
	dOmega : 2D array
		Solid angle [sr] per pixel: dOmega = sin(Theta)*dTheta/r (limit
		dTheta^2 at r=0). NaN beyond the field of view.
	"""
	Theta = theta_of_r(r)
	dTheta = dtheta_dr_of_r(r)

	dOmega = n.full_like(r, n.nan, dtype=float)
	on_axis = (r == 0)
	off_axis = ~on_axis
	dOmega[off_axis] = n.sin(Theta[off_axis]) * dTheta[off_axis] / r[off_axis]
	dOmega[on_axis] = dTheta[on_axis]**2

	beyond_fov = r > R
	Theta[beyond_fov] = n.nan
	dOmega[beyond_fov] = n.nan

	return Theta, dOmega


def mag_to_ucd(mag):
	"""
	Convert mag/arcsec^2 to luminance in microcandela per square meter,
	identical to metrics.py's mag_to_ucd() function.

	Parameters
	----------
	mag : float or array
		Brightness in mag/arcsec^2.

	Returns
	-------
	L : float or array
		Luminance [ucd/m^2].
	"""
	return 108.48*n.exp(20.7233-0.92104*mag)


def illuminance_horizontal(img, Theta, dOmega):
	"""
	Horizontal illuminance, identical to metrics.py's
	illuminance_horizontal() function.

	Parameters
	----------
	img : 2D array
		Fisheye image, masked with the 90-degree horizon mask, calibrated
		in mag/arcsec^2.
	Theta : 2D array
		Zenith angle [radians] at each pixel, same shape as img.
	dOmega : 2D array
		Solid angle [sr] per pixel, same shape as img.

	Returns
	-------
	E_h : float
		Horizontal illuminance [mlx].
	"""
	L = mag_to_ucd(img) #[ucd m-2]
	dE = L*n.cos(Theta)*dOmega
	E_h = n.nansum(dE)/1000 #Horizontal illuminance [mlx]
	return E_h


def precompute_vertical_geometry(Phi, step_deg=5):
	"""
	Precompute the azimuth sweep grid and cos(incidence) weighting array,
	identical to metrics.py's precompute_vertical_geometry() function.

	Parameters
	----------
	Phi : 2D array
		Azimuth [radians] at each pixel, compass convention (0=North,
		clockwise E/S/W).
	step_deg : number, optional
		Azimuth sweep interval in degrees. Defaults to 5.

	Returns
	-------
	azimuths_deg : 1D array
		Swept azimuths, 0 to 355 degrees.
	cos_incidence : 3D array
		cos(Phi-facing_azimuth) per swept azimuth, clipped >= 0. Shape:
		(n_az, ny, nx).
	"""
	azimuths_deg = n.arange(0, 360, step_deg)
	facing_azimuths = n.deg2rad(azimuths_deg)

	cos_incidence = n.cos(Phi[n.newaxis,:,:] - facing_azimuths[:,n.newaxis,n.newaxis])
	cos_incidence = n.clip(cos_incidence, 0, None)

	return azimuths_deg, cos_incidence


def illuminance_vertical_max(img, Theta, dOmega, azimuths_deg, cos_incidence):
	"""
	Maximum vertical illuminance over all swept facing directions,
	identical to metrics.py's illuminance_vertical_max() function.
	Matches the MAXVERT_MLX metric in the CCD database.

	Parameters
	----------
	img : 2D array
		Fisheye image, masked with the 90-degree horizon mask, calibrated
		in mag/arcsec^2.
	Theta : 2D array
		Zenith angle [radians] at each pixel, same shape as img.
	dOmega : 2D array
		Solid angle [sr] per pixel, same shape as img.
	azimuths_deg : 1D array
		From precompute_vertical_geometry().
	cos_incidence : 3D array
		From precompute_vertical_geometry().

	Returns
	-------
	E_v_max : float
		Maximum vertical illuminance [mlx].
	best_azimuth_deg : float
		Azimuth (compass, 0=North) at which the maximum occurs.
	"""
	L = mag_to_ucd(img) #[ucd m-2]

	base = L*n.sin(Theta)*dOmega
	dE = base[n.newaxis,:,:] * cos_incidence
	E_v_sweep = n.nansum(dE, axis=(1,2)) / 1000

	E_v_max = round(float(E_v_sweep.max()), 2)
	best_azimuth_deg = azimuths_deg[n.argmax(E_v_sweep)]

	return E_v_max, best_azimuth_deg


def ALR_and_mean(img, dOmega, natural_reference=NATURAL_REFERENCE):
	"""
	All-sky Light Pollution Ratio (ALR) and mean all-sky brightness,
	identical to metrics.py's ALR_and_mean() function.

	Parameters
	----------
	img : 2D array
		Fisheye image [mag/arcsec^2], terrain-masked.
	dOmega : 2D array
		Solid angle [sr] per pixel.
	natural_reference : number, optional
		Natural reference luminance [ucd/m^2]. Defaults to 250.

	Returns
	-------
	alr : float
		All-sky Light Pollution Ratio (dimensionless).
	mean_mag : float
		Solid-angle-weighted mean sky brightness [mag/arcsec^2].
	"""
	L = mag_to_ucd(img) #[ucd m-2]

	dOmega_masked = n.where(n.isnan(img), n.nan, dOmega)
	b_average = n.nansum(L*dOmega_masked)/n.nansum(dOmega_masked) #[ucd m-2]

	alr = round((b_average - natural_reference)/natural_reference, 2)
	mean_mag = round((20.7233 - n.log(b_average/108.48)) / 0.92104, 2)

	return alr, mean_mag


def sky_brightness_percentiles(img, dOmega, percentiles=(1, 50, 99)):
	"""
	Solid-angle-weighted sky brightness percentiles, identical to
	metrics.py's sky_brightness_percentiles() (sq_deg_areas omitted here,
	since only the fixed percentiles P1/P50/P99 are needed for this
	comparison).

	Parameters
	----------
	img : 2D array
		Fisheye image [mag/arcsec^2], terrain-masked.
	dOmega : 2D array
		Solid angle [sr] per pixel.
	percentiles : tuple of float, optional
		Cumulative sky-area percentiles, faintest to brightest. Defaults
		to (1, 50, 99).

	Returns
	-------
	values : dict
		Percentile -> brightness value [mag/arcsec^2], or None if no
		valid pixels exist.
	"""
	dOmega_masked = n.where(n.isnan(img), n.nan, dOmega)

	valid = ~n.isnan(img) & ~n.isnan(dOmega_masked)
	if not n.any(valid):
		return {p: None for p in percentiles}

	b = img[valid]
	w = dOmega_masked[valid]

	order = n.argsort(-b) #faintest first (descending mag/arcsec^2)
	b_sorted = b[order]
	w_sorted = w[order]
	cum_frac = n.cumsum(w_sorted) / n.sum(w_sorted) * 100

	values = {}
	for p in percentiles:
		values[p] = round(float(n.interp(p, cum_frac, b_sorted)), 2)

	return values


def synthetic_sqm(img, Theta, dOmega, sqm_v_offset=SQM_V_OFFSET):
	"""
	Synthetic Sky Quality Meter (SQM) reading, identical to metrics.py's
	synthetic_sqm() function.

	Parameters
	----------
	img : 2D array
		Fisheye image [mag/arcsec^2], masked with the 90-degree horizon
		mask.
	Theta : 2D array
		Zenith angle [radians] per pixel (angle from boresight, SQM
		pointed at zenith).
	dOmega : 2D array
		Solid angle [sr] per pixel.
	sqm_v_offset : number, optional
		SQM-V conversion factor [mag/arcsec^2]. Defaults to 0.07.

	Returns
	-------
	sqm_mag : float or None
		Synthetic SQM reading [mag/arcsec^2], or None if no valid pixels
		contribute.
	"""
	L = mag_to_ucd(img) #[ucd m-2]

	D = sqm_response(n.rad2deg(Theta))
	weight = n.where(n.isnan(img), n.nan, D*dOmega)

	denom = n.nansum(weight)
	if denom == 0 or n.isnan(denom):
		return None

	L_avg = n.nansum(L*weight) / denom #[ucd m-2]
	V_avg = (20.7233 - n.log(L_avg/108.48)) / 0.92104 #V-band-equivalent
	sqm_mag = round(V_avg + sqm_v_offset, 2)

	return sqm_mag


def get_fisheye_metrics(dataset, fname, camera, C, theta_of_r, dtheta_dr_of_r, R_cal):
	"""
	Compute all fisheye metrics for a single image, using the same
	geometry and math as metrics.py. Zenith, ALR, Mean, and percentiles
	use the dataset's terrain mask (mask.fit); horizontal/vertical
	illuminance and the synthetic SQM reading use a simple 90-degree
	horizon mask instead. The calibrated theta(r) relationship
	(theta_r.xlsx) is used throughout.

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
	theta_of_r, dtheta_dr_of_r : callable
		From load_theta_r_calibration().
	R_cal : float
		Maximum calibrated pixel radius (theta=90 degrees).

	Returns
	-------
	metrics : dict
		Maps each key in CCD_COLUMNS to its fisheye-computed value, or
		None throughout if the image could not be read.
	"""
	data_cal = DATA_ROOT+dataset+'/'
	fpath = data_cal+fname

	try:
		with fits.open(fpath, uint=False, memmap=False) as hdul:
			raw_img = hdul[0].data.astype(float, copy=True)
	except Exception as e:
		print('  Could not read %s: %s' %(fpath, e))
		return {key: None for key in CCD_COLUMNS}

	xc = C['Xcenter'][camera]
	yc = C['Ycenter'][camera]
	R = C['Radius'][camera]

	ny, nx = raw_img.shape
	x, y = n.meshgrid(n.arange(nx), n.arange(ny))
	r = n.sqrt((x-xc)**2 + (y-yc)**2)
	Phi = -n.arctan2(y-yc, x-xc) + n.pi/2

	#terrain mask, for zenith/ALR/Mean/percentiles
	try:
		with fits.open(data_cal+'mask.fit', uint=False, memmap=False) as hdul:
			terrain_mask = hdul[0].data.astype(float, copy=True)
	except Exception as e:
		print('  Could not read mask.fit for %s: %s' %(dataset, e))
		terrain_mask = n.ones_like(raw_img)

	img_terrain = raw_img * terrain_mask

	#simple 90-degree horizon mask, for illuminance/SQM
	horizon_mask = n.ones_like(raw_img)
	horizon_mask[r > R] = n.nan
	img_horizon = raw_img * horizon_mask

	Theta, dOmega = build_calibrated_geometry(r, theta_of_r, dtheta_dr_of_r, R_cal)
	azimuths_deg, cos_incidence = precompute_vertical_geometry(Phi, step_deg=VERTICAL_STEP_DEG)

	zenith_mag = zenith(img_terrain, r)

	horiz_illum = illuminance_horizontal(img_horizon, Theta, dOmega)
	horiz_illum = round(float(horiz_illum), 2) if horiz_illum is not None else None

	vert_illum_max, _ = illuminance_vertical_max(
		img_horizon, Theta, dOmega, azimuths_deg, cos_incidence)

	alr, mean_mag = ALR_and_mean(img_terrain, dOmega)

	star_pct = round(float(predict_visibility(alr, warn_outside_range=False))) \
			   if alr is not None and alr > 0 else None

	sqm_mag = synthetic_sqm(img_horizon, Theta, dOmega)

	pct = sky_brightness_percentiles(img_terrain, dOmega)

	return {
		'Zenith': zenith_mag,
		'Horizontal': horiz_illum,
		'Vertical Max': vert_illum_max,
		'ALR': alr,
		'Mean': mean_mag,
		'P1': pct[1],
		'P50': pct[50],
		'P99': pct[99],
		'Star': star_pct,
		'SQM_syn': sqm_mag,
	}


def save_with_autofit_columns(df, outpath, fisheye_cols=None, ccd_cols=None, uniform_width_cols=None):
	"""
	Write a DataFrame to an .xlsx file with each column's width set to fit
	its content, using openpyxl, with a minimum width equal to the width
	of the 'Date' column. Columns listed in fisheye_cols have their text
	colored dark orange; columns listed in ccd_cols have their text
	colored dark gray. Columns listed in uniform_width_cols are all set to
	the same width (the widest among them), so metric columns line up
	visually regardless of individual content length.

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
	uniform_width_cols : list of str, optional
		Column names to force to a single shared width (the max natural
		width among them).
	"""
	if fisheye_cols is None:
		fisheye_cols = []
	if ccd_cols is None:
		ccd_cols = []
	if uniform_width_cols is None:
		uniform_width_cols = []

	fisheye_font = Font(color=FISHEYE_TEXT_COLOR)
	ccd_font = Font(color=CCD_TEXT_COLOR)
	center_align = Alignment(horizontal='center', vertical='center')

	if 'Date' in df.columns:
		date_values = df['Date'].tolist()
		min_width = max([len('Date')] + [len(str(v)) for v in date_values]) + 2
	else:
		min_width = 10

	#compute the natural width of each column once, then find the max
	#among uniform_width_cols so they can all share that single width
	natural_widths = {}
	for col in df.columns:
		values = df[col].tolist()
		max_len = max([len(str(v)) for v in values]) if values else 0
		natural_widths[col] = max(max_len+2, min_width)

	if uniform_width_cols:
		shared_width = max(natural_widths[col] for col in uniform_width_cols if col in natural_widths)
		shared_width = shared_width / 2
	else:
		shared_width = None

	with pd.ExcelWriter(outpath, engine='openpyxl') as writer:
		df.to_excel(writer, index=False, sheet_name='Sheet1')
		ws = writer.sheets['Sheet1']

		for i, col in enumerate(df.columns, start=1):
			col_letter = get_column_letter(i)
			values = df[col].tolist()

			if col in uniform_width_cols:
				ws.column_dimensions[col_letter].width = shared_width
			else:
				ws.column_dimensions[col_letter].width = natural_widths[col]

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
	Computes all fisheye metrics for every row in
	concurrent_observations.xlsx, looks up the matching CCD metrics from
	the CCD database, and writes the combined table. See the script
	description for detail.
	"""

	#--------------------------------------------------------------------------#
	#							  Load input tables							   #
	#--------------------------------------------------------------------------#
	print('Loading concurrent observations from %s' %CONCURRENT_OBS_PATH)
	obs = pd.read_excel(CONCURRENT_OBS_PATH)

	#Date may be read in as a full Timestamp; keep only the date portion
	if 'Date' in obs.columns:
		obs['Date'] = pd.to_datetime(obs['Date']).dt.date

	print('Loading fisheye database from %s' %DATA_SUMMARY_PATH)
	log = pd.read_excel(DATA_SUMMARY_PATH, sheet_name='Log')

	#Dataset -> camera name lookup (e.g. 'Fish5'), from the Log sheet's
	#Camera column ('Fish<N>' convention used in imagecenter.csv)
	cameras = {row['Dataset']: 'Fish'+str(row['Camera']).strip()
			   for _, row in log.iterrows()}

	print('Loading CCD database from %s' %NPMAPS_PATH)
	npmaps = pd.read_excel(NPMAPS_PATH)

	#DSET may be read as text or as a float (e.g. 2.0 if the column had
	#any blank/NaN values, which forces pandas to upcast the whole column
	#to float64); normalize both sides to a clean integer-like string so
	#"2", "2.0", and 2 all match consistently
	def clean_dset(x):
		try:
			return str(int(float(x)))
		except (ValueError, TypeError):
			return str(x).strip()

	npmaps['DSET'] = npmaps['DSET'].apply(clean_dset)
	obs['CCD Dset'] = obs['CCD Dset'].apply(clean_dset)

	print('Loading fisheye camera geometry from %simagecenter.csv' %CALIBRATION)
	C = pd.read_csv(CALIBRATION+'imagecenter.csv', index_col=0)

	print('Loading theta(r) calibration from %stheta_r.xlsx' %CALIBRATION)
	theta_of_r, dtheta_dr_of_r, R_cal = load_theta_r_calibration(CALIBRATION+'theta_r.xlsx')

	#--------------------------------------------------------------------------#
	#		  Compute fisheye metrics and look up matching CCD metrics	   #
	#--------------------------------------------------------------------------#
	fisheye_results = {key: [] for key in CCD_COLUMNS}
	ccd_results = {key: [] for key in CCD_COLUMNS}

	#only track the CCD-only real SQM column if it actually exists in
	#this CCD database
	has_ccd_sqm = CCD_ONLY_COLUMN in npmaps.columns
	if has_ccd_sqm:
		ccd_sqm_values = []

	for _, row in obs.iterrows():

		dataset = row['Fisheye Dataset']
		fname = row['Fisheye File']
		ccd_dnight = row['CCD Dataset']
		ccd_dset = row['CCD Dset']

		print('Processing %s / %s ...' %(dataset, fname))

		#fisheye metrics
		camera = cameras.get(dataset)
		if camera is None:
			print('  No camera specified for %s, skipping.' %dataset)
			metrics = {key: None for key in CCD_COLUMNS}
		else:
			metrics = get_fisheye_metrics(
				dataset, fname, camera, C, theta_of_r, dtheta_dr_of_r, R_cal)
		for key in CCD_COLUMNS:
			fisheye_results[key].append(metrics[key])

		#CCD metrics -- match on both DNIGHT and DSET, since a single
		#DNIGHT can contain multiple dsets/sites
		ccd_match = npmaps[(npmaps['DNIGHT'] == ccd_dnight) &
							(npmaps['DSET'] == ccd_dset)]
		if ccd_match.empty:
			print('  No CCD match found for %s dset %s' %(ccd_dnight, ccd_dset))
			for key in CCD_COLUMNS:
				ccd_results[key].append(None)
			if has_ccd_sqm:
				ccd_sqm_values.append(None)
		else:
			rec = ccd_match.iloc[0]
			for key, ccd_col in CCD_COLUMNS.items():
				ccd_results[key].append(rec.get(ccd_col, None))
			if has_ccd_sqm:
				ccd_sqm_values.append(rec.get(CCD_ONLY_COLUMN, None))

	for key in CCD_COLUMNS:
		obs['F '+key] = fisheye_results[key]
		obs['C '+key] = ccd_results[key]

	#CCD-only real SQM reading, appended as the last column if present
	if has_ccd_sqm:
		obs['C '+CCD_ONLY_COLUMN] = ccd_sqm_values

	#ensure fisheye values are rounded to 2 decimal places in the output,
	#even if float precision drifted during DataFrame assignment
	for key in CCD_COLUMNS:
		if key == 'Star':
			obs['F Star'] = obs['F Star'].round(0)
			obs['C Star'] = obs['C Star'].round(0)
		else:
			obs['F '+key] = obs['F '+key].round(2)

	#round the time difference to the nearest whole minute
	if 'Time Difference (min)' in obs.columns:
		obs['Time Difference (min)'] = obs['Time Difference (min)'].round(0).astype('Int64')

	#--------------------------------------------------------------------------#
	#								Save results							   #
	#--------------------------------------------------------------------------#
	save_with_autofit_columns(
		obs, OUTPUT_PATH,
		fisheye_cols=['F '+key for key in CCD_COLUMNS],
		ccd_cols=['C '+key for key in CCD_COLUMNS] + (['C '+CCD_ONLY_COLUMN] if has_ccd_sqm else []),
		uniform_width_cols=METRIC_COLS)
	print('\nSaved metric comparison table to', OUTPUT_PATH)
	print(obs)


if __name__ == '__main__':
	main()