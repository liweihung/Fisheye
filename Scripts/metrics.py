#------------------------------------------------------------------------------#
#metrics.py
#
#NPS Night Skies Program
#
#Extracts night sky brightness metrics from calibrated fisheye images.
#
#Zenith brightness, ALR, Mean, and percentiles use the terrain mask
#(mask.fit), reflecting only the actually-visible sky. Horizontal and
#vertical illuminance use a simple 90-degree horizon mask instead,
#reflecting the full hemisphere.
#
#Input:
#   (1) theta_r.xlsx -- empirically measured theta(r) calibration table
#	(2) mask.fit -- terrain mask for this dataset
#
#Output:
#   (1) p.data_cal+'metrics.xlsx' -- Date, Time (LMT), File name, Zenith,
#       percentiles, Mean, ALR, Star (% of naked-eye-visible stars), SQM
#       (synthetic Sky Quality Meter reading), Horizontal/Vertical
#       Illuminance, azimuth, plus each brightness metric converted to
#       ucd/m^2, for every processed image
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
from star_visibility_model_AI import predict_visibility
from sqm_model_AI import sqm_response

#------------------------------------------------------------------------------#

def get_obs_date_time_lmt(hdr, file):
    """
	Extract observing date and local mean time from a FITS header.

	Parameters
	----------
	hdr : astropy FITS header
		Must contain DATE-OBS, SITELAT, SITELONG.
	file : str
		File name, used only in the warning if DATE-OBS is unparseable.

	Returns
	-------
	date : str or None
		Observing date (YYYY-MM-DD), or None if unparseable.
	hour : str or None
		Observing time (HH:MM) in local mean time, or None if unparseable.
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
	Median zenith brightness over a circular aperture at image center.

	Corresponds to "Zenith" in the online NPMaps database.

	Parameters
	----------
	img : 2D array
        Fisheye image, masked with the terrain mask.
	a : number, optional
        Aperture radius in pixels (~1 degree), comparable to CCD zenith
        aperture. 1 pix = 379 arcsec.

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
	measured calibration table, replacing the equidistant fisheye
	assumption (theta = pi/2 * r/R, constant dtheta/dr).

	theta_r.xlsx columns: Zenith_Angle (deg) | Delta_Pixel | Cumulative_Pixel,
	sampled every 3 degrees from 3 to 90. Cumulative_Pixel is r(theta);
	Delta_Pixel is the pixel distance covered by each 3-degree step.

	Parameters
	----------
	path : str
		Path to the theta_r.xlsx calibration file.

	Returns
	-------
	theta_of_r : callable
		Pixel radius [pix] -> zenith angle [radians], valid for r in [0, R].
	dtheta_dr_of_r : callable
		Pixel radius [pix] -> local d(theta)/dr [radians/pix].
	R : float
		Maximum calibrated pixel radius (theta=90 degrees).
	"""
	cal = pd.read_excel(path)

	theta_deg = cal['Zenith_Angle (deg)'].values
	r_cumulative = cal['Cumulative_Pixel'].values
	delta_pixel = cal['Delta_Pixel'].values
	step_deg = theta_deg[0] #assumes a uniform step, e.g. 3 degrees

	#prepend the origin (r=0 at theta=0), since the table starts at 3 deg
	theta_deg_full = n.concatenate(([0], theta_deg))
	r_full = n.concatenate(([0], r_cumulative))

	#theta(r)
	theta_of_r = interp1d(r_full, n.deg2rad(theta_deg_full),
						  bounds_error=False, fill_value=(0, n.nan),
						  kind='linear')

	#dtheta/dr(r): step_deg/Delta_Pixel per band, interpolated at each
	#band's midpoint radius
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
	Build per-pixel Theta, dTheta, and dOmega using the calibrated theta(r)
	relationship.

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
	dTheta : 2D array
		Local angular step [radians/pix].
	dOmega : 2D array
		Solid angle [sr] per pixel: dOmega = sin(Theta)*dTheta/r, i.e.
		sin(theta)*dtheta*dphi with dphi=1/r (one pixel of tangential arc
		at radius r). At r=0 (removable singularity) the limit dTheta^2 is
		used instead. NaN beyond the field of view.
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

	return Theta, dTheta, dOmega


def mag_to_ucd(mag):
    """
	Convert mag/arcsec^2 to luminance in microcandela per square meter,
	using the same conversion applied throughout this script (Duriscoe
	2016 equation 1).

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


def illuminance_horizontal(img):
    """
	Horizontal illuminance from a fisheye image.

	Parameters
	----------
	img : 2D array
        Fisheye image, masked with the 90-degree horizon mask (not the
        terrain mask), so illuminance reflects the full hemisphere.

	Returns
	-------
	E_h : float
		Horizontal illuminance [mlx].
	"""
    L = mag_to_ucd(img) # [ucd m-2]
    dE = L*n.cos(Theta)*dOmega
    E_h = round(n.nansum(dE)/1000, 2)
    return E_h


def precompute_vertical_geometry(Phi, step_deg=5):
    """
	Precompute the azimuth sweep grid and cos(incidence) weighting used by
	illuminance_vertical_max(), once per dataset (Phi is fixed by camera
	geometry).

	Parameters
	----------
	Phi : 2D array
		Azimuth [radians], compass convention (0=North, clockwise E/S/W).
	step_deg : number, optional
		Azimuth sweep interval in degrees. Defaults to 5.

	Returns
	-------
	azimuths_deg : 1D array
		Swept azimuths, 0 to 355 degrees.
	cos_incidence : 3D array
		cos(Phi-facing_azimuth) per swept azimuth, clipped >= 0 to exclude
		light from behind the surface. Shape: (n_az, ny, nx).
	"""
    azimuths_deg = n.arange(0, 360, step_deg)
    facing_azimuths = n.deg2rad(azimuths_deg)

    cos_incidence = n.cos(Phi[n.newaxis,:,:] - facing_azimuths[:,n.newaxis,n.newaxis])
    cos_incidence = n.clip(cos_incidence, 0, None)

    return azimuths_deg, cos_incidence


def illuminance_vertical_max(img, azimuths_deg, cos_incidence):
    """
	Vertical illuminance at every swept azimuth, and its maximum. Matches
	the MAXVERT_MLX metric in the CCD database.

	Parameters
	----------
	img : 2D array
		Fisheye image, masked with the 90-degree horizon mask.
	azimuths_deg : 1D array
		From precompute_vertical_geometry().
	cos_incidence : 3D array
		From precompute_vertical_geometry().

	Returns
	-------
	E_v_sweep : 1D array
		Vertical illuminance [mlx] at each azimuth.
	E_v_max : float
		Maximum vertical illuminance [mlx].
	best_azimuth_deg : float
		Azimuth (compass, 0=North) at which the maximum occurs.
	"""
    L = mag_to_ucd(img) #[ucd m-2]

    base = L*n.sin(Theta)*dOmega #per-pixel factor independent of azimuth
    dE = base[n.newaxis,:,:] * cos_incidence
    E_v_sweep = n.nansum(dE, axis=(1,2)) / 1000

    E_v_max = round(E_v_sweep.max(), 2)
    best_azimuth_deg = azimuths_deg[n.argmax(E_v_sweep)]

    return E_v_sweep, E_v_max, best_azimuth_deg


def ALR_and_mean(img, dOmega, natural_reference=250):
    """
	All-sky Light Pollution Ratio (ALR) and mean all-sky brightness,
	computed together since both derive from the same solid-angle-weighted
	average luminance.

	ALR (Duriscoe 2016) is total skyglow brightness divided by the natural
	dark sky reference (250 ucd/m^2). Mean is the same weighted average
	luminance, reported directly in mag/arcsec^2 via the inverse of
	Duriscoe (2016) equation 1, instead of as a ratio.

	Duriscoe's average all-sky brightness (Eq. 2) is a plain pixel mean,
	which only equals a true equal-area sky average in an equal-area
	projection -- not the case here. Each pixel is instead weighted by its
	solid angle dOmega. dOmega is re-masked to match img's NaN pattern
	(terrain-blocked pixels) before summing, so the numerator and
	denominator cover exactly the same pixels.

	Corresponds to "Mean all-sky" (mean_mag) in the online NPMaps database.

	Parameters
	----------
	img : 2D array
		Fisheye image [mag/arcsec^2], already multiplied by the terrain
		mask (NaN where terrain blocks the sky).
	dOmega : 2D array
		Solid angle [sr] per pixel, from build_calibrated_geometry(); NaN
		only beyond the field of view (not terrain-masked).
	natural_reference : number, optional
		Natural reference luminance [ucd/m^2]. Defaults to 250 (Duriscoe
		2016 median natural all-sky brightness, sky+stars).

	Returns
	-------
	alr : float
		All-sky Light Pollution Ratio (dimensionless).
	mean_mag : float
		Solid-angle-weighted mean sky brightness [mag/arcsec^2].
	"""
    L = mag_to_ucd(img) #[ucd m-2]

    #match dOmega's NaN mask to img's, so both sums cover the same pixels
    dOmega_masked = n.where(n.isnan(img), n.nan, dOmega)
    b_average = n.nansum(L*dOmega_masked)/n.nansum(dOmega_masked) #[ucd m-2]

    alr = round((b_average - natural_reference)/natural_reference, 2)
    mean_mag = round((20.7233 - n.log(b_average/108.48)) / 0.92104, 2)

    return alr, mean_mag


def sky_brightness_percentiles(img, dOmega, percentiles=(1, 50, 95, 99), sq_deg_areas=(1, 0.25)):
    """
	Solid-angle-weighted sky brightness percentiles.

	A plain pixel percentile treats every pixel as equally representative
	of the sky, which is wrong here since dOmega varies across the
	calibrated, non-equidistant projection. Pixels are instead sorted by
	brightness and weighted by solid angle, so percentiles reflect
	fraction of sky *area*, not pixel count.

	img is in mag/arcsec^2 (lower = brighter), so percentiles accumulate
	from the faintest pixel toward the brightest: p=95 means "95% of the
	visible sky is darker than this value." High percentiles isolate the
	brightest, smallest-area end of the sky.

	Duriscoe's "brightest square degree" and "brightest quarter square
	degree" are fixed-area thresholds (1 and 0.25 square degrees), not
	fixed percentiles -- his 99.995th/99.999th percentile figures only
	equal those areas because his hemisphere totals a fixed 20,500 square
	degrees. Since terrain can block part of the sky here, those target
	percentiles are instead derived from this image's actual visible area
	(sum of dOmega, converted to square degrees), so "brightest square
	degree" always means the brightest 1 square degree actually present,
	whatever fraction of the hemisphere that turns out to be.

	Corresponds to "Darkest" (P1), "Median" (P50), and "Brightest" (P99)
	in the online NPMaps database.

	Parameters
	----------
	img : 2D array
		Fisheye image [mag/arcsec^2], terrain-masked (NaN where blocked).
	dOmega : 2D array
		Solid angle [sr] per pixel, from build_calibrated_geometry(); NaN
		only beyond the field of view.
	percentiles : tuple of float, optional
		Fixed cumulative sky-area percentiles to compute directly (not
		area-based), faintest to brightest. Defaults to (1, 50, 95, 99).
	sq_deg_areas : tuple of float, optional
		Fixed sky areas [square degrees], brightest-end, to compute
		percentiles for (e.g. 1 -> "brightest square degree"). The target
		percentile for each is derived from this image's actual visible
		area, not assumed to be a fixed percentile. Defaults to (1, 0.25).

	Returns
	-------
	values : dict
		Maps each entry in percentiles to its brightness value
		[mag/arcsec^2], and each entry in sq_deg_areas (keyed as e.g.
		'brightest_1_sqdeg') to its brightness value, or None throughout
		if there are no valid (non-NaN) pixels.
	"""
    labels = list(percentiles) + ['brightest_%s_sqdeg' %a for a in sq_deg_areas]

    #match dOmega's NaN mask to img's, so only terrain-visible pixels
    #contribute to the weighted percentile
    dOmega_masked = n.where(n.isnan(img), n.nan, dOmega)

    valid = ~n.isnan(img) & ~n.isnan(dOmega_masked)
    if not n.any(valid):
        return {label: None for label in labels}

    b = img[valid]
    w = dOmega_masked[valid]

    #sort by brightness descending mag/arcsec^2 (i.e. ascending physical
    #brightness -- faintest pixel first), and build the cumulative
    #solid-angle fraction over this dataset's actual visible sky area,
    #accumulating from faintest to brightest
    order = n.argsort(-b)
    b_sorted = b[order]
    w_sorted = w[order]
    total_sr = n.sum(w_sorted) #total visible solid angle [sr]
    cum_frac = n.cumsum(w_sorted) / total_sr * 100 #percent of visible sky area darker than b_sorted[i]

    values = {}

    #fixed cumulative-area percentiles
    for p in percentiles:
        values[p] = round(float(n.interp(p, cum_frac, b_sorted)), 2)

    #fixed sky-area thresholds (e.g. "brightest 1 square degree"), with
    #the target percentile derived from this image's actual visible area
    #rather than assumed to be a fixed percentile
    total_sqdeg = total_sr * (180/n.pi)**2 #steradians -> square degrees
    for a in sq_deg_areas:
        label = 'brightest_%s_sqdeg' %a
        if a >= total_sqdeg:
            #requested area exceeds the entire visible sky; the threshold
            #is simply the faintest pixel present
            values[label] = round(float(b_sorted[-1]), 2)
        else:
            target_pct = (1 - a/total_sqdeg) * 100
            values[label] = round(float(n.interp(target_pct, cum_frac, b_sorted)), 2)

    return values


def synthetic_sqm(img, Theta, dOmega, sqm_v_offset=0.07):
    """
	Synthetic Sky Quality Meter (SQM) reading, following Cinzano (2005)
	Eq. (1): image luminance integrated over the sky, weighted by the
	SQM's angular sensitivity D(theta) (sqm_response()) and solid angle
	dOmega (which already includes sin(theta), so no separate term is
	needed).

	The weighted average is converted to a V-band-equivalent brightness
	via the inverse of Duriscoe (2016) equation 1, then shifted by
	sqm_v_offset (default 0.07 mag/arcsec^2, the SQM-V conversion factor
	from Cinzano 2005 Sec. 5, averaged over natural (+0.06) and polluted
	(+0.08) sky) to produce a value comparable to a real SQM reading.

	D(theta) tapers to zero beyond ~65 degrees from boresight, so pixels
	near the horizon contribute little regardless of the mask on img.

	Parameters
	----------
	img : 2D array
		Fisheye image [mag/arcsec^2], masked with the 90-degree horizon
		mask (not the terrain mask).
	Theta : 2D array
		Zenith angle [radians] per pixel (angle from boresight, assuming
		the SQM points at zenith), from build_calibrated_geometry().
	dOmega : 2D array
		Solid angle [sr] per pixel, from build_calibrated_geometry().
	sqm_v_offset : number, optional
		SQM-V conversion factor [mag/arcsec^2] added to the V-band
		weighted average. Defaults to 0.07.

	Returns
	-------
	sqm_mag : float or None
		Synthetic SQM reading [mag/arcsec^2], or None if no valid pixels
		contribute.
	"""
    L = mag_to_ucd(img) #[ucd m-2]

    #SQM angular sensitivity at each pixel's zenith angle
    D = sqm_response(n.rad2deg(Theta))

    #match D*dOmega's NaN mask to img's, so both sums cover the same pixels
    weight = n.where(n.isnan(img), n.nan, D*dOmega)

    denom = n.nansum(weight)
    if denom == 0 or n.isnan(denom):
        return None

    L_avg = n.nansum(L*weight) / denom #[ucd m-2]
    V_avg = (20.7233 - n.log(L_avg/108.48)) / 0.92104 #V-band-equivalent [mag/arcsec^2]
    sqm_mag = round(V_avg + sqm_v_offset, 2)

    return sqm_mag


def save_with_autofit_columns(df, outpath, min_width=10, decimal_cols=None):
    """
	Write a DataFrame to .xlsx with columns auto-sized to content, all
	cells center-justified, and decimal_cols formatted to 2 decimal places.

	Parameters
	----------
	df : DataFrame
	outpath : str
		Destination .xlsx path.
	min_width : number, optional
		Minimum column width regardless of content. Defaults to 10.
	decimal_cols : list of str, optional
		Columns to force to 2 decimal places.
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
            nrows = len(values)

            max_len = max([len(str(v)) for v in values]) if values else 0
            ws.column_dimensions[col_letter].width = max(max_len+2, min_width)

            for row in range(1, nrows+2):
                ws[f'{col_letter}{row}'].alignment = center

            if col in decimal_cols:
                for row in range(2, nrows+2):
                    ws[f'{col_letter}{row}'].number_format = '0.00'

#------------------------------------------------------------------------------#


def parse_lat_lon(coord):
    """Parse a lat/lon value that may be a number or a 'D M S' string."""
    if isinstance(coord, (float, int)):
        return float(coord)
    parts = coord.strip().split()
    if len(parts) == 3:
        d, m, s = map(float, parts)
        return d + m/60 + s/3600
    return float(coord)


#------------------------------------------------------------------------------#


#List of files to process
datasets = glob('../Data_processed/WHSA_20220227/process_input.py')
#		   glob('../Data_processed/ROMO_20241003B_testcopy/process_input.py')
#    	   glob'../Data_processed/GUIS_20231210B/process_input.py') #+ \
#		   glob('../Data_processed/EGLI*/process_input.py') + \
#          glob('../Data_processed/BOSE*/process_input.py')

for pi in datasets:
    shutil.copy(pi, 'process_input_metric.py')
    import process_input_metric as p
    importlib.reload(p)

    #Image center coordinates and FoV radius
    C = pd.read_csv('../Calibration/imagecenter.csv',index_col=0)
    xc = C['Xcenter'][p.camera]
    yc = C['Ycenter'][p.camera]
    R  = C['Radius'][p.camera]

    #Terrain mask (for zenith, ALR, Mean, and percentiles)
    with fits.open(p.data_cal+'mask.fit', uint=False, memmap=False) as hdul:
        mask = hdul[0].data.astype(float, copy=True)

    #Image geometry
    image_shape = mask.shape
    x, y = n.meshgrid(n.arange(image_shape[1]), n.arange(image_shape[0]))
    r = n.sqrt((x-xc)**2 + (y-yc)**2)
    Phi = -n.arctan2(y-yc, x-xc) + n.pi/2

    #Simple 90-degree horizon mask (for horizontal/vertical illuminance),
    #reflecting the full hemisphere including terrain foreground
    horizon_mask = n.ones(image_shape)
    horizon_mask[r > R] = n.nan

    #Calibrated, non-equidistant theta(r) geometry
    theta_r_file = p.calibration+'theta_r.xlsx'
    theta_of_r, dtheta_dr_of_r, R_cal = load_theta_r_calibration(path=theta_r_file)
    Theta, dTheta, dOmega = build_calibrated_geometry(r, theta_of_r, dtheta_dr_of_r, R_cal)
    azimuths_deg, cos_incidence = precompute_vertical_geometry(Phi, step_deg=5)

    #Images to process
    images = glob(p.data_cal+'Light*.fit') + glob(p.data_cal+'img*sky*.fit')
    reference_files = {os.path.normpath(f) for f in glob(p.data_cal+p.reference)}

    data_rows = []

    for file in images:

        print(f"Processing {file} ...")
        with fits.open(file, uint=False, memmap=False) as hdul:
            hdr = hdul[0].header
            raw_img = hdul[0].data.astype(float, copy=True)
        img_terrain = raw_img * mask
        img_horizon = raw_img * horizon_mask

        date, hour = get_obs_date_time_lmt(hdr, file)

        zenith_mag = zenith(img_terrain)
        print("Zenith:", zenith_mag)

        horizontal_illum = illuminance_horizontal(img_horizon)
        print("Horizontal Illuminance:", horizontal_illum)

        E_v_sweep, E_v_max, best_azimuth_deg = illuminance_vertical_max(
            img_horizon, azimuths_deg, cos_incidence)
        print(f"Max vertical illuminance is at azimuth {best_azimuth_deg}°: {E_v_max:.2f} mlx")

        alr, mean_mag = ALR_and_mean(img_terrain, dOmega)
        print("ALR:", alr, " Mean:", mean_mag)

        #Percentage of naked-eye-visible stars, predicted from ALR
        star_pct = round(float(predict_visibility(alr, warn_outside_range=False))) if alr is not None and alr > 0 else None
        print("Star:", star_pct)

        #Synthetic SQM reading, from the 90-degree horizon-masked image
        sqm_mag = synthetic_sqm(img_horizon, Theta, dOmega)
        print("SQM:", sqm_mag)

        pct = sky_brightness_percentiles(img_terrain, dOmega)
        print("Percentiles (1/50/95/99/brightest 1 sqdeg/brightest 0.25 sqdeg):",
              pct[1], pct[50], pct[95], pct[99],
              pct['brightest_1_sqdeg'], pct['brightest_0.25_sqdeg'])

        def to_ucd(mag):
            """Convert a mag/arcsec^2 value to ucd/m^2, or None."""
            if mag is None:
                return None
            return round(mag_to_ucd(mag))

        #brightness metrics [mag/arcsec^2] to convert, keyed by output
        #column prefix
        mag_values = {
            'Zenith': zenith_mag,
            'P1': pct[1],
            'P50': pct[50],
            'P95': pct[95],
            'P99': pct[99],
            'Bright1': pct['brightest_1_sqdeg'],
            'Bright025': pct['brightest_0.25_sqdeg'],
            'Mean': mean_mag,
        }
        ucd_values = {key: to_ucd(mag) for key, mag in mag_values.items()}

        row = {
            'Date': date,
            'Time': hour,
            'File': os.path.basename(file),
            'Zenith': zenith_mag,
            'P1': pct[1], #Darkest
            'P50': pct[50], #Median
            'P95': pct[95],
            'P99': pct[99], #Brightest 
            'Bright1': pct['brightest_1_sqdeg'],
            'Bright025': pct['brightest_0.25_sqdeg'],
            'Mean': mean_mag, #Mean all-sky
            'ALR': alr,									
            'Star': star_pct, #% of naked-eye-visible stars, from ALR
            'Horizontal': horizontal_illum, #Horizontal
            'Vertical': E_v_max, #Max Vertical
            'Azimuth': best_azimuth_deg,
            'SQM_syn': sqm_mag, #Synthetic Sky Quality Meter reading
        }

        #ucd values, all together
        for key, ucd in ucd_values.items():
            row[key+'_ucd'] = ucd

        data_rows.append(row)

        #Save the reference image's full azimuth sweep
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

    df_metrics = pd.DataFrame(data_rows)
    save_with_autofit_columns(
        df_metrics, p.data_cal+'metrics.xlsx',
        decimal_cols=['Zenith', 'Horizontal', 'Vertical', 'ALR', 'Mean', 'SQM_syn',
                      'P1', 'P50', 'P95', 'P99', 'Bright1', 'Bright025'])
    print('Saved metrics summary to', p.data_cal+'metrics.xlsx')