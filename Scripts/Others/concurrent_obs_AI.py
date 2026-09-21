#-----------------------------------------------------------------------------#
#list_bestpick_observing_times.py
#
#NPS Night Skies Program
#
#Last updated: 2026/08/27
#
#This script lists, for each hand-picked fisheye/CCD "best pick" concurrent
#comparison pair, the CCD dataset name, fisheye dataset name, fisheye image
#file, park/site name, observing date (from the fisheye observation), CCD
#mid-obs time (LMT), fisheye observing time (LMT, time-of-day only, to match
#the CCD time format), and the time difference between the two. Only rows
#where the fisheye and CCD observing times differ by MAX_TIME_DIFF_MIN
#minutes or less are kept, and only the top N_BEST_PER_NIGHT entries with
#the smallest time difference are kept per night of observation (grouped by
#park + date, so multi-part datasets like INDU_20240802A/B or
#GUIS_20230315A/B are treated as one night). Fisheye local time is derived
#from each image's FITS header DATE-OBS (UTC) using the same
#TimezoneFinder/pytz approach as projection.py's get_local_time_from_utc().
#Output rows are sorted by CCD dataset name.
#
#This script lives in Scripts/Others/, one level deeper than the other
#pipeline scripts in Scripts/, so its relative paths climb one extra level.
#
#Note: reading both the fisheye and CCD .xlsx databases requires the
#openpyxl package.
#
#Input:
#   (1) ../../Performance/Concurrent_observations/
#       Data_summary_fisheye_20260827.xlsx -- fisheye database (Log sheet),
#       for Park, Site, Latitude, Longitude, and Dataset name
#	(2) ../../Performance/Concurrent_observations/
#	    Data_summary_CCD_modified.xlsx -- CCD database, for MID_TIME_LMT
#	    and MID_DATE_LMT
#	(3) Calibrated fisheye images: data_cal+'Light*.fit' and
#	    data_cal+'img*sky*.fit' for each dataset listed in DATASETS
#
#Output:
#   (1) ../../Performance/Concurrent_observations/concurrent_observations.xlsx
#       -- one row per fisheye image (top N_BEST_PER_NIGHT per park+date),
#       with the matching CCD dataset info, both observing times, and
#       their difference in minutes
#
#History:
#	Li-Wei Hung -- Created
#
#-----------------------------------------------------------------------------#
import pandas as pd
import pytz

from astropy.io import fits
from datetime import datetime, timedelta
from glob import glob
from openpyxl.styles import Alignment
from openpyxl.utils import get_column_letter
from timezonefinder import TimezoneFinder
import os

#-----------------------------------------------------------------------------#
#                        Dataset pairs to compare                             #
#-----------------------------------------------------------------------------#
#Each entry: fisheye dataset folder name, matching CCD DNIGHT, and CCD DSET
#number (since a single DNIGHT can contain multiple dsets/sites). These are
#the "best pick" concurrent, co-located pairs identified from the fisheye
#and CCD databases below.

DATASETS = [
	#fisheye_dataset,   ccd_dnight,     ccd_dset
	('BOSE_20221020',   'USFW221021',   1),
	('BOSE_20231210A',  'USFW231211',   1),
	('BOSE_20231210B',  'USFW231211',   1),
	('CEBR_20230714',   'CEBR230715_1', 2),
	('GUIS_20221019',   'GUIS221020',   2),
	('GUIS_20221021A',  'GUIS221022A',  2),
	('GUIS_20221021B',  'GUIS221022C',  2),
	('GUIS_20230315A',  'GUIS230316',   1),
	('GUIS_20230315B',  'GUIS230316',   1),
	('NIOB_20240830A',  'NIOB240830',   4),
	('NIOB_20240830C',  'NIOB240831',   3),
	('ROMO_20241003B',  'ROMO241004',   1),
	('WHSA_20220227',   'WHSA220227',   3),
]

#-----------------------------------------------------------------------------#
#                              File locations                                 #
#-----------------------------------------------------------------------------#
#Root folder containing each dataset's processed data, e.g.
#DATA_ROOT+'NIOB_20240830C/'. This script lives in Scripts/Others/, so an
#extra '../' is needed to reach Data_processed/ relative to the Scripts/
#folder. Adjust to match your directory layout.
DATA_ROOT = '../../Data_processed/'

#Fisheye database (Log sheet) -- for Park, Site, Latitude, and Longitude.
#Reading this requires the openpyxl package.
DATA_SUMMARY_PATH = '../../Performance/Concurrent_observations/Data_summary_fisheye_20260827.xlsx'

#CCD database -- for MID_TIME_LMT and MID_DATE_LMT. Reading this requires
#the openpyxl package.
NPMAPS_PATH = '../../Performance/Concurrent_observations/Data_summary_CCD_modified.xlsx'

#Output file
OUTPUT_PATH_XLSX = '../../Performance/Concurrent_observations/concurrent_observations.xlsx'

#Only keep rows whose fisheye/CCD observing times differ by no more than
#this many minutes
MAX_TIME_DIFF_MIN = 180

#Only keep this many entries with the smallest time difference, per night of
#observation (grouped by park + date, since e.g. INDU_20240802A/B and
#GUIS_20230315A/B are the same night split across multiple datasets)
N_BEST_PER_NIGHT = 3

#Specific fisheye image files to exclude from consideration entirely,
#keyed by dataset. File names must match exactly (including extension) as
#they appear in the dataset's data_cal folder.
EXCLUDED_FILES = {
	'GUIS_20221021A': {
		'img-0003-sky-V.fit',
		'img-0005-sky-V.fit',
	},
	'NIOB_20240830C': {
		'Light_V_%d.fit' %i for i in range(1, 11)
	},
	'ROMO_20241003B': {
		'Light_RainbowROMO_30.0s_Bin4_V_20241004-003742_0001.fit',
		'Light_RainbowROMO_30.0s_Bin4_V_20241004-004224_0003.fit',
	},
	'WHSA_20220227': {
		'img-0005-sky-V.fit',
		'img-0006-sky-V.fit',
		'img-0014-sky-V.fit',
	},
}

#Columns to center-align in the output .xlsx
CENTER_ALIGN_COLS = ['CCD Dset', 'Date', 'CCD Mid Obs Time (LMT)',
					 'Fisheye Obs Time (LMT)', 'Time Difference (min)']

#-----------------------------------------------------------------------------#

def get_local_time_from_utc(lat, lon, utc_time):
	"""
	Convert a UTC datetime to local mean time given a latitude/longitude.
	Identical approach to projection.py's get_local_time_from_utc().

	Parameters
	----------
	lat, lon : float
		Latitude and longitude in decimal degrees.
	utc_time : datetime
		Naive datetime interpreted as UTC.

	Returns
	-------
	local_time : datetime or None
		Timezone-aware local datetime, or None if the timezone could not be
		determined for the given coordinates.
	"""
	tf = TimezoneFinder()
	timezone_str = tf.timezone_at(lng=lon, lat=lat)
	if timezone_str is None:
		return None
	timezone = pytz.timezone(timezone_str)
	utc_time = utc_time.replace(tzinfo=pytz.utc)
	local_time = utc_time.astimezone(timezone)
	return local_time


def lmt_decimal_to_hhmm(t):
	"""
	Convert a decimal-hour local mean time (as used in the CCD database's
	MID_TIME_LMT column, e.g. 23.09) to an "HH:MM" string. Values may
	exceed 24 to indicate times after local midnight relative to the
	observing date.

	Parameters
	----------
	t : float or None
		Decimal-hour local mean time.

	Returns
	-------
	str or None
		Time formatted as "HH:MM", or None if t is missing.
	"""
	if t is None or pd.isna(t):
		return None
	th = float(t) % 24
	hh = int(th)
	mm = int(round((th-hh)*60))
	if mm == 60:
		mm = 0
		hh = (hh+1) % 24
	return '%02d:%02d' %(hh, mm)


def ccd_datetime(ccd_date_str, ccd_time_decimal):
	"""
	Reconstruct a full datetime for the CCD observation from its date
	(MID_DATE_LMT) and decimal-hour local mean time (MID_TIME_LMT, e.g.
	23.09). Decimal hours >= 24 roll over to the next calendar day,
	matching the convention of recording post-midnight times relative to
	the start of the observing night.

	MID_DATE_LMT may arrive either as a string (e.g. '8/30/24', when the
	source spreadsheet stores dates as text) or as a pandas Timestamp
	(when Excel auto-detects the column as a date), depending on how the
	source file was created; both forms are handled here.

	Parameters
	----------
	ccd_date_str : str, datetime, or pandas.Timestamp
		Observing date, either as an 'M/D/YY' or 'M/D/YYYY' string, or as
		an already-parsed datetime/Timestamp.
	ccd_time_decimal : float
		Decimal-hour local mean time.

	Returns
	-------
	dt : datetime or None
		Reconstructed local datetime, or None if inputs are missing/invalid.
	"""
	if ccd_date_str is None or pd.isna(ccd_date_str):
		return None
	if ccd_time_decimal is None or pd.isna(ccd_time_decimal):
		return None

	if isinstance(ccd_date_str, (pd.Timestamp, datetime)):
		base_date = datetime(ccd_date_str.year, ccd_date_str.month, ccd_date_str.day)
	else:
		try:
			base_date = datetime.strptime(str(ccd_date_str).strip(), '%m/%d/%y')
		except ValueError:
			try:
				base_date = datetime.strptime(str(ccd_date_str).strip(), '%m/%d/%Y')
			except ValueError:
				return None

	total_hours = float(ccd_time_decimal)
	day_offset = int(total_hours // 24)
	hours_in_day = total_hours % 24
	hh = int(hours_in_day)
	mm = int(round((hours_in_day-hh)*60))
	if mm == 60:
		mm = 0
		hh = (hh+1) % 24

	dt = base_date.replace(hour=hh, minute=mm)
	if day_offset:
		dt += timedelta(days=day_offset)

	return dt


def get_fisheye_files(dataset):
	"""
	Return all fisheye image files for a dataset, matching either of the two
	naming conventions used across this pipeline, excluding any files
	listed in EXCLUDED_FILES for this dataset.

	Parameters
	----------
	dataset : str
		Fisheye dataset folder name, e.g. 'NIOB_20240830C'.

	Returns
	-------
	files : list of str
		Sorted list of matching file paths, with excluded files removed.
	"""
	data_cal = DATA_ROOT+dataset+'/'
	files = glob(data_cal+'Light*.fit') + glob(data_cal+'img*sky*.fit')

	excluded = EXCLUDED_FILES.get(dataset, set())
	if excluded:
		files = [f for f in files if os.path.basename(f) not in excluded]

	return sorted(files)


def save_with_autofit_columns(df, outpath):
	"""
	Write a DataFrame to an .xlsx file with each column's width set to fit
	its content (based on the longest value in the column, including the
	header), using openpyxl, capped at the width of the 'Fisheye Dataset'
	column so no column exceeds it. Columns listed in CENTER_ALIGN_COLS are
	center-justified.

	Parameters
	----------
	df : DataFrame
		Table to write.
	outpath : str
		Destination .xlsx file path.
	"""
	center = Alignment(horizontal='center', vertical='center')

	#column width cap, based on the 'Fisheye Dataset' column's own content
	if 'Fisheye Dataset' in df.columns:
		max_width = max([len('Fisheye Dataset')] +
						[len(str(v)) for v in df['Fisheye Dataset']]) + 2
	else:
		max_width = None

	with pd.ExcelWriter(outpath, engine='openpyxl') as writer:
		df.to_excel(writer, index=False, sheet_name='Sheet1')
		ws = writer.sheets['Sheet1']

		for i, col in enumerate(df.columns, start=1):
			col_letter = get_column_letter(i)

			width = max([len(str(col))] + [len(str(v)) for v in df[col]]) + 2
			if max_width is not None:
				width = min(width, max_width)
			ws.column_dimensions[col_letter].width = width

			if col in CENTER_ALIGN_COLS:
				for row in range(1, len(df)+2):
					ws[f'{col_letter}{row}'].alignment = center


def main():
	"""
	Lists, for each best-pick dataset pair, every fisheye image's observing
	time (converted to LMT) alongside the matching CCD dataset's mid-obs
	time (LMT), park/site names, observing date, and the time difference
	between the two. See the script description for detail.
	"""

	#--------------------------------------------------------------------------#
	#							  Load reference databases					   #
	#--------------------------------------------------------------------------#
	print('Loading fisheye database from %s' %DATA_SUMMARY_PATH)
	log = pd.read_excel(DATA_SUMMARY_PATH, sheet_name='Log')

	print('Loading CCD database from %s' %NPMAPS_PATH)
	npmaps = pd.read_excel(NPMAPS_PATH)

	#DSET is sometimes read as text (e.g. "1" instead of 1) depending on how
	#the source spreadsheet stores it; normalize both sides to string so the
	#DATASETS list's integer dset numbers still match correctly
	npmaps['DSET'] = npmaps['DSET'].astype(str).str.strip()

	#--------------------------------------------------------------------------#
	#				 Build one row per fisheye image, per dataset			   #
	#--------------------------------------------------------------------------#
	rows = []
	for dataset, ccd_dnight, ccd_dset in DATASETS:

		print('Processing %s ...' %dataset)

		#fisheye Log entry -- for Park, Site, Latitude, Longitude
		log_match = log[log['Dataset'] == dataset]
		if log_match.empty:
			print('  No Log entry found for %s, skipping.' %dataset)
			continue
		log_rec = log_match.iloc[0]
		park = log_rec.get('Park', None)
		site = log_rec.get('Site', None)
		lat = log_rec.get('Latitude', None)
		lon = log_rec.get('Longitude', None)

		#CCD record -- for MID_TIME_LMT and MID_DATE_LMT
		ccd_match = npmaps[(npmaps['DNIGHT'] == ccd_dnight) &
							(npmaps['DSET'] == str(ccd_dset))]
		if ccd_match.empty:
			print('  No CCD match found for %s dset %s' %(ccd_dnight, ccd_dset))
			ccd_time_lmt = None
			ccd_dt = None
		else:
			ccd_rec = ccd_match.iloc[0]
			ccd_time_decimal = ccd_rec.get('MID_TIME_LMT', None)
			ccd_time_lmt = lmt_decimal_to_hhmm(ccd_time_decimal)
			ccd_dt = ccd_datetime(ccd_rec.get('MID_DATE_LMT', None), ccd_time_decimal)

		#all fisheye image files for this dataset
		files = get_fisheye_files(dataset)
		if len(files) == 0:
			print('  No fisheye image files found in %s%s/' %(DATA_ROOT, dataset))
			continue

		for f in files:
			try:
				with fits.open(f, uint=False, memmap=False) as hdul:
					date_obs = hdul[0].header.get('DATE-OBS', None)
			except Exception as e:
				print('  Could not read header of %s: %s' %(f, e))
				date_obs = None

			#parse the UTC timestamp from the header, trying both the
			#with- and without-fractional-seconds formats
			utc_time = None
			for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S'):
				try:
					utc_time = datetime.strptime(date_obs, fmt)
					break
				except (TypeError, ValueError):
					continue

			#convert to local mean time using the site's lat/lon
			fisheye_lmt = None
			fisheye_local_dt = None
			if utc_time is not None and lat is not None and lon is not None:
				local_time = get_local_time_from_utc(lat, lon, utc_time)
				if local_time is not None:
					fisheye_local_dt = local_time.replace(tzinfo=None)
					fisheye_lmt = local_time.strftime('%H:%M')

			#time difference between CCD mid-obs time and this fisheye image,
			#in minutes (positive = fisheye is later than CCD)
			time_diff_min = None
			if fisheye_local_dt is not None and ccd_dt is not None:
				time_diff_min = round(
					(fisheye_local_dt-ccd_dt).total_seconds()/60, 1)

			#only keep entries within the allowed time difference
			if time_diff_min is None or abs(time_diff_min) > MAX_TIME_DIFF_MIN:
				continue

			#observing date taken from the fisheye observation itself
			fisheye_date = fisheye_local_dt.date() if fisheye_local_dt else None

			rows.append({
				'CCD Dataset': ccd_dnight,
				'CCD Dset': ccd_dset,
				'Fisheye Dataset': dataset,
				'Fisheye File': f[len(DATA_ROOT+dataset+'/'):],
				'Park': park,
				'Site': site,
				'Date': fisheye_date,
				'CCD Mid Obs Time (LMT)': ccd_time_lmt,
				'Fisheye Obs Time (LMT)': fisheye_lmt,
				'Time Difference (min)': time_diff_min,
			})

	#--------------------------------------------------------------------------#
	#		 Keep only the top N_BEST_PER_NIGHT entries per park+date			   #
	#--------------------------------------------------------------------------#
	df_all = pd.DataFrame(rows)
	if not df_all.empty:
		df_all['AbsTimeDiff'] = df_all['Time Difference (min)'].abs()
		df_all.sort_values('AbsTimeDiff', inplace=True)
		df_all = df_all.groupby(['Park', 'Date'], as_index=False, sort=False) \
						.head(N_BEST_PER_NIGHT)
		df_all.drop(columns='AbsTimeDiff', inplace=True)

		#rank entries by CCD dataset name
		df_all.sort_values('CCD Dataset', inplace=True)
		df_all.reset_index(drop=True, inplace=True)

	#--------------------------------------------------------------------------#
	#								Save results							   #
	#--------------------------------------------------------------------------#
	save_with_autofit_columns(df_all, OUTPUT_PATH_XLSX)
	print('\nSaved observing times table to', OUTPUT_PATH_XLSX)
	print(df_all)


if __name__ == '__main__':
	main()