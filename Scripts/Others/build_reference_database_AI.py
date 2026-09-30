#------------------------------------------------------------------------------#
#build_reference_database_AI.py
#
#NPS Night Skies Program
#
#Builds and maintains a master database of reference-image night sky
#brightness metrics, aggregated from the per-dataset metrics.xlsx files
#produced by metrics.py.
#
#For every dataset under DATA_ROOT, this script:
#   (1) reads its metrics.xlsx,
#   (2) selects the single row flagged Reference == 'Yes',
#   (3) attaches provenance (dataset name, camera, processor, site, lat/lon,
#       reference image pattern, and the time this database row was built),
#   (4) writes the combined table to three outputs that are always kept in
#       sync:
#         - reference_metrics.db   (SQLite -- the system of record)
#         - reference_metrics.csv  (plain text, git-diffable, opens anywhere)
#         - reference_metrics.xlsx (formatted, for collaborators)
#
#The SQLite table is keyed on the dataset name, so re-running after
#reprocessing a dataset cleanly replaces that dataset's row rather than
#duplicating it. Because the per-dataset metrics.xlsx files are themselves
#the source of truth, the table is rebuilt in full on each run.
#
#Provenance is read from the process_input.py that main.py copies into each
#dataset's processed-data folder. If that copy is missing, the row is still
#written using whatever provenance is available from metrics.xlsx.
#
#Datasets with no metrics.xlsx, no Reference == 'Yes' row, or more than one
#Reference == 'Yes' row are reported and handled gracefully (see below).
#
#Input:
#   (1) DATA_ROOT+'*/metrics.xlsx'      -- per-dataset metrics (metrics.py)
#   (2) DATA_ROOT+'*/process_input.py'  -- per-dataset provenance (optional)
#
#Output:
#   (1) OUTPUT_DIR+'reference_metrics.db'   -- SQLite database
#   (2) OUTPUT_DIR+'reference_metrics.csv'  -- CSV view
#   (3) OUTPUT_DIR+'reference_metrics.xlsx' -- Excel view
#
#History:
#	Li-Wei Hung -- Created
#
#------------------------------------------------------------------------------#
import importlib.util
import os
import sqlite3
import sys
from datetime import datetime, timezone
from glob import glob

import pandas as pd

from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

#------------------------------------------------------------------------------#
#                              Configuration                                   #
#------------------------------------------------------------------------------#
#Root folder containing each dataset's processed data, e.g.
#DATA_ROOT+'ANDE_20250915/metrics.xlsx'. This script lives in
#Scripts/Others/, so an extra '../' is needed to reach Data_processed/
#relative to the Scripts/ folder. Adjust to your layout.
DATA_ROOT = '../../Data_processed/'

#Where the master database and its companion views are written (one level
#above Scripts/, i.e. the project root).
OUTPUT_DIR = '../../'

#Base name (without extension) for the three synced outputs.
DB_BASENAME = 'reference_metrics'

#SQLite table name.
TABLE_NAME = 'reference_metrics'

#Provenance fields to pull from each dataset's process_input.py, if present.
#Keyed by the attribute name in process_input.py; the value is the column
#name to use in the database (renamed where needed to avoid collisions --
#e.g. process_input's 'reference' image glob would clash with metrics.py's
#'Reference' Yes/No flag, since SQLite column names are case-insensitive).
PROVENANCE_FIELDS = {
    'camera':     'camera',
    'processor':  'processor',
    'longitude':  'longitude',
    'latitude':   'latitude',
    'park_name':  'park_name',
    'site_name':  'site_name',
    'observers':  'observers',
    'reference':  'reference_pattern',
}

#Columns to force to 2 decimal places in the Excel view.
DECIMAL_COLS = ['Zenith', 'P1', 'P50', 'P95', 'P99', 'Bright1', 'Bright025',
                'Mean', 'ALR', 'Horizontal', 'Vertical', 'SQM_syn']

#------------------------------------------------------------------------------#


def load_process_input(dataset_dir):
	"""
	Load the process_input.py copied into a dataset's processed-data folder
	and return the requested provenance fields as a dict.

	The module is loaded under a unique name each call so that datasets do
	not clobber one another in sys.modules.

	Parameters
	----------
	dataset_dir : str
		Path to the dataset's processed-data folder.

	Returns
	-------
	dict
		Maps each database column name in PROVENANCE_FIELDS to its value
		(or None if the process_input.py or the source attribute is
		missing).
	"""
	prov = {col: None for col in PROVENANCE_FIELDS.values()}
	pi_path = os.path.join(dataset_dir, 'process_input.py')
	if not os.path.exists(pi_path):
		return prov

	try:
		mod_name = 'pi_' + os.path.basename(os.path.normpath(dataset_dir))
		spec = importlib.util.spec_from_file_location(mod_name, pi_path)
		module = importlib.util.module_from_spec(spec)
		spec.loader.exec_module(module)
		for attr, col in PROVENANCE_FIELDS.items():
			prov[col] = getattr(module, attr, None)
	except Exception as e:
		print('  Warning: could not read provenance from %s (%s)' % (pi_path, e))

	return prov


def get_reference_row(metrics_path, dataset):
	"""
	Read a dataset's metrics.xlsx and return its reference-image row.

	Parameters
	----------
	metrics_path : str
		Path to the dataset's metrics.xlsx.
	dataset : str
		Dataset name, used only in warning messages.

	Returns
	-------
	dict or None
		The reference row as a dict, or None if the file has no row flagged
		Reference == 'Yes'. If multiple rows are flagged, the first is used
		and a warning is printed.
	"""
	df = pd.read_excel(metrics_path)

	if 'Reference' not in df.columns:
		print('  Warning: %s has no Reference column, skipping.' % dataset)
		return None

	ref = df[df['Reference'].astype(str).str.strip().str.lower() == 'yes']

	if ref.empty:
		print('  Warning: %s has no reference image (Reference == Yes), skipping.' % dataset)
		return None
	if len(ref) > 1:
		print('  Warning: %s has %d reference rows; using the first.' % (dataset, len(ref)))

	return ref.iloc[0].to_dict()


def build_dataframe():
	"""
	Aggregate the reference-image metrics from every dataset under DATA_ROOT
	into a single DataFrame, with provenance columns prepended.

	Returns
	-------
	DataFrame
		One row per dataset that has a valid reference-image metrics row,
		indexed by dataset name.
	"""
	metrics_files = sorted(glob(os.path.join(DATA_ROOT, '*', 'metrics.xlsx')))

	if not metrics_files:
		print('No metrics.xlsx files found under %s' % DATA_ROOT)
		return pd.DataFrame()

	build_time = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
	rows = []

	for metrics_path in metrics_files:
		dataset_dir = os.path.dirname(metrics_path)
		dataset = os.path.basename(os.path.normpath(dataset_dir))
		print('Processing %s ...' % dataset)

		try:
			ref_row = get_reference_row(metrics_path, dataset)
		except Exception as e:
			print('  Warning: could not read %s (%s), skipping.' % (metrics_path, e))
			continue
		if ref_row is None:
			continue

		prov = load_process_input(dataset_dir)

		#provenance first, then the metrics, then bookkeeping
		row = {'dataset': dataset}
		row.update(prov)
		row.update(ref_row)
		row['source_file'] = metrics_path
		row['db_build_time'] = build_time
		rows.append(row)

	if not rows:
		print('No reference-image rows found in any dataset.')
		return pd.DataFrame()

	df = pd.DataFrame(rows)

	#every row in this database is a reference image by construction, so the
	#per-image 'Reference' Yes/No flag from metrics.xlsx is redundant here
	#(and its name collides with the 'reference' provenance field under
	#SQLite's case-insensitive column naming); drop it
	df = df.drop(columns=['Reference'], errors='ignore')

	#stable, readable column order: dataset + provenance, then everything
	#else in its existing order, with bookkeeping columns last
	prov_cols = list(PROVENANCE_FIELDS.values())
	lead = ['dataset'] + prov_cols
	tail = ['source_file', 'db_build_time']
	#keep only columns that actually exist, preserve order, avoid duplicates
	ordered = [c for c in lead if c in df.columns]
	ordered += [c for c in df.columns if c not in lead + tail]
	ordered += [c for c in tail if c in df.columns]
	df = df[ordered].set_index('dataset').sort_index()

	return df


def save_sqlite(df, path):
	"""Write the DataFrame to a SQLite table (rebuilds the table each run)."""
	conn = sqlite3.connect(path)
	try:
		df.to_sql(TABLE_NAME, conn, if_exists='replace', index=True)
	finally:
		conn.close()
	print('Saved SQLite database to %s (table: %s)' % (path, TABLE_NAME))


def save_csv(df, path):
	"""Write the DataFrame to a CSV view."""
	df.to_csv(path, index=True)
	print('Saved CSV view to %s' % path)


def save_excel(df, path):
	"""
	Write the DataFrame to an Excel view with columns auto-sized to content,
	all cells center-justified, and DECIMAL_COLS formatted to 2 decimals.
	"""
	df_out = df.reset_index()  #make 'dataset' a visible column

	with pd.ExcelWriter(path, engine='openpyxl') as writer:
		df_out.to_excel(writer, index=False, sheet_name='Reference Metrics')
		ws = writer.sheets['Reference Metrics']
		center = Alignment(horizontal='center', vertical='center')
		header_font = Font(bold=True)

		for i, col in enumerate(df_out.columns, start=1):
			col_letter = get_column_letter(i)
			values = df_out[col].tolist()
			max_len = max([len(str(col))] + [len(str(v)) for v in values])
			ws.column_dimensions[col_letter].width = max(max_len + 2, 10)

			ws['%s1' % col_letter].font = header_font
			for row in range(1, len(values) + 2):
				ws['%s%d' % (col_letter, row)].alignment = center

			if col in DECIMAL_COLS:
				for row in range(2, len(values) + 2):
					ws['%s%d' % (col_letter, row)].number_format = '0.00'

	print('Saved Excel view to %s' % path)


def main():
	"""
	Aggregate reference-image metrics from all datasets and write the SQLite
	database and its CSV and Excel companion views. See the script
	description for detail.
	"""
	df = build_dataframe()
	if df.empty:
		print('Nothing to write; exiting.')
		return

	os.makedirs(OUTPUT_DIR, exist_ok=True)
	base = os.path.join(OUTPUT_DIR, DB_BASENAME)

	save_sqlite(df, base + '.db')
	save_csv(df, base + '.csv')
	save_excel(df, base + '.xlsx')

	print('\nDatabase built: %d datasets.' % len(df))


if __name__ == '__main__':
	main()