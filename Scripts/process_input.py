#Datasets
datanight = 'CORI_20240925E'
camera = 'Fish5'
processor = 'L Hung'

#Folders (do not change)
calibration = '../Calibration/'
data_cal = '../Data_processed/'+datanight+'/'
data_raw = '../Data_raw/'+datanight+'/'

#Calibration files
linearity = calibration + 'linearity_fish5_20241009.txt'
flatv     = calibration + 'flat_v_fish5_20241016.fit'
flatb     = calibration + ''
mask      = calibration + 'mask_Fish5_1178_775_760.fit'

#Measure zeropoint, extinction coefficient, and center on the reference image?
measure_reference = True			#[True/False] Solve the reference image?
reference = 'Light*V*3.fit'		#Reference image
apikey = 'cllxijkpvxsibace' 		#Astrometry API key

#Select zeropoint to use
use_default_zeropoint = False	#[True/False] If False, use measured zeropoint

#File names (do not change)
fn_linearity = linearity[len(calibration):]
fn_mask = mask[len(calibration):]
