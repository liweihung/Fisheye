#Datasets
datanight = 'BOSE_20231210'
camera = 'Fish2'
processor = 'L Hung'

#Folders (do not change)
calibration = '../Calibration/'
data_cal = '../Data_processed/'+datanight+'/'
data_raw = '../Data_raw/'+datanight+'/'

#Calibration files
linearity = calibration + 'linearity_fish2_20230201.txt'
flatv     = calibration + 'flat_v_fish2_20240131.fit'
flatb     = calibration + ''
mask      = calibration + 'mask_Fish2_1179_780_756.fit'

#Measure zeropoint, extinction coefficient, and center on the reference image?
measure_reference = True			#[True/False] Solve the reference image?
reference = 'Light*0003.fit'		#Reference image
apikey = 'cllxijkpvxsibace' 		#Astrometry API key

#Select zeropoint to use
use_default_zeropoint = False	#[True/False] If False, use measured zeropoint

#File names (do not change)
fn_linearity = linearity[len(calibration):]
fn_mask = mask[len(calibration):]
