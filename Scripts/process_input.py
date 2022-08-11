#Datasets
datanight = 'ALAM_20220225'
UTCoffset = -7
camera = 'Fish2'
processor = 'L_Hung'

#Folders (do not change)
calibration = '../Calibration/'
data_cal = '../Data_processed/'+datanight+'/'
data_raw = '../Data_raw/'+datanight+'/'

#Calibration files
linearity = calibration + 'linearity_fish1_20200722.txt'
flatv     = calibration + 'flat_v_fish2_20210624.fit'
flatb     = calibration + ''
mask      = calibration + 'flat_v_fish2_20210624_mask.fit'

#Measure zeropoint, extinction coefficient, and center on the reference image?
measure_reference = False				#[True/False] Solve the reference image?
reference = 'img-0006-sky-V.fit'		#Reference image
apikey = 'cllxijkpvxsibace' 			#Astrometry API key

#Select zeropoint to use
use_default_zeropoint = True	#[True/False] If False, use measured zeropoint

#File names (do not change)
fn_linearity = linearity[len(calibration):]
fn_mask = mask[len(calibration):]
