#Datasets
datanight = 'PNGL_20201117A'
lat, long = 40.816692, -104.601119 #observing location [degree]
camera = 'Fish1'
processor = 'L_Hung'

#Folders
calibration = '../Calibration/'
data_cal = '../Data_processed/'+datanight+'/'
data_raw = '../Data_raw/'+datanight+'/'

#Calibration files
linearity = calibration + 'ASI_6200_Linearity.txt'
flatv     = calibration + 'Sigma_3.5_Master_Flat.fit'
flatb     = calibration + ''
mask      = calibration + 'Sigma_3.5_Master_Flat_mask.fit'

#Measure zeropoint, extinction coefficient, and center on the reference image?
measure_reference = False				#[True/False]Slove the reference image?
reference = 'img-005-sky-V.fit'			#Reference image
apikey = 'cllxijkpvxsibace' 			#Astrometry API key

#Select zeropoint to use
use_default_zeropoint = False	#[True/False] If False, use measured zeropoint

#File names (do not change)
fn_linearity = linearity[len(calibration):]
fn_mask = mask[len(calibration):]
