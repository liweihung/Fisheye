#Datasets
datanight = 'PNGL20200427'
lat, long = 40.696034, -104.600301 #observing location [degree]
camera = 'Fish1'
fn_linearity = 'ASI_6200_Linearity.txt'
fn_flatv = 'flat_4x4_rebinned_20200623.fit'
fn_flatb = ''
fn_mask = 'flat_4x4_rebinned_20200623_mask.fit'
fn_zeropoint = 'zeropoint.txt'
use_default_zeropoint = True	#[True/False] If False, use measured zeropoint

processor = 'L_Hung'

#Measure zeropoint, extinction coefficient, and center on the reference image?
measure_reference = False				#[True/False]Slove the reference image?
reference = '40_sec_V_light_4x4.fit'	#Reference image
apikey = 'cllxijkpvxsibace' 			#Astrometry API key

#File folders
calibration = 'c:/users/lhung/Research/Monitoring/Calibration_files/'
data_cal = 'c:/users/lhung/Research/Monitoring/Data_calibrated/'+datanight+'/'
data_raw = 'c:/users/lhung/Research/Monitoring/Data_raw/'+datanight+'/'

#Calibration file paths
linearity = calibration+fn_linearity
flatv = calibration+fn_flatv
flatb = calibration+fn_flatb
mask = calibration+fn_mask
if use_default_zeropoint:
	zeropoint = calibration+fn_zeropoint
else: 
	zeropoint = data_cal+fn_zeropoint




