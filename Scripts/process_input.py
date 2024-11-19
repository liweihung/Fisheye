#Datasets
datanight = 'ROMO_20241003'
camera = 'Fish5'
processor = 'L Hung'

#Update headers
update_headers = True #[True/False] 
longitude =  -105.663598734252  #in decimal degrees, negative sign for west
latitude  =  40.400346402882  #in decimal degrees
park_name = 'Rocky Mountain National Park'   # full name, eg. 'Rocky Mountain National Park'
site_name = 'Rainbow Curve'   # full name, eg. 'Parking Lot A'
observers = 'B Banet, J White'   # follow this format: 'B Banet, J White, L Hung'

#Calibration files
Linearity = 'linearity_fish5_20241009.txt'
Flat_V = 'flat_v_fish5_20241016.fit'
Flatb_B = ''
Mask = 'mask_Fish5_1178_775_760.fit'

#Measure zeropoint, extinction coefficient, and center on the reference image?
measure_reference = True			#[True/False] Solve the reference image?
reference = 'Light*V*3.fit'		    #Reference image
apikey = 'cllxijkpvxsibace' 		#Astrometry API key

#Select zeropoint to use
use_default_zeropoint = False	#[True/False] If False, use measured zeropoint


#--------------------------------------------------------------------------#
#      Folders and filepaths (do not change anything below this line)      #
#--------------------------------------------------------------------------#
calibration  = '../Calibration/'
data_cal     = '../Data_processed/'+datanight+'/'
data_raw     = '../Data_raw/'+datanight+'/'

linearity    = calibration + Linearity
flatv        = calibration + Flat_V
flatb        = calibration + Flatb_B
mask         = calibration + Mask

fn_linearity = linearity[len(calibration):]
fn_mask      = mask[len(calibration):]


