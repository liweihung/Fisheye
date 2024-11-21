#Datasets
datanight = ''  # ROMO_20241003 for example
camera = ''     # Fish5 for example
processor = ''  # L Hung for example

#Update headers
update_headers = #[True/False] 
longitude =      # in decimal degrees, negative sign for west
latitude  =      # in decimal degrees
park_name = ''   # full name, eg. 'Rocky Mountain National Park'
site_name = ''   # full name, eg. 'Parking Lot A'
observers = ''   # follow this format: 'B Banet, J White, L Hung'

#Calibration files
Linearity = ''
Flat_V = ''
Flatb_B = ''
Mask = ''

#Measure zeropoint, extinction coefficient, and center on the reference image?
measure_reference = 	#[True/False] Solve the reference image?
reference = ''		    #Reference image for example Light*V*3.fit
apikey = '' 		    #Astrometry API key cllxijkpvxsibace

#Select zeropoint to use
use_default_zeropoint = 	#[True/False] If False, use measured zeropoint


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


