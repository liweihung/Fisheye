#Datasets
datanight = ''  # 'ROMO_20241003' for example
camera = ''     # 'Fish5' for example
processor = ''  # 'L Hung' for example

#Update headers
update_headers = False  # [True/False]
longitude = -0.0         # in decimal degrees, negative sign for west
latitude  = 0.0         # in decimal degrees
park_name = ''   # full name, eg. 'Rocky Mountain National Park'
site_name = ''   # full name, eg. 'Parking Lot A'
observers = ''   # follow this format: 'B Banet, J White, L Hung'

#Photometric band(s) to process
bands = ('V',)   # which bands to process this run

#Calibration files
Linearity = ''
Flat_V = ''
Flat_B = ''
Flat_R = ''
Mask = ''

#Measure zeropoint, extinction coefficient, and center on the reference image?
measure_reference = False  # [True/False] Solve the reference image?
reference = ''		    #Reference image for example 'Light*V*3.fit'
apikey = '' 		    #Astrometry API key 

#Select zeropoint to use
use_default_zeropoint = True  # [True/False] If False, use measured zeropoint


#--------------------------------------------------------------------------#
#      Folders and filepaths (do not change anything below this line)      #
#--------------------------------------------------------------------------#
calibration  = '../Calibration/'
data_cal     = '../Data_processed/'+datanight+'/'
data_raw     = '../Data_raw/'+datanight+'/'

linearity    = calibration + Linearity
flatv        = calibration + Flat_V
flatb        = calibration + Flat_B
flatr        = calibration + Flat_R
mask         = calibration + Mask

fn_linearity = linearity[len(calibration):]
fn_mask      = mask[len(calibration):]

# Bands that the pipeline is actually capable of calibrating
supported_bands = ('V',)
