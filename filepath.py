#Datasets
datanight = 'PNGL20200427'
linearity_curve = 'ASI_6200_Linearity.txt'
flatv = 'flat_4x4_rebinned_20200623.fit'
flatb = ''
processor = 'L_Hung'

#File folders
calibration = 'c:/users/lhung/Research/Monitoring/Calibration_files/'
data_cal = 'c:/users/lhung/Research/Monitoring/Data_calibrated/'+datanight+'/'
data_raw = 'c:/users/lhung/Research/Monitoring/Data_raw/'+datanight+'/'

#Calibration file paths
linearity = calibration+linearity_curve
flatV = calibration+flatv
mask = calibration+flatv[:-4]+'_mask.fit'

#Astrometry API key
apikey = 'cllxijkpvxsibace'
