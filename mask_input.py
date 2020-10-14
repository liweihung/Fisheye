#input file and fitting parameters for fitting the fisheye mask (mask.py)

#file
folder = 'c:/users/lhung/Research/Monitoring/Calibration_files/'
#filein = folder+'flat_2x2_20200623.fit'
filein = folder+'flat_4x4_rebinned_20200623.fit' #input file: usually a flat
fileout = filein[:-4]+'_mask.fit'				 #output file: fisheye mask

#fitting parameters:
t = 0.1 			 #bright/light pixel threshold 
p0 = [1186,773,747]  #initial guess: x_center, y_center, radius 
s = 1e-1 			 #stepsize
