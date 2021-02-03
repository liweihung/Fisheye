#input file and fitting parameters for fitting the fisheye mask (mask.py)

#file
folder = 'c:/users/lhung/Research/Fisheye/Calibration/'
filein = folder+'Sigma_3.5_Master_Flat.fit' 	#input file: usually a flat
fileout = filein[:-4]+'_mask.fit'				#output file: fisheye mask

#brightness cutoff:
t = 14000		 #bright/light pixel threshold 
