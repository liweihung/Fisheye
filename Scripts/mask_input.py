#input file and fitting parameters for fitting the fisheye mask (mask.py)

#file
folder = 'c:/users/lhung/Research/Fisheye/Calibration/'
filein = folder+'flat_2x2_20200623.fit' 	#input file: usually a flat
fileout = filein[:-4]+'_mask.fit'				#output file: fisheye mask

#brightness cutoff:
t = 0.5			 #bright/light pixel threshold 
