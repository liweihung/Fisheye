#input file and fitting parameters for fitting the fisheye mask (mask.py)

#file
folder = 'c:/users/lhung/Research/Monitoring/Calibration_files/'
filein = folder+'SIgma_3.5_Master_Flat.fit' 	#input file: usually a flat
fileout = filein[:-4]+'_mask.fit'				#output file: fisheye mask

#brightness cutoff:
t = 3000.			 #bright/light pixel threshold 
