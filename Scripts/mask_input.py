#input file and fitting parameters for fitting the fisheye mask (mask.py)

#file
folder = '../Calibration/'
filein = folder+'flat_v_fish2_20210624.fit' 	#input file: usually a flat
fileout = filein[:-4]+'_mask.fit'				#output file: fisheye mask

#brightness cutoff:
t = 0.05		 #bright/light pixel threshold 
