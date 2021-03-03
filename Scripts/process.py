#-----------------------------------------------------------------------------#
#main.py
#
#NPS Night Skies Program
#
#Last updated: 2020/10/13
#
#This is the main pipeline for processing fisheye images. If measure_reference 
#is True, photometric and positional calibration constants will be derived from 
#the reference image. If not, the default values will be used.
#
#Input: 
#   (1) scripts listed under the local sources
#
#Output:
#   (1) see the documentation for the individual script called
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
# Local Source
import astrometry
import process_input
import median_filter
import photometric_calibration
import photometry
import positional_calibration
import projection
import reduction

#------------------------------------------------------------------------------#
reduction.main()						#data reduction
astrometry.main()						#astrometric measurement
photometry.main()						#Gaussian PSF photometry
photometric_calibration.main()			#absolute photometric calibration
positional_calibration.main()			#positional calibration
median_filter.main()					#median filter
projection.main()						#plotting

	
	
