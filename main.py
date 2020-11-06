#-----------------------------------------------------------------------------#
#main.py
#
#NPS Night Skies Program
#
#Last updated: 2020/10/13
#
#This script is the main pipeline to run to process the fisheye images.
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
import centering
import filepath
import median_filter
import photometric_calibration
import photometry
import projection
import reduction

#------------------------------------------------------------------------------#
reduction.main()

#process the referece image to get the pointing and instrumental zeropoint
if filepath.measure_reference:
	astrometry.main()
	photometry.main()
	
photometric_calibration.main()
centering.main()
median_filter.main()
projection.main()

	
	
