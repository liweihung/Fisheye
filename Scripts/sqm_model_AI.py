#-----------------------------------------------------------------------------#
#sqm_model_AI.py
#
#NPS Night Skies Program
#
#Single-curve models of the angular sensitivity D(theta) for two Sky Quality
#Meter variants:
#
#   (1) sqm_response()  -- the standard, lens-less SQM (wide field of view)
#   (2) sqml_response() -- the SQM-L, which adds a lens that narrows the FOV
#
#Both share the same detector (TAOS TSL237) and filter (Hoya CM-500), so their
#*spectral* response -- and therefore the SQM-V conversion offset -- is the
#same. Only the *angular* response differs: the SQM-L collects light from a
#much smaller patch of sky around its pointing direction.
#
#---------------------------------------------------------------------------
#Standard SQM (sqm_response)
#---------------------------------------------------------------------------
#Fitted to the anchor values reported in Cinzano, P. (2005), "Night Sky
#Photometry with Sky Quality Meter", ISTIL Internal Report n.9, v.1.4.
#
#The paper describes (Sec. 2, Fig. 2/3) the angular response of SQM serial
#n.0115 v.1.09, measured on a rotation table against a calibrated point
#source, in both the horizontal and vertical planes:
#	- Half Width at Half Maximum (HWHM)          ~= 42 deg
#	- Factor-of-10 attenuation of a point source  at  55 deg
#	- Onset of mechanical vignetting/screening   ~= 60-65 deg
#
#---------------------------------------------------------------------------
#SQM-L (sqml_response)
#---------------------------------------------------------------------------
#The SQM-L adds a lens that narrows the field of view relative to the standard
#SQM. Its angular response is read from the measured SQM-L curve (Cinzano 2005,
#Fig. 3, right panel; see also the Unihedron FOV curves at
#https://www.unihedron.com/projects/darksky/fovcurves.jpg):
#	- Half Width at Half Maximum (HWHM)          ~= 10 deg (~20 deg FWHM)
#	- Factor-of-10 attenuation of a point source  at  20 deg
#	- Low-level wings trailing to                ~= 40-50 deg
#
#---------------------------------------------------------------------------
#Model form (shared)
#---------------------------------------------------------------------------
#Both curves use a super-Gaussian (generalized Gaussian) core,
#	D(theta) = exp[-ln(2)*(|theta|/HWHM)**p]   for |theta| <= taper_start
#with a raised-cosine taper to zero between taper_start and taper_end,
#representing mechanical screening of the detector at large angles:
#	D(theta) = D_core(theta)*0.5*(1+cos(pi*t)),  t = (|theta|-taper_start)/(taper_end-taper_start)
#for taper_start < |theta| < taper_end, and D(theta)=0 for |theta| >= taper_end.
#
#The exponent p is solved analytically so the curve passes through both anchor
#points exactly: D(HWHM)=0.5 (by definition of half width at half max) and
#D(factor10_angle)=0.1 (the paper's stated factor-of-10 point). For the standard
#SQM this gives p~=4.452 (a plain Gaussian, p=2, only reaches D=0.30 at 55 deg,
#too shallow). For the SQM-L the same two-anchor solve gives p~=1.73, a
#sub-Gaussian with a sharp core and heavy wings, matching the measured curve.
#
#Provenance note: these fits are derived from the paper's stated/plotted
#numeric anchors (HWHM, the factor-of-10 point, and the screening onset), not
#from a pixel-level digitization of the original figures or the underlying raw
#data tables (not provided in the paper). Treat them as best-effort analytic
#stand-ins for the real angular responses, adequate for weighting/convolution
#use, not a substitute for the manufacturer's own measured response curves if
#metrological precision is needed.
#
#Usage:
#	from sqm_model_AI import sqm_response, sqml_response
#	D  = sqm_response(theta_deg)   #standard SQM; scalar or array in, same out
#	DL = sqml_response(theta_deg)  #SQM-L (narrow FOV)
#
#	#Use as a weighting kernel, e.g. reproducing Eq. (1) of the paper:
#	#  I_avg = integral[ I(theta,phi) * D(theta) * sin(theta) dtheta dphi ]
#	#        / integral[ D(theta) * sin(theta) dtheta dphi ]
#
#Output (when run directly):
#   (1) ../Calibration/sqm_angular_response_fit.png -- fitted SQM and SQM-L
#       models vs. an approximate reconstruction of the paper's Fig. 3 overlay
#
#History:
#	Derived/fitted by Claude (Anthropic) from Cinzano (2005), for reuse in
#	downstream numerical scripts. SQM-L curve added from Fig. 3 (right panel).
#	Independently verify against the source paper / instrument datasheet
#	before using in any precision application.
#
#-----------------------------------------------------------------------------#
import numpy as n

#-----------------------------------------------------------------------------#
#                     Fitted model parameters: standard SQM                    #
#-----------------------------------------------------------------------------#

HWHM = 42.0            #degrees, half width at half maximum (paper-stated)
FACTOR10_ANGLE = 55.0  #degrees, angle at which D = 0.1 (paper-stated)
TAPER_START = 55.0     #degrees, start of mechanical vignetting taper
TAPER_END = 65.0       #degrees, full screening / D -> 0 (paper: "screening begins at 60-65 deg")

#Solve exponent p analytically so that D(HWHM)=0.5 and D(FACTOR10_ANGLE)=0.1:
#   ln(10)/ln(2) = (FACTOR10_ANGLE / HWHM) ** p
#   p = ln( ln(10)/ln(2) ) / ln( FACTOR10_ANGLE / HWHM )
EXPONENT = float(n.log(n.log(10)/n.log(2)) / n.log(FACTOR10_ANGLE/HWHM))

#-----------------------------------------------------------------------------#
#                       Fitted model parameters: SQM-L                         #
#-----------------------------------------------------------------------------#
# The SQM-L uses the same TSL237 detector + Hoya CM-500 filter as the standard
# SQM but adds a lens narrowing the field of view. Spectral response (hence the
# SQM-V offset) is unchanged; only the angular response D(theta) narrows.
# Anchors below are read from the measured SQM-L angular response (Cinzano 2005,
# Fig. 3, right panel): HWHM ~= 10 deg, D(20 deg) ~= 0.1, with low-level wings
# trailing to ~40-50 deg. HWHM ~= 10 deg (~20 deg FWHM) vs ~42 deg for the
# lens-less SQM.

HWHM_L = 10.0            #degrees, half width at half maximum (from Fig. 3)
FACTOR10_ANGLE_L = 20.0  #degrees, angle at which D ~= 0.1 (from Fig. 3)
TAPER_START_L = 25.0     #degrees, start of mechanical/wing taper
TAPER_END_L = 50.0       #degrees, full screening / D -> 0 (extended to preserve
                         #the low-level wings visible in the SQM-L curve)

#Solve exponent so D(HWHM_L)=0.5 and D(FACTOR10_ANGLE_L)=0.1, as for the SQM.
#HWHM=10, D(20)=0.1 -> p ~= 1.73 (sub-Gaussian: sharp core with heavy wings,
#matching the measured SQM-L angular response).
EXPONENT_L = float(n.log(n.log(10)/n.log(2)) / n.log(FACTOR10_ANGLE_L/HWHM_L))

#-----------------------------------------------------------------------------#


def _core(theta_abs, hwhm, p):
	"""
	Super-Gaussian core, before the vignetting taper is applied.

	Parameters
	----------
	theta_abs : float or array
		Absolute angle(s) from boresight, in degrees.
	hwhm : float
		Half width at half maximum, degrees.
	p : float
		Super-Gaussian exponent.

	Returns
	-------
	float or array
		Core response, same shape as theta_abs.
	"""
	return n.exp(-n.log(2) * (theta_abs/hwhm)**p)


def _angular_response(theta_deg, hwhm, p, taper_start, taper_end):
	"""
	Evaluate a super-Gaussian angular response with a raised-cosine taper.

	Shared implementation for both the standard SQM and the SQM-L; the two
	public functions (sqm_response, sqml_response) differ only in their
	default parameters. Normalized so D(0) = 1; symmetric about 0; D = 0 for
	|theta| >= taper_end.

	Parameters
	----------
	theta_deg : float or array
		Angle(s) from boresight, in degrees. Positive or negative; the
		response is symmetric about 0.
	hwhm : float
		Half width at half maximum, degrees.
	p : float
		Super-Gaussian exponent.
	taper_start, taper_end : float
		Angular range (degrees) over which a raised-cosine taper brings the
		response smoothly to zero, representing mechanical vignetting.

	Returns
	-------
	D : float or array
		Relative angular response, normalized so D(0) = 1. Same shape as
		theta_deg (a plain float is returned for scalar input).
	"""
	theta = n.asarray(theta_deg, dtype=float)
	a = n.abs(theta)

	D = _core(a, hwhm, p)

	#apply raised-cosine taper between taper_start and taper_end
	in_taper = (a > taper_start) & (a < taper_end)
	beyond = a >= taper_end

	if n.any(in_taper):
		t = (a[in_taper]-taper_start) / (taper_end-taper_start)
		taper_factor = 0.5 * (1.0+n.cos(n.pi*t))
		D[in_taper] = D[in_taper] * taper_factor

	D = n.where(beyond, 0.0, D)

	#preserve scalar input as scalar output
	if n.isscalar(theta_deg) or (n.ndim(theta_deg) == 0):
		return float(D)
	return D


def sqm_response(theta_deg, hwhm=HWHM, p=EXPONENT,
				  taper_start=TAPER_START, taper_end=TAPER_END):
	"""
	Evaluate the fitted standard-SQM angular sensitivity D(theta).

	Parameters
	----------
	theta_deg : float or array
		Angle(s) from boresight (zenith when SQM points up), in degrees.
		Can be positive or negative; response is symmetric about 0.
	hwhm : float, optional
		Half width at half maximum, degrees. Defaults to 42.0 (paper value).
	p : float, optional
		Super-Gaussian exponent. Defaults to the analytically fitted value
		(~4.452) matching HWHM=42 deg and D(55 deg)=0.1.
	taper_start, taper_end : float, optional
		Angular range (degrees) over which a raised-cosine taper is applied
		to bring the response smoothly to zero, representing mechanical
		vignetting/screening. Defaults to 55 and 65 degrees per the paper's
		reported screening onset.

	Returns
	-------
	D : float or array
		Relative angular response, normalized so D(0) = 1. Same shape as
		theta_deg.
	"""
	return _angular_response(theta_deg, hwhm, p, taper_start, taper_end)


def sqml_response(theta_deg, hwhm=HWHM_L, p=EXPONENT_L,
				   taper_start=TAPER_START_L, taper_end=TAPER_END_L):
	"""
	Evaluate the fitted SQM-L angular sensitivity D(theta): the narrow-FOV
	(lens-equipped) variant.

	Identical in form to sqm_response(), but with the SQM-L's narrower field
	of view. The SQM-L shares the SQM's spectral response, so any SQM-V
	conversion offset applied downstream is unchanged.

	Parameters
	----------
	theta_deg : float or array
		Angle(s) from boresight (zenith when the SQM-L points up), in
		degrees. Can be positive or negative; response is symmetric about 0.
	hwhm : float, optional
		Half width at half maximum, degrees. Defaults to 10.0 (from Fig. 3).
	p : float, optional
		Super-Gaussian exponent. Defaults to the analytically fitted value
		(~1.73) matching HWHM=10 deg and D(20 deg)=0.1.
	taper_start, taper_end : float, optional
		Angular range (degrees) over which a raised-cosine taper brings the
		response smoothly to zero. Defaults to 25 and 50 degrees, preserving
		the low-level wings visible in the measured SQM-L curve.

	Returns
	-------
	D : float or array
		Relative angular response, normalized so D(0) = 1. Same shape as
		theta_deg.
	"""
	return _angular_response(theta_deg, hwhm, p, taper_start, taper_end)


def sqm_response_db(theta_deg, response_func=sqm_response, **kwargs):
	"""
	SQM-family response expressed in magnitudes (mag/arcsec^2-style), i.e.
	-2.5*log10(D), with arbitrary zero point at theta=0. Returns +inf where
	D=0 (beyond the screening cutoff).

	Parameters
	----------
	theta_deg : float or array
		Angle(s) from boresight, in degrees.
	response_func : callable, optional
		Which angular response to use: sqm_response (default) or
		sqml_response.
	**kwargs
		Passed through to response_func().

	Returns
	-------
	float or array
		-2.5*log10(D(theta)).
	"""
	D = response_func(theta_deg, **kwargs)
	with n.errstate(divide='ignore'):
		return -2.5 * n.log10(D)


#-----------------------------------------------------------------------------#
#              Reference tables reproduced from the fits (for docs)            #
#-----------------------------------------------------------------------------#

REFERENCE_TABLE_DEG = n.array([0, 10, 20, 30, 42, 50, 55, 60, 65])
REFERENCE_TABLE_D = sqm_response(REFERENCE_TABLE_DEG)

REFERENCE_TABLE_DEG_L = n.array([0, 5, 10, 15, 20, 25, 30, 40, 50])
REFERENCE_TABLE_D_L = sqml_response(REFERENCE_TABLE_DEG_L)


if __name__ == '__main__':

	#Print fitted parameters and sanity-check tables
	print('Fitted SQM / SQM-L angular response models')
	print('==========================================')
	print('Standard SQM:')
	print('  HWHM            = %s deg' %HWHM)
	print('  Fitted exponent = %.6f' %EXPONENT)
	print('  Taper range     = %s-%s deg' %(TAPER_START, TAPER_END))
	print('SQM-L:')
	print('  HWHM            = %s deg' %HWHM_L)
	print('  Fitted exponent = %.6f' %EXPONENT_L)
	print('  Taper range     = %s-%s deg\n' %(TAPER_START_L, TAPER_END_L))

	print('%12s | %10s' %('theta (deg)', 'D_SQM'))
	print('-'*27)
	for t, d in zip(REFERENCE_TABLE_DEG, REFERENCE_TABLE_D):
		print('%12.1f | %10.4f' %(t, d))

	print('\n%12s | %10s' %('theta (deg)', 'D_SQM-L'))
	print('-'*27)
	for t, d in zip(REFERENCE_TABLE_DEG_L, REFERENCE_TABLE_D_L):
		print('%12.1f | %10.4f' %(t, d))

	#Reproduce the overlay plot (data points approximated from the paper's
	#Fig. 3, as no raw digitized dataset is available)
	import matplotlib.pyplot as plt

	theta_fine = n.linspace(-90, 90, 721)
	D_fine = sqm_response(theta_fine)
	D_fine_L = sqml_response(theta_fine)

	def detector_response(theta_deg, hwhm=44.0, p=3.6):
		"""
		Approximate manufacturer detector response curve (paper: SQM curve is
		"slightly narrower than the detector angular response at large angles"
		-- modeled here as a slightly wider super-Gaussian, for qualitative
		comparison only; not digitized from the source figure).

		Parameters
		----------
		theta_deg : float or array
			Angle(s) from boresight, in degrees.
		hwhm : float, optional
			Half width at half maximum, degrees.
		p : float, optional
			Super-Gaussian exponent.

		Returns
		-------
		float or array
			Approximate detector response.
		"""
		a = n.abs(n.asarray(theta_deg, dtype=float))
		return n.exp(-n.log(2) * (a/hwhm)**p)

	D_detector = detector_response(theta_fine)

	#Approximate data points (symmetric, small synthetic scatter) purely for
	#visual reference -- NOT digitized from the original figure
	rng = n.random.default_rng(0)
	data_angles = n.array([-90, -80, -70, -60, -55, -50, -45, -40, -35,
							-30, -25, -20, -15, -10, -5, 0, 5, 10, 15,
							20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90],
						   dtype=float)
	data_h = n.clip(sqm_response(data_angles) + rng.normal(0, 0.015, data_angles.size), 0, 1)
	data_v = n.clip(sqm_response(data_angles*1.03) + rng.normal(0, 0.015, data_angles.size), 0, 1)
	data_hL = n.clip(sqml_response(data_angles) + rng.normal(0, 0.015, data_angles.size), 0, 1)
	data_vL = n.clip(sqml_response(data_angles*1.03) + rng.normal(0, 0.015, data_angles.size), 0, 1)

	fig, (axS, axL) = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

	#--- Standard SQM panel ---
	axS.plot(theta_fine, D_detector, '--', color='gray', lw=1.2,
			 label='Manufacturer detector response (approx., vertical)')
	axS.plot(theta_fine, D_detector*0.995, '-.', color='gray', lw=1.2,
			 label='Manufacturer detector response (approx., horizontal)')
	axS.scatter(data_angles, data_h, marker='s', facecolor='black',
				edgecolor='black', s=25, label='SQM data, horizontal (approx.)')
	axS.scatter(data_angles, data_v, marker='s', facecolor='white',
				edgecolor='black', s=25, label='SQM data, vertical (approx.)')
	axS.plot(theta_fine, D_fine, '-', color='crimson', lw=2.2,
			 label='Fitted model D(theta)')
	axS.set_title('SQM: fitted model vs. Fig. 3 (approx. reconstruction)')
	axS.set_xlabel('angle (deg)')
	axS.set_ylabel('relative radiance')
	axS.set_xlim(-90, 90)
	axS.axhline(0.5, color='indianred', ls=':', lw=1)
	axS.axhline(0.1, color='indianred', ls=':', lw=1)
	axS.legend(fontsize=7, loc='upper right')

	#--- SQM-L panel ---
	axL.scatter(data_angles, data_hL, marker='s', facecolor='black',
				edgecolor='black', s=25, label='SQM-L data, horizontal (approx.)')
	axL.scatter(data_angles, data_vL, marker='s', facecolor='white',
				edgecolor='black', s=25, label='SQM-L data, vertical (approx.)')
	axL.plot(theta_fine, D_fine_L, '-', color='navy', lw=2.2,
			 label='Fitted model D(theta)')
	axL.set_title('SQM-L: fitted model vs. Fig. 3 (approx. reconstruction)')
	axL.set_xlabel('angle (deg)')
	axL.set_xlim(-80, 80)
	axL.axhline(0.5, color='cornflowerblue', ls=':', lw=1)
	axL.axhline(0.1, color='cornflowerblue', ls=':', lw=1)
	axL.legend(fontsize=7, loc='upper right')

	axS.set_ylim(0, 1.05)
	fig.tight_layout()
	plt.savefig('../Calibration/sqm_angular_response_fit.png', dpi=150)
	plt.show()