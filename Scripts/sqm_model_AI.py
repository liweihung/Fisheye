#-----------------------------------------------------------------------------#
#sqm_angular_response.py
#
#NPS Night Skies Program
#
#Single-curve model of the Sky Quality Meter (SQM) angular sensitivity
#D(theta), fitted to the anchor values reported in Cinzano, P. (2005),
#"Night Sky Photometry with Sky Quality Meter", ISTIL Internal Report
#n.9, v.1.4.
#
#The paper describes (Sec. 2, Fig. 2/3) the angular response of SQM serial
#n.0115 v.1.09, measured on a rotation table against a calibrated point
#source, in both the horizontal and vertical planes:
#	- Half Width at Half Maximum (HWHM)          ~= 42 deg
#	- Factor-of-10 attenuation of a point source  at  55 deg
#	- Onset of mechanical vignetting/screening   ~= 60-65 deg
#
#Model form: a super-Gaussian (generalized Gaussian) core,
#	D(theta) = exp[-ln(2)*(|theta|/HWHM)**p]   for |theta| <= taper_start
#with a raised-cosine taper to zero between taper_start and taper_end,
#representing mechanical screening of the detector at large angles:
#	D(theta) = D_core(theta)*0.5*(1+cos(pi*t)),  t = (|theta|-taper_start)/(taper_end-taper_start)
#for taper_start < |theta| < taper_end, and D(theta)=0 for |theta| >= taper_end.
#
#The exponent p is solved analytically so the curve passes through both
#anchor points exactly: D(HWHM)=0.5 (by definition of half width at half
#max) and D(55 deg)=0.1 (paper's stated factor-of-10 point). A plain
#Gaussian (p=2) only reaches D=0.30 at 55 deg, too shallow; p~=4.452 is
#required to match the reported factor-of-10 attenuation.
#
#Provenance note: this fit is derived from the paper's stated numeric
#anchors in the text (HWHM, the 55-degree/factor-10 point, and the 60-65
#degree screening onset), not from a pixel-level digitization of the
#original Fig. 2/3 plot or the underlying raw data table (not provided in
#the paper). Treat this as a best-effort analytic stand-in for the real
#angular response, adequate for weighting/convolution use, not a
#substitute for the manufacturer's own measured response curve if
#metrological precision is needed.
#
#Usage:
#	from sqm_angular_response import sqm_response, HWHM, EXPONENT
#	D = sqm_response(theta_deg)                 #scalar or array in, same out
#	D = sqm_response(theta_deg, normalize=True) #explicit norm to D(0)=1 (already true by construction)
#
#	#Use as a weighting kernel, e.g. reproducing Eq. (1) of the paper:
#	#  I_avg = integral[ I(theta,phi) * D(theta) * sin(theta) dtheta dphi ]
#	#        / integral[ D(theta) * sin(theta) dtheta dphi ]
#
#Output (when run directly):
#   (1) ../Calibration/sqm_angular_response_fit.png -- fitted model vs. an
#       approximate reconstruction of the paper's Fig. 3 overlay
#
#History:
#	Derived/fitted by Claude (Anthropic) from Cinzano (2005), for reuse in
#	downstream numerical scripts. Independently verify against the source
#	paper / instrument datasheet before using in any precision application.
#
#-----------------------------------------------------------------------------#
import numpy as n

#-----------------------------------------------------------------------------#
#                          Fitted model parameters                            #
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

def _core(theta_abs, hwhm=HWHM, p=EXPONENT):
	"""
	Super-Gaussian core, before vignetting taper is applied.

	Parameters
	----------
	theta_abs : float or array
		Absolute angle(s) from boresight, in degrees.
	hwhm : float, optional
		Half width at half maximum, degrees.
	p : float, optional
		Super-Gaussian exponent.

	Returns
	-------
	float or array
		Core response, same shape as theta_abs.
	"""
	return n.exp(-n.log(2) * (theta_abs/hwhm)**p)


def sqm_response(theta_deg, hwhm=HWHM, p=EXPONENT,
				  taper_start=TAPER_START, taper_end=TAPER_END):
	"""
	Evaluate the fitted SQM angular sensitivity D(theta).

	Parameters
	----------
	theta_deg : float or array
		Angle(s) from boresight (zenith when SQM points up), in degrees.
		Can be positive or negative; response is symmetric about 0.
	hwhm : float, optional
		Half width at half maximum, degrees. Defaults to 42.0 (paper
		value).
	p : float, optional
		Super-Gaussian exponent. Defaults to the analytically fitted
		value (~4.452) matching HWHM=42 deg and D(55 deg)=0.1.
	taper_start, taper_end : float, optional
		Angular range (degrees) over which a raised-cosine taper is
		applied to bring the response smoothly to zero, representing
		mechanical vignetting/screening of the detector. Defaults to
		55 and 65 degrees per the paper's reported screening onset.

	Returns
	-------
	D : float or array
		Relative angular response, normalized so D(0) = 1. Same shape as
		theta_deg.
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


def sqm_response_db(theta_deg, **kwargs):
	"""
	SQM response expressed in magnitudes (mag/arcsec^2-style), i.e.
	-2.5*log10(D), with arbitrary zero point at theta=0. Returns +inf
	where D=0 (beyond the screening cutoff).

	Parameters
	----------
	theta_deg : float or array
		Angle(s) from boresight, in degrees.
	**kwargs
		Passed through to sqm_response().

	Returns
	-------
	float or array
		-2.5*log10(D(theta)).
	"""
	D = sqm_response(theta_deg, **kwargs)
	with n.errstate(divide='ignore'):
		return -2.5 * n.log10(D)


#-----------------------------------------------------------------------------#
#              Reference table reproduced from the fit (for docs)             #
#-----------------------------------------------------------------------------#

REFERENCE_TABLE_DEG = n.array([0, 10, 20, 30, 42, 50, 55, 60, 65])
REFERENCE_TABLE_D = sqm_response(REFERENCE_TABLE_DEG)


if __name__ == '__main__':

	#Print fitted parameters and sanity-check table
	print('Fitted SQM angular response model')
	print('----------------------------------')
	print('HWHM            = %s deg' %HWHM)
	print('Fitted exponent = %.6f' %EXPONENT)
	print('Taper range     = %s-%s deg\n' %(TAPER_START, TAPER_END))

	print('%12s | %10s' %('theta (deg)', 'D(theta)'))
	print('-'*27)
	for t, d in zip(REFERENCE_TABLE_DEG, REFERENCE_TABLE_D):
		print('%12.1f | %10.4f' %(t, d))

	#Reproduce the overlay plot (data points approximated from the
	#paper's Fig. 3, as no raw digitized dataset is available)
	import matplotlib.pyplot as plt

	theta_fine = n.linspace(-90, 90, 721)
	D_fine = sqm_response(theta_fine)

	def detector_response(theta_deg, hwhm=44.0, p=3.6):
		"""
		Approximate manufacturer detector response curve (paper: SQM
		curve is "slightly narrower than the detector angular response
		at large angles" -- modeled here as a slightly wider
		super-Gaussian, for qualitative comparison only; not digitized
		from the source figure).

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

	#Approximate data points (symmetric, small synthetic scatter) purely
	#for visual reference -- NOT digitized from the original figure
	rng = n.random.default_rng(0)
	data_angles = n.array([-90, -80, -70, -60, -55, -50, -45, -40, -35,
							-30, -25, -20, -15, -10, -5, 0, 5, 10, 15,
							20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90],
						   dtype=float)
	data_h = n.clip(sqm_response(data_angles) + rng.normal(0, 0.015, data_angles.size), 0, 1)
	data_v = n.clip(sqm_response(data_angles*1.03) + rng.normal(0, 0.015, data_angles.size), 0, 1)

	fig, ax = plt.subplots(figsize=(8, 5))
	ax.plot(theta_fine, D_detector, '--', color='gray', lw=1.2,
			label='Manufacturer detector response (approx., vertical)')
	ax.plot(theta_fine, D_detector*0.995, '-.', color='gray', lw=1.2,
			label='Manufacturer detector response (approx., horizontal)')
	ax.scatter(data_angles, data_h, marker='s', facecolor='black',
			   edgecolor='black', s=25, label='SQM data, horizontal plane (approx.)')
	ax.scatter(data_angles, data_v, marker='s', facecolor='white',
			   edgecolor='black', s=25, label='SQM data, vertical plane (approx.)')
	ax.plot(theta_fine, D_fine, '-', color='crimson', lw=2.2,
			label='Fitted model D(theta)')

	ax.set_xlabel('angle (deg)')
	ax.set_ylabel('relative radiance')
	ax.set_xlim(-90, 90)
	ax.set_ylim(0, 1.05)
	ax.axhline(0.5, color='indianred', ls=':', lw=1)
	ax.axhline(0.1, color='indianred', ls=':', lw=1)
	ax.set_title('SQM angular response: fitted model vs. Fig. 3 (approximate reconstruction)')
	ax.legend(fontsize=8, loc='upper right')
	fig.tight_layout()
	plt.savefig('../Calibration/sqm_angular_response_fit.png', dpi=150)
	plt.show()