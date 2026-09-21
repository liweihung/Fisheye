#-----------------------------------------------------------------------------#
#star_visibility_model.py
#
#NPS Night Skies Program
#
#Self-contained model for percentage of naked-eye-visible stars as a
#function of All-Sky Light Pollution Ratio (ALR). Import predict_visibility()
#in any other script to evaluate the model; no external data file needed.
#
#Model form (fit in log10(x)/log10(y) space):
#   log10(y) = log10(d) + (log10(a)-log10(d)) / (1 + 10**(b*(log10(x)-log10(c))))
#where:
#   x = ALR (dimensionless, must be > 0)
#   y = % of naked-eye-visible stars (0-100)
#   a = upper asymptote, % visible as x -> 0 (fixed at 100, physical ceiling)
#   d = lower asymptote, % visible as x -> infinity
#   c = inflection point, in ALR units
#   b = steepness of the transition
#
#Fit to 649 field observations from the NPS Natural Sounds & Night Skies
#Division dataset (NPMaps.xlsx, columns ALR_POS and VISSTARS_PCT), via
#nonlinear least squares (scipy.optimize.curve_fit). RMSE in log10(y)
#space ~= 0.045. Fit only over ALR in [0.01, 66.83]; predictions outside
#this range are unvalidated extrapolations, and predict_visibility() warns
#when called there.
#
#Usage:
#	from star_visibility_model import predict_visibility
#	predict_visibility(0.5)          #-> ~89.1 (single value)
#	predict_visibility([0.1, 1, 10]) #-> array of values (vectorized)
#
#History:
#	Li-Wei Hung -- Created
#
#-----------------------------------------------------------------------------#
import numpy as n

#-----------------------------------------------------------------------------#
#                          Fitted model parameters                            #
#-----------------------------------------------------------------------------#

MODEL_PARAMS = {
	'a': 100.0,   #upper asymptote (%), fixed by physical constraint
	'd': 4.9688,  #lower asymptote (%)
	'c': 18.7140, #inflection point (ALR units)
	'b': 0.8366,  #steepness
}

VALID_RANGE = {
	'x_min_observed': 0.01,
	'x_max_observed': 66.83,
}

MODEL_METADATA = {
	'model_name': 'star_visibility_vs_alr_logistic',
	'description': ('Percentage of naked-eye stars visible as a function of '
					 'All-Sky Light Pollution Ratio (ALR_POS), fit to 649 '
					 'NPS Natural Sounds & Night Skies Division field '
					 'observations.'),
	'model_form': ('log10(y) = log10(d) + (log10(a)-log10(d)) / '
				   '(1 + 10**(b*(log10(x)-log10(c))))'),
	'variables': {
		'x': 'ALR_POS -- All-Sky Light Pollution Ratio (dimensionless, x > 0)',
		'y': 'VISSTARS_PCT -- percentage of naked-eye-visible stars (0-100)',
	},
	'fit_details': {
		'method': ('nonlinear least squares (scipy.optimize.curve_fit), '
				   'fit in log10(x)/log10(y) space'),
		'n_observations': 649,
		'rmse_log10_y_space': 0.045,
		'upper_asymptote_constraint': ('fixed at exactly 100.0 (physical '
										'ceiling: visibility cannot exceed '
										'all stars being seen); unconstrained '
										'fit had landed at ~100.2, a '
										'physically meaningless overshoot'),
	},
	'source_data': ('NPMaps.xlsx (NPS Natural Sounds & Night Skies Division), '
					 'columns ALR_POS and VISSTARS_PCT'),
	'date_fit': '2026-09-18',
}

#-----------------------------------------------------------------------------#

def predict_visibility(x, warn_outside_range=True):
	"""
	Predict the percentage of naked-eye-visible stars for a given All-Sky
	Light Pollution Ratio (ALR).

	Parameters
	----------
	x : float, list, or array
		ALR value(s). Must be > 0.
	warn_outside_range : bool, optional
		If True (default), prints a warning when any input falls outside
		the range the model was fit on ([0.01, 66.83]). Predictions
		outside this range are extrapolations, not validated against
		real data.

	Returns
	-------
	y : float or array
		Predicted percentage of stars visible (0-100), same shape as
		input (a plain float is returned for scalar input).
	"""
	x_arr = n.asarray(x, dtype=float)

	if n.any(x_arr <= 0):
		raise ValueError('ALR values must be strictly positive (the model '
						  'works in log10(x) space and is undefined for x <= 0).')

	if warn_outside_range:
		lo, hi = VALID_RANGE['x_min_observed'], VALID_RANGE['x_max_observed']
		out_of_range = (x_arr < lo) | (x_arr > hi)
		if n.any(out_of_range):
			print('Warning: some input value(s) fall outside the observed '
				  'fit range [%s, %s]. Predictions there are extrapolations, '
				  'not validated against real data.' %(lo, hi))

	a, d, c, b = MODEL_PARAMS['a'], MODEL_PARAMS['d'], MODEL_PARAMS['c'], MODEL_PARAMS['b']
	log_x = n.log10(x_arr)
	log_a, log_d, log_c = n.log10(a), n.log10(d), n.log10(c)
	log_y = log_d + (log_a-log_d) / (1 + 10**(b*(log_x-log_c)))
	y = 10**log_y

	if y.shape == ():
		return float(y)
	return y


def model_info():
	"""
	Return the full model metadata dict (parameters, fit details, etc).

	Returns
	-------
	dict
		MODEL_METADATA merged with 'parameters' (MODEL_PARAMS) and
		'valid_range' (VALID_RANGE).
	"""
	return {**MODEL_METADATA, 'parameters': MODEL_PARAMS, 'valid_range': VALID_RANGE}


if __name__ == '__main__':

	#Quick self-check: reproduce known values from the original fit
	test_points = [0.01, 0.1, 1, 5, 10, 30, 70]
	print('ALR (x)   ->  % visible (y)')
	for xv in test_points:
		print('%7s   ->  %.1f' %(xv, predict_visibility(xv, warn_outside_range=False)))

	#Prompt the user for an ALR value and print the predicted visibility
	print()
	user_input = input('Enter an ALR value to predict % of visible stars: ')
	try:
		x_user = float(user_input)
		y_user = predict_visibility(x_user)
		print('ALR = %s  ->  %.1f%% of naked-eye stars visible' %(x_user, y_user))
	except ValueError as e:
		print('Invalid input:', e)