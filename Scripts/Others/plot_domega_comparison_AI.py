#-----------------------------------------------------------------------------#
#plot_domega_comparison.py
#
#NPS Night Skies Program
#
#Plots the per-pixel solid angle dOmega, and separately cos(theta)*dOmega
#(the horizontal illuminance weighting), as functions of pixel radius r,
#comparing the empirically calibrated theta(r) relationship (theta_r.xlsx)
#against the idealized equidistant fisheye assumption
#(theta = pi/2 * r/R, dOmega = pi/(2*R*r) * sin(pi*r/(2*R))).
#
#Input:
#   (1) theta_r.xlsx -- empirically measured theta(r) calibration table
#
#Output:
#   (1) domega_comparison.png -- dOmega comparison plot
#   (2) cos_theta_domega_comparison.png -- cos(theta)*dOmega comparison plot
#
#-----------------------------------------------------------------------------#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------#

#Path to the calibration file -- adjust as needed
THETA_R_PATH = '../../Calibration/theta_r.xlsx'

#Output image paths
OUTPUT_PATH_DOMEGA = '../../Performance/Metrics_calculation/domega_comparison.png'
OUTPUT_PATH_COS_DOMEGA = '../../Performance/Metrics_calculation/cos_theta_domega_comparison.png'

#-----------------------------------------------------------------------------#

def load_theta_r_calibration(path):
	"""
	Load the theta(r) calibration table and compute per-band dtheta/dr,
	theta, and solid angle density, identical to metrics.py's calibration
	logic.

	Parameters
	----------
	path : str
		Path to the theta_r.xlsx calibration file.

	Returns
	-------
	r_band_mid : 1D array
		Band midpoint pixel radii.
	theta_rad : 1D array
		Zenith angle [radians] at each band (same bands as r_band_mid).
	dOmega_at_bands : 1D array
		Per-pixel solid angle density [sr] at each band midpoint.
	R : float
		Maximum calibrated pixel radius (at theta=90 degrees).
	"""
	cal = pd.read_excel(path)

	theta_deg = cal['Zenith_Angle (deg)'].values
	r_cumulative = cal['Cumulative_Pixel'].values
	delta_pixel = cal['Delta_Pixel'].values
	step_deg = theta_deg[0]

	dtheta_dr_per_band = np.deg2rad(step_deg) / delta_pixel
	r_prev = np.concatenate(([0], r_cumulative[:-1]))
	r_band_mid = (r_prev + r_cumulative) / 2
	theta_rad = np.deg2rad(theta_deg)

	dOmega_at_bands = np.where(
		r_band_mid == 0,
		dtheta_dr_per_band**2,
		dtheta_dr_per_band * np.sin(theta_rad) / np.where(r_band_mid == 0, 1, r_band_mid)
	)

	R = r_cumulative[-1]

	return r_band_mid, theta_rad, dOmega_at_bands, R


def dOmega_equidistant(r, R):
	"""
	Per-pixel solid angle under the idealized equidistant fisheye
	assumption: theta = pi/2 * r/R, constant dtheta/dr.

	Parameters
	----------
	r : array-like
		Pixel radius values.
	R : float
		Fisheye field-of-view radius in pixels (theta=90 degrees).

	Returns
	-------
	array
		Per-pixel solid angle [sr] at each r.
	"""
	dtheta_dr = np.pi / (2*R)
	return np.where(r == 0, dtheta_dr**2, (np.pi/(2*R*r)) * np.sin(np.pi*r/(2*R)))


def theta_equidistant(r, R):
	"""
	Zenith angle under the idealized equidistant fisheye assumption.

	Parameters
	----------
	r : array-like
		Pixel radius values.
	R : float
		Fisheye field-of-view radius in pixels (theta=90 degrees).

	Returns
	-------
	array
		Zenith angle [radians] at each r.
	"""
	return (np.pi/2) * (r/R)


def plot_comparison(r_fine, curve1, curve2, ylabel, title, outpath):
	"""
	Plot two curves (calibrated vs. equidistant) against pixel radius and
	save the result as a static image.

	Parameters
	----------
	r_fine : array
		Pixel radius values, x-axis.
	curve1 : array
		Calibrated-lens curve values, same length as r_fine.
	curve2 : array
		Equidistant-assumption curve values, same length as r_fine.
	ylabel : str
		Y-axis label.
	title : str
		Plot title.
	outpath : str
		Destination image path.
	"""
	R = r_fine[-1]

	fig, ax = plt.subplots(figsize=(9, 6))
	ax.plot(r_fine, curve1, color='#2ca89a', linewidth=2.5,
			label='Calibrated (theta_r.xlsx)')
	ax.plot(r_fine, curve2, color='#e67e22', linewidth=2.5,
			linestyle='--', label='Equidistant assumption')

	ax.axvline(0, color='gray', linestyle=':', linewidth=1)
	ax.axvline(R, color='gray', linestyle=':', linewidth=1)
	ymax = ax.get_ylim()[1]
	ax.text(5, ymax*0.97, 'zenith', color='#f39c12', fontsize=9, va='top')
	ax.text(R-90, ymax*0.97, 'horizon (θ=90°)', color='#e74c3c', fontsize=9, va='top')

	ax.set_xlabel('Pixel radius r (pixels)', fontsize=12)
	ax.set_ylabel(ylabel, fontsize=12)
	ax.set_title(title, fontsize=13)
	ax.legend(fontsize=11)
	ax.grid(alpha=0.3)

	plt.tight_layout()
	plt.savefig(outpath, dpi=200)
	print('Saved plot to', outpath)


def main():
	"""
	Loads the theta(r) calibration, computes dOmega(r) and
	cos(theta)*dOmega(r) for both the calibrated lens and the equidistant
	assumption, and plots each pair together for comparison.
	"""
	r_band_mid, theta_band, dOmega_at_bands, R = load_theta_r_calibration(THETA_R_PATH)

	r_fine = np.linspace(0.001, R, 500)

	#dOmega comparison
	dOmega_calibrated = np.interp(r_fine, r_band_mid, dOmega_at_bands)
	dOmega_equi = dOmega_equidistant(r_fine, R)

	plot_comparison(
		r_fine, dOmega_calibrated*1e6, dOmega_equi*1e6,
		ylabel='dΩ (μsr per pixel)',
		title='dΩ vs. Pixel Radius: Calibrated Lens vs. Equidistant Assumption',
		outpath=OUTPUT_PATH_DOMEGA)

	#cos(theta)*dOmega comparison (horizontal illuminance weighting)
	theta_calibrated = np.interp(r_fine, r_band_mid, theta_band)
	cos_dOmega_calibrated = np.cos(theta_calibrated) * dOmega_calibrated

	theta_equi = theta_equidistant(r_fine, R)
	cos_dOmega_equi = np.cos(theta_equi) * dOmega_equi

	plot_comparison(
		r_fine, cos_dOmega_calibrated*1e6, cos_dOmega_equi*1e6,
		ylabel='cos(θ)·dΩ (μsr per pixel)',
		title='cos(θ)·dΩ vs. Pixel Radius: Calibrated Lens vs. Equidistant Assumption',
		outpath=OUTPUT_PATH_COS_DOMEGA)

	plt.show()


if __name__ == '__main__':
	main()