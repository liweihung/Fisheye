#-----------------------------------------------------------------------------#
#terrain_mask.py
#
#NPS Night Skies Program
#
#Detects the sky-terrain boundary in a fisheye image and builds a terrain
#mask. The fisheye image is first transformed to a rectangular (r, theta)
#projection, then the sky-terrain boundary is detected and refined in four
#steps:
#   (1) Find the first drop beyond a threshold in each column.
#   (2) Clip and median-filter the resulting boundary indices.
#   (3) Fine-tune each column's boundary to be near the median-filtered estimate.
#   (4) Median-filter again to smooth the final boundary.
#The refined rectangular boundary is then transformed back into the
#fisheye projection to produce the final mask.
#
#After a mask is built, the user is shown the result and asked whether it
#looks satisfactory. If not, the user can supply new boundary-detection
#parameters (t1, m1, t2, m2) and the mask is rebuilt, repeating until the
#user accepts the result.
#
#Input:
#   (1) MF_* -- median-filtered reference fisheye image
#	(2) ../Calibration/imagecenter.csv -- fisheye center coordinates and
#	    FOV radius, keyed by camera
#
#Output:
#   (1) p.data_cal+'mask.png' -- boundary detection process, grayscale
#   (2) p.data_cal+'mask_color.png' -- boundary detection process, in color
#   (3) p.data_cal+'mask.fit' -- terrain mask in fisheye view, for
#       scientific analysis (1 = sky, NaN = terrain/beyond FoV)
#   (4) p.data_cal+'mask_fisheye.fit' -- reference image with mask applied,
#       for making movies/gifs
#   (5) p.data_cal+'mask_fisheye.png' -- reference image with mask applied
#
#History:
#	Li-Wei Hung -- Created
#
#-----------------------------------------------------------------------------#
import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pandas as pd
import shutil

from astropy.io import fits
from glob import glob
from scipy.ndimage import median_filter

# Local Source
import colormaps
import process_input as p

#-----------------------------------------------------------------------------#

def fisheye_to_rectangular(fisheye_img, xc, yc, radius):
	"""
	Transform a fisheye image into a rectangular (r, theta) projection.

	Each output row y corresponds to a pixel radius r = y from the fisheye
	center, and each output column x corresponds to an azimuthal angle
	theta = x/radius (0 to 2*pi). Output pixels that fall outside the
	input image are left as NaN.

	Vectorized over the full output grid at once (rather than looping
	over every output pixel in Python), for speed on large images.

	Parameters
	----------
	fisheye_img : 2D array
		Fisheye image, shape (height, width).
	xc, yc : float
		Fisheye image center coordinates [pix].
	radius : int
		Fisheye field-of-view radius [pix]. Sets the output image size:
		height = radius+50, width = int(2*pi*radius).

	Returns
	-------
	output_img : 2D array
		Rectangular projection, shape (radius+50, int(2*pi*radius)). NaN
		wherever the corresponding fisheye coordinate falls outside
		fisheye_img.
	"""
	height, width = fisheye_img.shape

	output_height = radius+50
	output_width = int(2*np.pi*radius)

	#r corresponds to output row y; theta corresponds to output column x
	y_idx, x_idx = np.meshgrid(np.arange(output_height), np.arange(output_width), indexing='ij')
	theta = x_idx / radius
	r = y_idx

	#convert (r, theta) back to fisheye (x, y) pixel coordinates
	fisheye_x = (r*np.cos(theta) + xc).astype(int)
	fisheye_y = (r*np.sin(theta) + yc).astype(int)

	#only sample fisheye_img where the computed coordinates are in bounds
	in_bounds = (fisheye_x >= 0) & (fisheye_x < width) & (fisheye_y >= 0) & (fisheye_y < height)

	output_img = np.full((output_height, output_width), np.nan)
	output_img[in_bounds] = fisheye_img[fisheye_y[in_bounds], fisheye_x[in_bounds]]

	return output_img


def terrain_boundary_mask(rectangular_img, R, t1=0.6, m1=80, t2=0.5, m2=10):
	"""
	Detect the sky-terrain boundary in a rectangular (r, theta) image,
	using a multi-step thresholding and filtering approach:
		1. Mask out everything below the first NaN in each column.
		2. Detect the first sharp vertical drop beyond threshold t1.
		3. Refine the boundary near the median-filtered estimate, using
		   the finer threshold t2.
		4. Smooth the final boundary with a second median filter.

	Parameters
	----------
	rectangular_img : 2D array
		Rectangular (r, theta) image, NaN indicating invalid/masked
		regions. Modified in place (values below the first NaN in each
		column are set to NaN).
	R : int or float
		Fisheye field-of-view radius [pix], used to clip the detected
		boundary indices to a valid range.
	t1 : float, optional
		Threshold for initial boundary detection. Defaults to 0.6.
	m1 : int, optional
		Window size for the first median filter. Defaults to 80.
	t2 : float, optional
		Threshold for fine-tuning the boundary near the median-filtered
		estimate. Defaults to 0.5.
	m2 : int, optional
		Window size for the second median filter. Defaults to 10.

	Returns
	-------
	list of 1D arrays
		[jump_indices, jump_indices_mf, jump_indices2, jump_indices_mf2]:
		the initial, median-filtered, refined, and final smoothed
		boundary index for each column.
	"""
	#default boundary: the first NaN in each column
	nan_mask = np.isnan(rectangular_img)
	jump_indices = np.argmax(nan_mask, axis=0)
	for col in range(rectangular_img.shape[1]):
		rectangular_img[jump_indices[col]:, col] = np.nan

	#--------------------------------------------------------------------------#
	#						  Finding the sky-terrain boundary				   #
	#--------------------------------------------------------------------------#

	#Step 1: first drop beyond threshold t1, per column
	differences = np.abs(np.diff(rectangular_img, n=5, axis=0))

	for col in range(differences.shape[1]):
		jump_index = np.where(differences[:, col] > t1)[0]
		if jump_index.size > 0:
			jump_indices[col] = jump_index[0]

	#Step 2: clip and median-filter
	jump_indices = np.clip(jump_indices, 0, R)
	jump_indices_mf = median_filter(jump_indices, size=m1)

	#Step 3: fine-tune near the median-filtered estimate, using threshold t2
	jump_indices2 = jump_indices_mf.copy()

	for col in range(differences.shape[1]):
		jump_index = np.where(differences[:, col] > t2)[0]
		if jump_index.size > 0:
			a = np.min(np.abs(jump_index-jump_indices_mf[col]))
			b = np.abs(R-jump_indices_mf[col])
			if a < b:
				w = np.argmin(np.abs(jump_index-jump_indices_mf[col]))
				jump_indices2[col] = jump_index[w]

	#Step 4: final smoothing
	jump_indices_mf2 = median_filter(jump_indices2, size=m2)

	return [jump_indices, jump_indices_mf, jump_indices2, jump_indices_mf2]


def mask_fisheye_image(fisheye_img, rectangular_mask, xc, yc, radius):
	"""
	Mask a fisheye image using a rectangular (r, theta) terrain mask.

	Converts each fisheye pixel's (x, y) coordinate to polar (r, theta),
	masks out pixels beyond the field-of-view radius, then masks out
	terrain within the field of view based on rectangular_mask.

	Vectorized over the full image at once (rather than looping over every
	unmasked pixel in Python), for speed on large images.

	Parameters
	----------
	fisheye_img : 2D array
		Fisheye image to mask.
	rectangular_mask : 2D array
		Rectangular (r, theta) mask (NaN = terrain/masked, elsewhere =
		valid sky), from terrain_boundary_mask().
	xc, yc : float
		Fisheye image center coordinates [pix].
	radius : int
		Fisheye field-of-view radius [pix].

	Returns
	-------
	masked_fisheye : 2D array
		fisheye_img with terrain and beyond-FoV pixels set to NaN.
	"""
	masked_fisheye = fisheye_img.copy()
	height, width = fisheye_img.shape

	x, y = np.meshgrid(np.arange(width), np.arange(height))
	R_grid = np.sqrt((x-xc)**2 + (y-yc)**2).astype(int)
	Theta = np.arctan2(y-yc, x-xc) #-pi to pi
	Theta = (np.where(Theta < 0, Theta+2*np.pi, Theta) * radius).astype(int)

	#beyond the fisheye's FoV
	beyond_fov = R_grid > radius
	masked_fisheye[beyond_fov] = np.nan

	#terrain within the FoV: clip indices into the rectangular mask's
	#valid range before indexing, then look up validity per-pixel
	within_fov = ~beyond_fov
	r_idx = np.clip(R_grid, 0, rectangular_mask.shape[0]-1)
	theta_idx = np.clip(Theta, 0, rectangular_mask.shape[1]-1)
	is_terrain = np.isnan(rectangular_mask[r_idx, theta_idx])

	masked_fisheye[within_fov & is_terrain] = np.nan

	return masked_fisheye


def plot(image, jump_ind, color=False, show=False, save_path=None):
	"""
	Plot the original image alongside four progressively masked versions,
	one per boundary-detection stage, to visualize the terrain boundary
	detection process.

	Parameters
	----------
	image : 2D array
		Input rectangular (r, theta) image.
	jump_ind : list of 1D arrays
		Four boundary index arrays (one per column), from
		terrain_boundary_mask(); each produces one masked subplot.
	color : bool, optional
		If True, use the 'NPS_mag' colormap with masked pixels shown in
		black. If False, use grayscale. Defaults to False.
	show : bool, optional
		If True, display the plot interactively. Defaults to False.
	save_path : str, optional
		Destination image path. Defaults to p.data_cal+'mask_color.png'
		or p.data_cal+'mask.png', depending on color.
	"""
	fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(10.5, 9), sharex=True, sharey=True)
	plt.tight_layout()

	#build the original plus four progressively masked images
	masked_images = [image.copy()]
	for idx in jump_ind:
		masked = image.copy()
		masked[np.arange(image.shape[0])[:, None] >= idx] = np.nan
		masked_images.append(masked)

	if color:
		cmap = copy.copy(mpl.colormaps["NPS_mag"])
		cmap.set_bad(color='black')
		vmin, vmax = 14, 24
	else:
		cmap = 'gray'
		vmin = vmax = None

	for ax, img in zip(axs, masked_images):
		ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)

	filename = save_path or (p.data_cal+"mask_color.png" if color else p.data_cal+"mask.png")
	plt.savefig(filename, dpi=300)

	if show:
		plt.show(block=False)


def build_and_save_mask(Xc, Yc, R, fisheye_image, ref_img, t1=0.6, m1=80, t2=0.5, m2=10):
	"""
	Build the terrain mask from the reference fisheye image using the
	given boundary-detection parameters, display the detection-process
	plots, and save all mask output files.

	Parameters
	----------
	Xc, Yc : float
		Fisheye image center coordinates [pix].
	R : int or float
		Fisheye field-of-view radius [pix].
	fisheye_image : 2D array
		Median-filtered reference fisheye image.
	ref_img : str
		Path to the reference image file, used for the BASEIMG header.
	t1, m1, t2, m2 : optional
		Boundary-detection parameters, passed through to
		terrain_boundary_mask(). See that function for details.
	"""
	print('Building mask with t1=%s, m1=%s, t2=%s, m2=%s' %(t1, m1, t2, m2))

	rectangular_image = fisheye_to_rectangular(fisheye_image, Xc, Yc, R)
	jump_indices = terrain_boundary_mask(rectangular_image, R, t1=t1, m1=m1, t2=t2, m2=m2)

	#build the mask in fisheye view from the final smoothed boundary
	masked = rectangular_image.copy()
	masked[np.arange(rectangular_image.shape[0])[:, None] >= jump_indices[3]] = np.nan
	masked_fisheye = mask_fisheye_image(fisheye_image, masked, Xc, Yc, R)

	mask = np.ones_like(masked_fisheye)
	mask[np.isnan(masked_fisheye)] = np.nan

	#--------------------------------------------------------------------------#
	#								Display the images						   #
	#--------------------------------------------------------------------------#

	plot(rectangular_image, jump_indices, color=False, show=False)
	plot(rectangular_image, jump_indices, color=True, show=True)

	#--------------------------------------------------------------------------#
	#								   Save the mask						   #
	#--------------------------------------------------------------------------#

	hdu = fits.PrimaryHDU()
	hdu.header['XCENTER'] = Xc
	hdu.header['YCENTER'] = Yc
	hdu.header['RADIUS'] = R
	hdu.header['BASEIMG'] = ref_img[18:]

	hdu.data = mask
	hdu.writeto(p.data_cal+'mask.fit', overwrite=1) #for scientific analysis

	hdu.data = masked_fisheye
	hdu.writeto(p.data_cal+'mask_fisheye.fit', overwrite=1) #for making a gif

	plt.figure()
	plt.imshow(masked_fisheye, cmap='gray')
	plt.savefig(p.data_cal+'mask_fisheye.png', dpi=300)
	plt.show(block=False)


def prompt_for_new_parameters(t1, m1, t2, m2):
	"""
	Explain what each boundary-detection parameter controls, print the
	values just used, and ask the user for new values to try.

	Parameters
	----------
	t1, m1, t2, m2 : float or int
		The parameter values just used, shown to the user as defaults.

	Returns
	-------
	t1, m1, t2, m2 : float or int
		The new parameter values (unchanged from the input if the user
		presses Enter without typing a value).
	"""
	print('\nParameters used: t1=%s, m1=%s, t2=%s, m2=%s' %(t1, m1, t2, m2))
	print('  t1 -- initial detection threshold (Step 1). Lower it if the '
		  'boundary is missed entirely (terrain is not masked out); raise '
          'it if noise/stars are mistaken for the boundary.')
	print('  m1 -- smoothing window before fine-tuning (Step 2). Increase '
		  'it if the boundary is jagged or jumps around; decrease it if '
		  'real terrain features (peaks, ridgelines) are being smoothed '
		  'away.')
	print('  t2 -- fine-tuning threshold (Step 3). Lower it if it missed '
          'masking out the terrian. Raise it if noise are mistaken for '
          'the boundary and the mask is too aggressive.')
	print('  m2 -- final smoothing window (Step 4). Increase it for a '
		  'smoother final boundary; decrease it to preserve sharper '
		  'terrain detail.\n')

	def ask(name, current, cast):
		raw = input('New %s [%s]: ' %(name, current)).strip()
		if raw == '':
			return current
		try:
			return cast(raw)
		except ValueError:
			print('  Could not parse "%s", keeping %s.' %(raw, current))
			return current

	t1 = ask('t1', t1, float)
	m1 = ask('m1', m1, int)
	t2 = ask('t2', t2, float)
	m2 = ask('m2', m2, int)

	return t1, m1, t2, m2


def main():
	"""
	Detects the sky-terrain boundary for this dataset's reference image
	and builds/saves the terrain mask. See the script description for
	detail. All inputs are defined and passed through process_input.

	If p.data_cal+'mask.fit' already exists, the user is asked whether to
	regenerate it; answering no skips processing entirely. After building
	a mask, the user is shown the result and asked whether it looks good;
	if not, they can supply new t1/m1/t2/m2 values and the mask is
	rebuilt, repeating until the user is satisfied.
	"""

	#--------------------------------------------------------------------------#
	#					  Check for an existing mask							   #
	#--------------------------------------------------------------------------#

	mask_path = p.data_cal+'mask.fit'
	if os.path.exists(mask_path):
		response = input(
			'The terrain mask for %s already exists. Regenerate it? [y/N]: ' %p.datanight).strip().lower()
		if response not in ('y', 'yes'):
			print('Keeping existing mask. Skipping.')
			return

	print('Generating terrain mask for %s' % p.datanight)

	#--------------------------------------------------------------------------#
	#						 Read in files and image preparation			   #
	#--------------------------------------------------------------------------#

	#image center coordinates and FOV radius
	C = pd.read_csv('../Calibration/imagecenter.csv', index_col=0)
	Xc = C['Xcenter'][p.camera]
	Yc = C['Ycenter'][p.camera]
	R = C['Radius'][p.camera]

	#median-filtered reference fisheye image
	ref_img = glob(p.data_cal+'MF_'+p.reference)[0]
	with fits.open(ref_img, uint=False, memmap=False) as hdul:
		fisheye_image = hdul[0].data.astype(float, copy=True)

	#--------------------------------------------------------------------------#
	#				 Build the mask, review, and refine if needed			   #
	#--------------------------------------------------------------------------#

	t1, m1, t2, m2 = 0.6, 80, 0.5, 10 #boundary-detection defaults

	while True:
		build_and_save_mask(Xc, Yc, R, fisheye_image, ref_img, t1=t1, m1=m1, t2=t2, m2=m2)

		response = input('\nIs this mask satisfactory? [Y/n]: ').strip().lower()
		if response in ('', 'y', 'yes'):
			print('Mask accepted.')
			break

		t1, m1, t2, m2 = prompt_for_new_parameters(t1, m1, t2, m2)


if __name__ == '__main__':

	#copy this dataset's process_input.py so it can be imported as p
	f = glob('../Data_processed/ANDE_20250915/process_input.py')
	shutil.copy(f[0], 'process_input_mask.py')
	import process_input_mask as p

	main()