#-----------------------------------------------------------------------------#
#terrain_mask.py
#
#NPS Night Skies Program
#
#Last updated: 2025/09/05
#
#This script detects the sky-terrain boundary and creats a terrain mask. The 
#script first transforms a fisheye image to a rectangular image. Then it applies 
#the following four steps to detect and refine the sky terrain boundary: 
#(1) Find the location of the first drop beyond the threshold in each column. 
#(2) Clip and median filter the jump_indices
#(3) Fine tune by finding the drop near the median filtered location
#(4) Final median filter to smooth the sky-terrain boundary 
#The rectangular mask is then transformed into the fisheye projection.
#
#Input: 
#   (1) median-filtered reference fisheye image
#	(2) imagecenter.csv
#
#Output:
#   (1) mask.png -- black and white image showing the boundary detection process
#   (2) mask_color.png -- showing the boundary detection process in color
#   (3) mask.fit -- mask in fisheye view for scientific analysis
#   (4) mask_fisheye.fit -- reference fisheye image with mask for making movies
#   (5) mask_fisheye.png -- reference fisheye image with mask 
#
#History:
#	Li-Wei Hung -- Created 
#
#-----------------------------------------------------------------------------#
import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from glob import glob
from scipy.ndimage import median_filter

# Local Source
import colormaps
import process_input as p

#-----------------------------------------------------------------------------#
def fisheye_to_rectangular(fisheye_img, xc, yc, radius):
    """
    Convert a fisheye image to a rectangular projection.

    This function transforms a fisheye image into a rectangular representation 
    based on the specified center coordinates and radius. It maps pixels 
    from the fisheye image to a new output image using polar coordinates, 
    where each point (y, x) in the output image is derived from its 
    corresponding polar coordinates (r, theta).

    :param fisheye_img: A 2D numpy array representing the fisheye image.
                        The shape of the array should be (height, width).
    :param xc: The x-coordinate of the center of the fisheye image (typically 
               the horizontal center).
    :param yc: The y-coordinate of the center of the fisheye image (typically 
               the vertical center).
    :param radius: The radius for the output rectangular mapping. It determines 
                   the height of the output image and the scale of the mapping.

    :return: A 2D numpy array representing the rectangular projection of the
             fisheye image. The shape of the output image is (radius + 50, 
             int(2 * np.pi * radius)).

    :rtype: numpy.ndarray
    """
    
    height, width = fisheye_img.shape
    theta_max = 2*np.pi  # Full range for theta
        
    output_height = radius+50
    output_width = int(2*np.pi*radius)  # Adjust width

    output_img = np.full((output_height, output_width),np.nan)

    # Loop through each pixel in the output image
    for y in range(output_height):
        for x in range(output_width):
            theta = x / radius  # Scale x to theta (0 to 2pi)
            r = y  # y corresponds to r
            
            # Convert polar coordinates (r, theta) back to fisheye coordinates
            fisheye_x = int((r * np.cos(theta)) + xc)
            fisheye_y = int((r * np.sin(theta)) + yc)
            
            if (0 <= fisheye_x < width) and (0 <= fisheye_y < height):
                pixel_value = fisheye_img[fisheye_y, fisheye_x]
                output_img[y, x] = pixel_value

    return output_img


def terrain_boundary_mask(rectangular_img, t1=0.6, t2=0.5):
    """
    Detects the sky-terrain boundary in a 2D rectangular image using a multi-step 
    thresholding and filtering approach.

    This function identifies the vertical boundary between sky and terrain by:
    1. Masking out regions below the first NaN in each column.
    2. Detecting sharp vertical drops in pixel values using a first threshold (t1).
    3. Refining the boundary using a median filter and a second, lower threshold (t2).
    4. Smoothing the final boundary with another median filter.

    Args:
        rectangular_img (np.ndarray): 2D array representing the image, where NaNs 
                                      indicate invalid or masked regions.
        t1 (float, optional): Threshold for initial boundary detection using 
                              higher sensitivity. Defaults to 0.6.
        t2 (float, optional): Threshold for fine-tuning the boundary near the 
                              median-filtered estimate. Defaults to 0.5.

    Returns:
        list of np.ndarray: A list containing four 1D arrays (one per column):
            - `jump_indices`: Initial boundary estimate based on NaNs and t1.
            - `jump_indices_mf`: Median-filtered version of `jump_indices`.
            - `jump_indices2`: Refined boundary using t2 near the median-filtered estimate.
            - `jump_indices_mf2`: Final smoothed boundary after second median filter.
    """

    # Set the default drop indices to be the first nan in each column 
    nan_mask = np.isnan(rectangular_img)
    jump_indices = np.argmax(nan_mask, axis=0)
    for col in range(rectangular_img.shape[1]):
        rectangular_img[jump_indices[col]:,col] = np.nan
        
    #----------------------------------------------------------------------#
    #   			   Finding the sky-terrain boundary                    #
    #----------------------------------------------------------------------#
    
    # Step 1: Find the first drop beyond the threshold in each column
    differences = np.abs(np.diff(rectangular_img, n=5, axis=0))
    threshold1 = t1
    
    for col in range(differences.shape[1]):
        
        jump_index = np.where(differences[:, col] > threshold1)[0]
        
        if jump_index.size > 0:
            jump_indices[col] = jump_index[0]
            
    # Step 2: Clip and median filter the jump_indices
    jump_indices = np.clip(jump_indices, 0, R)
    jump_indices_mf = median_filter(jump_indices, size=80)
    
    # Step 3: Fine tune by finding the drop near the median filtered location
    threshold2 = t2
    jump_indices2 = jump_indices_mf.copy()

    for col in range(differences.shape[1]):
        
        jump_index = np.where(differences[:, col] > threshold2)[0]

        if jump_index.size > 0:
            a = np.min(np.abs(jump_index-jump_indices_mf[col]))
            b = np.abs(R-jump_indices_mf[col])

            if a < b:
                w = np.argmin(np.abs(jump_index-jump_indices_mf[col]))
                jump_indices2[col] = jump_index[w]
    
    # Step 4: Final median filter to smooth the sky-terrain boundary 
    jump_indices_mf2 = median_filter(jump_indices2, size=10)
    
    return [jump_indices, jump_indices_mf, jump_indices2, jump_indices_mf2]
    
    
def mask_fisheye_image(fisheye_img, rectangular_mask, xc, yc, radius):
    """
    Mask a fisheye image based on a rectangular mask and the given center and 
    radius coornidates of the original fisyeye image. 

    This function first converts the xy coordinates in the fisheye image to 
    polar coordinates (r, theta), then mask out pixels beyond the field of view 
    radius and terrian given the input retangular terrian mask. 
    
    :param fisheye_img: A 2D numpy array representing the fisheye image.
    :param rectangular_mask: A 2D numpy array representing the rectangular mask.
    :param xc: The x-coordinate of the center of the fisheye image.
    :param yc: The y-coordinate of the center of the fisheye image.
    :param radius: The radius for the fisheye image's field of view.

    :return: A 2D numpy array representing the masked fisheye image.
    :rtype: numpy.ndarray
    """

    masked_fisheye = fisheye_img.copy()
    
    # Convert the xy coordinate in the fisheye image to polar r and theta
    x, y = np.meshgrid(*map(np.arange, list(reversed(fisheye_image.shape))))  
    R_grid = np.sqrt((x-xc)**2 + (y-yc)**2).astype(int)
    Theta = np.arctan2(y-yc, x-xc) # -pi to pi
    Theta = (np.where(Theta < 0, Theta + 2*np.pi, Theta)*radius).astype(int) 
    
    # Set pixels beyond the fisheye's FoV to nan
    masked_fisheye[np.where(R_grid>radius)] = np.nan 

    # Mask out the terrian within the FoV
    indices = np.where(~np.isnan(masked_fisheye))   

    for i,j in zip(*indices):
        r = R_grid[i,j]
        theta = Theta[i,j]
        if np.isnan(rectangular_mask[r,theta]):
            masked_fisheye[i,j] = np.nan
            
    return masked_fisheye


def plot(image, jump_ind, color=False, show=False, save_path=None):
    """
    Plot a sequence of masked images to visualize terrain boundaries at different stages.

    This function displays five vertically stacked subplots:
    - The original image
    - Four progressively masked versions based on jump indices

    Args:
        image (ndarray): 2D input image array (e.g., terrain or elevation data).
        jump_ind (list or array-like): List of 4 row indices indicating where masking begins 
                                       in each column for the four masked subplots.
        color (bool, optional): If True, use the 'NPS_mag' colormap with masked pixels shown in black.
                                If False, use grayscale. Defaults to False.
        show (bool, optional): If True, display the plot interactively. Defaults to False.
        save_path (str, optional): Full path to save the output image. If None, saves to a default 
                                   location based on the `color` flag. Defaults to None.

    Returns:
        None
    """

    fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(10.5, 9), sharex=True, sharey=True)
    plt.tight_layout()
    
    # Create masked images
    masked_images = [image.copy()]
    for idx in jump_ind:
        masked = image.copy()
        masked[np.arange(image.shape[0])[:, None] >= idx] = np.nan
        masked_images.append(masked)

    # Set colormap
    if color:
        cmap = copy.copy(mpl.colormaps["NPS_mag"])
        cmap.set_bad(color='black')
        vmin, vmax = 14, 24
    else:
        cmap = 'gray'
        vmin = vmax = None

    # Plot images
    for ax, img in zip(axs, masked_images):
        ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)

    # Save figure
    if save_path:
        filename = save_path
    else:
        filename = p.data_cal+"mask_color.png" if color else p.data_cal+"mask.png"
    plt.savefig(filename, dpi=300)

    # Show figure
    if show:
        plt.show(block=False)
        

#--------------------------------------------------------------------------#
#				    Read in files and image preparation                    #
#--------------------------------------------------------------------------#

#Read in the image center coordinates, and FOV radius
C = pd.read_csv('../Calibration/imagecenter.csv',index_col=0)
Xc = C['Xcenter'][p.camera]
Yc = C['Ycenter'][p.camera] 
R = C['Radius'][p.camera]

# Load grayscale FITS image
ref_img = glob(p.data_cal+'MF_'+p.reference)[0] # Median filtered reference image
fisheye_image = fits.open(ref_img, uint=False)[0].data


#--------------------------------------------------------------------------#
#					       	   Build the mask         					   #
#--------------------------------------------------------------------------#

# Transform the fisheye image to a rectangular image
rectangular_image = fisheye_to_rectangular(fisheye_image, Xc, Yc, R)

# Finding the sky-terrain boundary 
jump_indices = terrain_boundary_mask(rectangular_image)

# Creat the mask in fisheye view
masked = rectangular_image.copy()
masked[np.arange(rectangular_image.shape[0])[:,None] >= jump_indices[3]] = np.nan
masked_fisheye = mask_fisheye_image(fisheye_image, masked, Xc, Yc, R)
mask = np.ones_like(masked_fisheye)
mask[np.where(np.isnan(masked_fisheye))] = np.nan

#--------------------------------------------------------------------------#
#					       	Display the images	    					   #
#--------------------------------------------------------------------------#

# Plot
plot(rectangular_image, jump_indices, color=False, show=True)
plot(rectangular_image, jump_indices, color=True, show=True)


#--------------------------------------------------------------------------#
#					          Save the mask          					   #
#--------------------------------------------------------------------------#

#Save the mask as a fit file
hdu = fits.PrimaryHDU()
hdu.header['XCENTER'] = Xc
hdu.header['YCENTER'] = Yc
hdu.header['RADIUS'] = R
hdu.header['BASEIMG'] = ref_img[18:]
hdu.data = mask
hdu.writeto(p.data_cal+'mask.fit',overwrite=1) #for scientific analysis 

#Save the mask overlaid on the reference fisheye image
hdu.data = masked_fisheye
hdu.writeto(p.data_cal+'mask_fisheye.fit', overwrite=1) #for making a gif

plt.figure()
plt.imshow(masked_fisheye, cmap='gray')
plt.savefig(p.data_cal+'mask_fisheye.png', dpi=300)
plt.show(block=False)