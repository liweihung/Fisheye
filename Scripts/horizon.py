import numpy as np
import cv2
from astropy.io import fits

def detect_horizon_fits(fits_file):
    # Load the FITS file
    with fits.open(fits_file) as hdul:
        # Assuming the image is in the first HDU
        image_data = hdul[0].data

    # Check if the image data is loaded
    if image_data is None:
        print("Error: Image data not found in FITS file.")
        return

    # Handle NaN values by replacing them with the mean of the non-NaN values
    nan_mask = np.isnan(image_data)
    mean_value = np.nanmean(image_data)  # Calculate the mean ignoring NaNs
    image_data[nan_mask] = mean_value    # Replace NaNs with the mean value

    # Normalize the image to 0-255 and convert to uint8
    image_data = (image_data - np.min(image_data)) / (np.max(image_data) - np.min(image_data)) * 255
    image_gray = image_data.astype(np.uint8)

    # Apply Gaussian blur to reduce noise
    blurred = cv2.GaussianBlur(image_gray, (5, 5), 0)

    # Perform edge detection using Canny
    edges = cv2.Canny(blurred, 50, 150)

    # Use Hough Transform to detect lines
    lines = cv2.HoughLinesP(edges, 1, np.pi / 180, threshold=100, minLineLength=100, maxLineGap=10)

    # Create a copy of the image to draw lines
    output_image = cv2.cvtColor(image_gray, cv2.COLOR_GRAY2BGR)

    # Initialize variables to store the best horizon line
    best_horizon_line = None
    max_y = -1  # To keep track of the highest line found

    # If lines are detected, find the horizontal line that represents the horizon
    if lines is not None:
        for line in lines:
            x1, y1, x2, y2 = line[0]
            # Calculate the slope of the line
            slope = (y2 - y1) / (x2 - x1) if (x2 - x1) != 0 else float('inf')

            # Check if the line is approximately horizontal
            if -0.1 < slope < 0.1:  # Adjust this threshold as needed for horizontal lines
                # Check if this line is the highest found so far
                if y1 > max_y:
                    max_y = y1
                    best_horizon_line = (x1, y1, x2, y2)

    # Draw the best horizon line on the output image
    if best_horizon_line is not None:
        x1, y1, x2, y2 = best_horizon_line
        cv2.line(output_image, (x1, y1), (x2, y2), (0, 255, 0), 2)

    # Display the original image and the output image with detected horizon
    cv2.imshow('Original Image', image_gray)
    cv2.imshow('Detected Horizon', output_image)

    # Wait for a key press and close the windows
    cv2.waitKey(0)
    cv2.destroyAllWindows()


from glob import glob
img = "../Data_processed/ROMO_20241003A/MF_Light_Manyparks_30.0s_Bin4_V_20241003-214913_0004.fit"

detect_horizon_fits(img)