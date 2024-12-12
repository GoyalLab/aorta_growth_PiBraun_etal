import cv2
import numpy as np
from skimage import color
from skimage.segmentation import clear_border

def watershed(imagePath, distTransformMultiplyer):
    """
    Performs watershed segmentation on an image.

    Parameters:
    - imagePath (str): The file path of the input image.
    - distTransformMultiplyer (float): A multiplier for the distance transform 
      to determine the sure foreground threshold.

    Returns:
    - markersNew (ndarray): The updated markers array after applying watershed 
      segmentation, where different segments are labeled with distinct integers.
    """
    
    # Read the image from the given file path
    img1 = cv2.imread(imagePath)

    # Convert the image to grayscale
    img = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)
    
    # Apply Otsu's thresholding to obtain a binary image
    ret1, thresh = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    # Create a kernel for morphological operations
    kernel = np.ones((3, 3), np.uint8)

    # Perform morphological opening to remove small objects
    opening = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel, iterations=2)

    # Clear objects connected to the border
    opening = clear_border(opening)

    # Dilate the binary image to get the sure background
    sure_bg = cv2.dilate(opening, kernel, iterations=10)

    # Convert the opening image to an 8-bit image
    opening8Bit = opening.astype(np.uint8)

    # Compute the distance transform
    dist_transform = cv2.distanceTransform(opening8Bit, cv2.DIST_L2, 3)

    # Create a binary image for sure foreground based on a threshold
    ret2, sure_fg = cv2.threshold(dist_transform, distTransformMultiplyer * dist_transform.max(), 255, 0)
    sure_fg = np.uint8(sure_fg)  # Ensure the sure_fg is of uint8 type
    
    # Identify unknown regions by subtracting sure foreground from sure background
    unknown = cv2.subtract(sure_bg, sure_fg)

    # Label the sure foreground markers
    ret3, markers = cv2.connectedComponents(sure_fg)

    # Increment the markers to ensure sure_fg has positive values
    markers = markers + 10
    
    # Mark unknown regions (those not in sure_fg) as zero
    markers[unknown == 255] = 0

    # Apply the watershed algorithm to segment the image
    markers = cv2.watershed(img1, markers)

    # Mark boundaries of segments in the original image
    img1[markers == -1] = [0, 255, 255]  # Set boundary pixels to cyan

    # Convert marker labels to an RGB image for visualization
    img2 = color.label2rgb(markers, bg_label=0)

    # Connect touching clusters by modifying the markers
    markersNew = np.copy(markers)
    markersNew[markersNew == 10] = 0  # Set initial marker values to zero for connected components

    # Identify the markers that are defined as boundaries
    minus1 = np.where(markersNew == -1)
    
    # Iterate through boundary pixels to connect clusters
    for i in range(len(minus1)):
        minus1y = minus1
        minus1x = minus1[i]
        
        # Ensure the current point is not on the image border
        if (minus1y != 0) and (minus1y != (markersNew.shape - 1)) and \
           (minus1x != 0) and (minus1x != (markersNew.shape - 1)):
            # Find unique marker values in neighboring pixels
            numbers = np.unique([
                markersNew[minus1y - 1, minus1x], 
                markersNew[minus1y + 1, minus1x], 
                markersNew[minus1y, minus1x - 1], 
                markersNew[minus1y, minus1x + 1]
            ])
            
            # Exclude the background and boundary marker values
            uniqueMasks = set(numbers) - set([-1, 0])
            masks = list(uniqueMasks)

            # If more than one unique mask is found, update the markers
            if len(uniqueMasks) > 1:
                masks = list(uniqueMasks)
                markersNew = np.where(markersNew == masks[1], masks, markersNew)
            
            # If the point is not connected to the background (0), assign it to a mask
            if 0 not in numbers:
                markersNew[minus1y, minus1x] = masks[0]

    return markersNew  # Return the updated markers for further processing