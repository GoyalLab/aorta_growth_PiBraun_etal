import numpy as np
import cv2

def nucleiCheck(cytoPath, nucleiPath, borderpixel, accuracy):
    """
    Check and verify the presence of nuclei within cytoplasm masks based on specified accuracy and border criteria.

    Parameters:
    - cytoPath (str): Path to the cytoplasm mask file (numpy format).
    - nucleiPath (str): Path to the nuclei mask file (numpy format).
    - borderpixel (int): Number of pixels from the image border to exclude in analyses.
    - accuracy (float): Minimum required coverage ratio of nucleus pixels in the cell for the nucleus to be considered inside.

    Returns:
    - cytoOutlines (ndarray): Updated outlines indicating the presence of nucleie within the cytoplasm.
    - cytoMasks (ndarray): Updated masks indicating which cells have nuclei inside.
    - imgCyto (ndarray): Original cytoplasm image array.
    """
    
    # Load nuclei and cytoplasm segmentation data
    nucleiSegmentation = np.load(nucleiPath, allow_pickle=True).item()
    cytoSegmentation = np.load(cytoPath, allow_pickle=True).item()
    nucleiMasks = nucleiSegmentation
    cytoMasks = cytoSegmentation
    cytoOutlines = cytoSegmentation['outlines']
    imgCyto = cytoSegmentation['img']

    # Generate excluded pixel ranges based on borderpixel parameter
    excludeX = np.arange(imgCyto.shape - borderpixel + 1, imgCyto.shape + 1)
    downValues = np.arange(0, borderpixel)
    excludeX = np.concatenate((downValues, excludeX))
    excludeY = np.arange(imgCyto.shape - borderpixel + 1, imgCyto.shape + 1)
    excludeY = np.concatenate((downValues, excludeY))

    # Check for each unique cell mask in cytoplasm
    uniqueCells = np.unique(cytoMasks)
    for c in uniqueCells:
        if c != 0:
            # Get coordinates of current cell in cytoplasm
            cellCoord = np.where(cytoMasks == c)
            # Check if coordinates are within the excluded borders
            checkY = np.isin(cellCoord, excludeY)
            checkX = np.isin(cellCoord, excludeX)
            # If part of the cell is at the image border, exclude it from analysis
            if np.any(checkY) or np.any(checkX):
                cytoMasks = 0
                cytoOutlines = 0
            else:
                # Check which nuclei are associated with this cell
                nuclei = nucleiMasks
                uniqueNuclei = np.unique(nuclei)
                
                # Remove cells with insufficient nuclei or too many unique nuclei
                if len(np.unique(nuclei)) <= 1 or len(np.unique(nuclei)) > 3:
                    cytoMasks = 0
                    cytoOutlines = 0
                else:
                    coverage = 0.0  # Initialize coverage variable to check nucleus pixel coverage
                    for n in uniqueNuclei:
                        if n != 0:
                            nucleiCoord = np.where(nucleiMasks == n)  # Get coordinates of the current nucleus
                            # Count how many pixels of the current nucleus are found in the cell mask
                            cytoNumber = cytoMasks[nucleiCoord, nucleiCoord]
                            cytoNumberNuclei = np.bincount(cytoNumber)
                            # Calculate coverage ratio of nucleus pixels in the cell size
                            if (cytoNumberNuclei / len(cytoNumber) > coverage):
                                coverage = cytoNumberNuclei / len(cytoNumber)

                    # Mask the cell if nucleus coverage is below acceptable accuracy threshold
                    if coverage < accuracy:
                        cytoMasks = 0
                        cytoOutlines[cellCoord, cellCoord] = 0

    return cytoOutlines, cytoMasks, imgCyto  # Returning the updated outlines, masks, and the original image

def nucleiCheckTA(cytoPath, nucleiPath, borderpixel, accuracy):
    """
    Perform a modified nuclei check focusing on the cytoplasm images and their components with additional tallying.

    Parameters:
    - cytoPath (str): Path to the cytoplasm image (as an image file).
    - nucleiPath (str): Path to the nuclei mask file (numpy format).
    - borderpixel (int): Number of pixels from the image border to exclude from analysis.
    - accuracy (float): Minimum required coverage ratio of nucleus pixels in cells for positive identification.

    Returns:
    - cytoMasks (ndarray): Updated masks indicating which cells have valid nuclei inside.
    """
    
    # Initialize counters for exclusion types
    be, ze, he, re = 0, 0, 0, 0
    
    # Load cytoplasm segmentation and convert to binary mask
    cytoSegm = cv2.imread(cytoPath)
    SegmCopy = cytoSegm.copy()  # Make a copy of the cytoplasm image
    SegmCopy = SegmCopy  # Take one channel (assuming single-channel processing for cytoplasm)
    SegmCopy = np.where(SegmCopy == 0, 255, 0)  # Create a binary mask
    SegmCopy = SegmCopy.astype('uint8')  # Ensure it's in the proper datatype
    
    # Find connected components in the cytoplasm mask
    ret3, markers = cv2.connectedComponents(SegmCopy, connectivity=4)
    
    # Load nuclei segmentation data
    nucleiSegmentation = np.load(nucleiPath, allow_pickle=True).item()
    cytoMasks = markers  # Initialize cytoplasm masks from connected components
    nucleiMasks = nucleiSegmentation
    
    # Generate excluded pixel ranges based on borderpixel parameter
    excludeX = np.arange(cytoSegm.shape - borderpixel + 1, cytoSegm.shape + 1)
    downValues = np.arange(0, borderpixel)
    excludeX = np.concatenate((downValues, excludeX))
    excludeY = np.arange(cytoSegm.shape - borderpixel + 1, cytoSegm.shape + 1)
    excludeY = np.concatenate((downValues, excludeY))

    # Check for each unique cell in cytoplasm
    uniqueCells = np.unique(cytoMasks)
    for c in uniqueCells:
        if c != 0:
            cellCoord = np.where(cytoMasks == c)
            checkY = np.isin(cellCoord, excludeY)
            checkX = np.isin(cellCoord, excludeX)

            # Mark and remove cells that are touching the border
            if np.any(checkY) or np.any(checkX):
                cytoMasks = 0
                be += 1  # Increase count for border-excluded
            else:
                nuclei = nucleiMasks
                uniqueNuclei = np.unique(nuclei)

                # Ensure there is at least one nucleus
                if len(np.unique(nuclei)) <= 1:
                    cytoMasks = 0
                    ze += 1  # Increment exclusion count
                elif len(np.unique(nuclei)) >= 3:
                    cytoMasks = 0
                    he += 1  # Increment exclusion count
                else:
                    coverage = 0.0  # Initialize coverage variable
                    cellCounter = 0  # Initialize counter for nucleus entries
                    for n in uniqueNuclei:
                        if n != 0:
                            nucleiCoord = np.where(nucleiMasks == n)
                            # Count how many pixels of the current nucleus are found in the cell mask
                            cytoNumber = cytoMasks
                            cytoNumberNuclei = np.bincount(cytoNumber)
                            # Determine the maximum percentage of how many pixels of the nucleus are within the cell size
                            if (cytoNumberNuclei / len(cytoNumber) > coverage):
                                coverage = cytoNumberNuclei / len(cytoNumber)
                                cellCounter += 1

                    # Mark cells as invalid based on coverage and nucleus conditions
                    if ((coverage < accuracy) or (cellCounter > 1)):
                        cytoMasks = 0
                        re += 1  # Increment exclusion count when conditions are not met
    return cytoMasks  # Return the updated cytoplasm mask

def nucleiCheckTA1(cytoPath, nucleiPath, borderpixel, accuracy):
    """
    Perform a nuclei check with additional verification of the nucleus populations within cytoplasm images.

    Parameters:
    - cytoPath (str): Path to the cytoplasm image (image file).
    - nucleiPath (str): Path to the nuclei mask file (numpy format).
    - borderpixel (int): Pixels from the border to exclude from analysis.
    - accuracy (float): Minimum required coverage ratio for nuclei identification in the cytoplasm.

    Returns:
    - cytoMasks (ndarray): Updated cytoplasm masks indicating valid nuclei.
    """
    
    # Initialize counters for exclusions
    be, ze, he, re = 0, 0, 0, 0
    
    # Load cytoplasm image and convert to binary mask
    cytoSegm = cv2.imread(cytoPath)
    SegmCopy = cytoSegm.copy()  # Copy the cytoplasm image
    SegmCopy = SegmCopy  # Select the channel
    SegmCopy = np.where(SegmCopy == 0, 255, 0)  # Create a binary mask
    SegmCopy = SegmCopy.astype('uint8')  # Set to correct data type

    # Identify connected components in the binary mask
    ret3, markers = cv2.connectedComponents(SegmCopy, connectivity=4)
    
    # Load the nuclei segmentation data
    nucleiSegmentation = np.load(nucleiPath, allow_pickle=True).item()
    cytoMasks = markers  # Initialize cytoplasm masks
    nucleiMasks = nucleiSegmentation['masks']
    
    # Generate excluded pixel ranges based on borderpixel
    excludeX = np.arange(cytoSegm.shape[1] - borderpixel + 1, cytoSegm.shape[1] + 1)
    downValues = np.arange(0, borderpixel)
    excludeX = np.concatenate((downValues, excludeX))
    excludeY = np.arange(cytoSegm.shape[0] - borderpixel + 1, cytoSegm.shape[0] + 1)
    excludeY = np.concatenate((downValues, excludeY))

    # Check for each unique cell in the cytoplasm
    uniqueCells = np.unique(cytoMasks)
    for c in uniqueCells:
        if c != 0:
            cellCoord = np.where(cytoMasks == c)
            checkY = np.isin(cellCoord[0], excludeY)
            checkX = np.isin(cellCoord[1], excludeX)

            # Mark and remove cells touching the border
            if np.any(checkY) or np.any(checkX):
                cytoMasks = 0
                be += 1  # Increment border exclusions
            else:
                nuclei = nucleiMasks
                uniqueNuclei = np.unique(nuclei)

                # Exclude cells with insufficient or problematic nuclei
                if len(np.unique(nuclei)) <= 1:
                    cytoMasks = 0
                    ze += 1  # Increment exclusion count
                elif len(np.unique(nuclei)) >= 3:
                    cytoMasks = 0
                    he += 1  # Increment exclusion count
                else:
                    coverage = 0.0  # Initialize coverage
                    cellCounter = 0  # Initialize counter for nuclei
                    for n in uniqueNuclei:
                        if n != 0:
                            nucleiCoord = np.where(nucleiMasks == n)
                            # Count pixel occurrences of the current nucleus within the cytoplasm
                            cytoNumber = cytoMasks
                            cytoNumberNuclei = np.bincount(cytoNumber)[c]
                            # Check maximal nucleus pixel coverage
                            if (cytoNumberNuclei / len(cytoNumber) > coverage):
                                coverage = cytoNumberNuclei / len(cytoNumber)
                                cellCounter += 1

                    # Mark cells as invalid based on coverage and nucleus count
                    if ((coverage < accuracy) or (cellCounter > 1)):
                        cytoMasks = 0
                        re += 1  # Increment exclusion count
    return cytoMasks  # Return the updated cytoplasm mask

def nucleiCheckTANOQC(cytoPath, nucleiPath, borderpixel, accuracy):
    """
    Perform a nuclei check without quality control.

    Parameters:
    - cytoPath (str): Path to the cytoplasm image (image file).
    - nucleiPath (str): Path to the nuclei mask file (numpy format).
    - borderpixel (int): Pixels from the image border to exclude from analysis.
    - accuracy (float): (Not used; placeholder).

    Returns:
    - cytoMasks (ndarray): Updated masks of cytoplasm.
    """

    # Initialize counters for exclusions
    be, ze, he, re = 0, 0, 0, 0
    
    # Load cytoplasm image and create a binary mask
    cytoSegm = cv2.imread(cytoPath)
    SegmCopy = cytoSegm.copy()  # Make a copy
    SegmCopy = SegmCopy[:, :, 0]  # Take one channel
    SegmCopy = np.where(SegmCopy == 0, 255, 0)  # Convert to binary
    SegmCopy = SegmCopy.astype('uint8')  # Set data type

    # Identify connected components
    ret3, markers = cv2.connectedComponents(SegmCopy, connectivity=4)
    
    # Load nuclei segmentation data
    nucleiSegmentation = np.load(nucleiPath, allow_pickle=True).item()
    cytoMasks = markers  # Assign masks based on connected components
    return cytoMasks  # Return the modified cytoplasm masks