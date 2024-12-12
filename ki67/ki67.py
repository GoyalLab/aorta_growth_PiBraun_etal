import pandas as pd
import numpy as np
import cv2
import fnmatch
from os.path import isfile, join
from skimage.measure import regionprops_table
from os import listdir
import skimage

def intensitiesNucleus(row, output, nucleiMasks, redImage, redImageNorm):
    """
    Extracts intensity values from the red image related to a nucleus based on its mask.
    
    Parameters:
    - row (pd.Series): A row containing properties of a nucleus, including the centroid and label.
    - output (str): Determines which intensity value to return ('median', 'medianNorm', 'medianOutside', 'medianOutsideNorm').
    - nucleiMasks (ndarray): Array containing the masks of the nuclei.
    - redImage (ndarray): The raw red image from which intensities are extracted.
    - redImageNorm (ndarray): The normalized version of the red image.
    
    Returns:
    - float: The requested intensity value based on the specified output parameter.
    """
    
    # Create an empty image for the current nucleus mask
    Img = np.zeros(nucleiMasks.shape, np.uint8)
    center = (int(row), int(row))  # Get centroid coordinates
    radius = int((row/2) + 2)  # Define the radius based on nucleus size

    # Create a circular mask around the nucleus
    circle = cv2.circle(Img, center, radius, 255, -1)

    nucleus = row['label']  # Get the label of the nucleus
    outside = np.where(nucleiMasks == nucleus, 0, circle)  # Create mask for outside the nucleus
    inside = np.where(nucleiMasks == nucleus, circle, 0)  # Create mask for inside the nucleus

    # Calculate median intensities for regions
    medianOutside = np.median(redImage)
    medianOutsideNorm = np.median(redImageNorm[np.where(outside == 255)])
    median = np.median(redImage)
    medianNorm = np.median(redImageNorm[np.where(nucleiMasks == nucleus)])
    
    # Return the requested intensity value
    if output == 'median':
        return median
    elif output == 'medianNorm':
        return medianNorm
    elif output == 'medianOutside':
        return medianOutside
    elif output == 'medianOutsideNorm':
        return medianOutsideNorm
    else:
        return np.inf  # Return infinity for unexpected output types

def ki67SingalExtraction(folder, redPath, savePath):
    """
    Extracts Ki67 signal intensities from nuclei masks contained in a specified folder.
    
    Parameters:
    - folder (str): The path to the folder containing the nuclei segmentation files.
    - redPath (str): The path to the folder containing the red channel images.
    - savePath (str): The path to save the extracted intensity results.
    
    Returns:
    - None: Saves the output DataFrame as a CSV file.
    """
    
    properties = ('label', 'major_axis_length', 'minor_axis_length', 'orientation', 'centroid')
    df = pd.DataFrame()  # Initialize an empty DataFrame to hold results
    filesNuclei = [f for f in listdir(folder) if isfile(join(folder, f)) and "_seg.npy" in f and "._" not in f]
    
    # Extract metadata from the folder path
    aorta = folder.split("/")
    age = folder.split("/")[-2]
    original = folder.split("/")[-3]

    # Process each nuclei segmentation file
    for fileNuclei in filesNuclei:
        # Load the nuclei masks from the segmentation file
        nuclei = np.load(join(folder, fileNuclei), allow_pickle=True).item()
        nucleiMasks = nuclei['masks']  # Extract masks array
        masks_array = np.asarray(nucleiMasks)
        
        # Construct the file path for the corresponding red image
        image = fileNuclei.split("_seg.npy").split('_')[-1]
        fileRed = f"ki67Image_{age}_ki67_{image}.tif"
        folderRed = redPath
        
        # Read the red image
        redImage = skimage.io.imread(join(folderRed, fileRed))
        
        # Create a masked red image where nuclei are present
        nucleiMaskRedInt = np.where(nucleiMasks > 0, redImage, 0)
        medianBackground = np.median(redImage[np.where(nucleiMasks == 0)])  # Median of the background
        redImageAdjusted = np.where(redImage <= medianBackground, medianBackground, redImage)  # Adjust the red image
        redImageNorm = cv2.normalize(redImageAdjusted, None, alpha=0, beta=1, norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_64F)  # Normalize the image

        # Get properties of the nuclei
        prop_dict = regionprops_table(masks_array, properties=properties)
        df_help = pd.DataFrame(prop_dict)

        if df_help.shape:  # If there are nuclei found
            # Calculate additional properties based on intensity analysis
            df_help['fractionMajorMinor'] = df_help['major_axis_length'] / df_help['minor_axis_length']
            # Apply intensity calculations for each nucleus
            df_help = df_help.apply(lambda row: intensitiesNucleus(row, 'median', nucleiMasks, redImage, redImageNorm), axis=1)
            df_help = df_help.apply(lambda row: intensitiesNucleus(row, 'medianNorm', nucleiMasks, redImage, redImageNorm), axis=1)
            df_help = df_help.apply(lambda row: intensitiesNucleus(row, 'medianOutside', nucleiMasks, redImage, redImageNorm), axis=1)
            df_help = df_help.apply(lambda row: intensitiesNucleus(row, 'medianOutsideNorm', nucleiMasks, redImage, redImageNorm), axis=1)
            df_help['fractionNorm'] = df_help.apply(lambda row: fractionNorm(row), axis=1)
            df_help['fraction'] = (df_help['median'] / df_help['medianOutside'])
            df_help['original'] = original  # Track origin information
            df_help['age'] = age  # Track age information
            df_help['aorta'] = aorta  # Track aorta information
            df_help['filenameNuc'] = fileNuclei  # Track filename for nuclei segmentation
            df_help['filenameKI67'] = fileRed  # Track filename for Ki67 image
            
            # Concatenate the results into the main DataFrame
            df = pd.concat((df, df_help), axis=0, ignore_index=True)

    # Save the results to a CSV file
    filename = f"coordinatesKI67Intensity_{original}_{age}_{aorta}_TiledWithOutside.csv"
    df.to_csv(join(savePath, filename), index=False)

def fractionNorm(row):
    """
    Normalizes the median intensity of the nucleus against the median intensity outside.
    
    Parameters:
    - row (pd.Series): Contains median intensities for the nucleus and outside.
    
    Returns:
    - float: Normalized fraction of the nucleus median intensity to outside median intensity.
    """
    
    if row == 0.000000:
        return row / (row + 0.00012)  # Avoid division by zero
    else:
        return row / row['medianOutsideNorm']  # Normalization calculation

def thresholding(row, threshold):
    """
    Applies thresholds to the median intensity and normalized fraction.
    
    Parameters:
    - row (pd.Series): Contains median and fraction values.
    - threshold (str): Type of thresholding to perform ('thresholdIntensity' or 'thresholdFraction').
    
    Returns:
    - int: 1 if the row passes the threshold, 0 otherwise.
    """
    
    if threshold == "thresholdIntensity":
        thresholdIntensity = thresholds['thresholdIntensity1'].loc.values
        return 1 if row['medianNorm'] > thresholdIntensity else 0
    elif threshold == "thresholdFraction":
        thresholdFraction = thresholds['thresholdFraction1'].loc[
            (thresholds['age'] == row['age']) & (thresholds['aorta'] == row['aorta'])
        ].values
        return 1 if row['fractionNormNew'] > thresholdFraction else 0
    else:
        return 99999999999  # Return a high value for unrecognized thresholds

def threshold(row):
    """
    Computes minimum threshold value from two threshold types.
    
    Parameters:
    - row (pd.Series): Contains two threshold values.
    
    Returns:
    - int: Minimum of the two threshold types.
    """
    
    return np.min([row["thresholdIntensity"], row["thresholdFraction"]])

def xStitched(row, stitched_image):
    """
    Computes the x-coordinate for an image in a stitched layout.
    
    Parameters:
    - row (pd.Series): Contains the image number and centroid x-coordinate.
    - stitched_image (PIL.Image): The stitched image to determine dimensions.
    
    Returns:
    - float: Adjusted x-coordinate for the stitched image.
    """
    
    widthParam = ((row - 1) % 12)  # Determine position along the width
    widthStep = int(stitched_image.size[0] / 12)  # Calculate step size in the x direction
    return row['centroid-1'] + widthParam * widthStep

def yStitched(row, stitched_image):
    """
    Computes the y-coordinate for an image in a stitched layout.
    
    Parameters:
    - row (pd.Series): Contains the image number and centroid y-coordinate.
    - stitched_image (PIL.Image): The stitched image to determine dimensions.
    
    Returns:
    - float: Adjusted y-coordinate for the stitched image.
    """
    
    heightParam = int((row['imageNumber'] - 1) / 12)  # Determine position along the height
    heightStep = int(stitched_image.size[1] / 12)  # Calculate step size in the y direction
    return row['centroid-0'] + heightParam * heightStep