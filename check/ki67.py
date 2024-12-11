import pandas as pd
import numpy as np
import cv2
import fnmatch
from os.path import isfile, join
from skimage.measure import regionprops_table
from os import listdir
import skimage

def intensitiesNucleus(row, output, nucleiMasks, redImage, redImageNorm):
    Img = np.zeros(nucleiMasks.shape, np.uint8)
    center = (int(row['centroid-1']), int(row['centroid-0']))
    radius = int((row['major_axis_length']/2) + 2)
    circle = cv2.circle(Img,center,radius,255,-1)
    nucleus = row['label']
    outside = np.where(nucleiMasks == nucleus, 0, circle)
    inside = np.where(nucleiMasks == nucleus, circle, 0)

    medianOutside = np.median(redImage[np.where(outside == 255)])
    medianOutsideNorm = np.median(redImageNorm[np.where(outside == 255)])
    median = np.median(redImage[np.where(nucleiMasks == nucleus)])
    medianNorm = np.median(redImageNorm[np.where(nucleiMasks == nucleus)])
    if output == 'median':
        return median
    elif output == 'medianNorm':
        return medianNorm
    elif output == 'medianOutside':
        return medianOutside
    elif output == 'medianOutsideNorm':
        return medianOutsideNorm
    else:
        return np.inf

def ki67SingalExtraction(folder, redPath, savePath):
    properties = ('label', 'major_axis_length', 'minor_axis_length', 'orientation', 'centroid')
    df = pd.DataFrame()
    filesNuclei = [f for f in listdir(folder) if isfile(join(folder, f)) if "_seg.npy" in f if "._" not in f]
    aorta = folder.split("/")[-1]
    age = folder.split("/")[-2]
    original = folder.split("/")[-3]

    for fileNuclei in filesNuclei:
        nuclei = np.load(join(folder,fileNuclei), allow_pickle = True).item()
        nucleiMasks = nuclei['masks']
        masks_array = np.asarray(nucleiMasks)
        image = fileNuclei.split("_seg.npy")[0].split('_')[-1]
        fileRed = "ki67Image_" + age + "_ki67_" + image + ".tif"
        folderRed = redPath
        redImage = skimage.io.imread(join(folderRed, fileRed))
        nucleiMaskRedInt = np.where(nucleiMasks > 0, redImage, 0 )
        medianBackground = np.median(redImage[np.where(nucleiMasks == 0)])
        redImageAdjusted = np.where(redImage <= medianBackground, medianBackground, redImage)
        redImageNorm = cv2.normalize(redImageAdjusted, None, alpha=0, beta=1, norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_64F)

        prop_dict = regionprops_table(masks_array, properties = properties)
        df_help = pd.DataFrame(prop_dict)

        if df_help.shape[0]:
            df_help['fractionMajorMinor'] = df_help['major_axis_length']/df_help['minor_axis_length']
            df_help['median'] = df_help.apply(lambda row: intensitiesNucleus(row, 'median', nucleiMasks, redImage, redImageNorm), axis=1)
            df_help['medianNorm'] = df_help.apply(lambda row: intensitiesNucleus(row, 'medianNorm', nucleiMasks, redImage, redImageNorm), axis=1)
            df_help['medianOutside'] = df_help.apply(lambda row: intensitiesNucleus(row, 'medianOutside', nucleiMasks, redImage, redImageNorm), axis=1)
            df_help['medianOutsideNorm'] = df_help.apply(lambda row: intensitiesNucleus(row, 'medianOutsideNorm', nucleiMasks, redImage, redImageNorm), axis=1)
            df_help['fractionNorm'] = df_help.apply(lambda row: fractionNorm(row), axis=1)
            df_help['fraction'] = (df_help['median']/df_help['medianOutside'])
            df_help['original'] = original
            df_help['age'] = age
            df_help['aorta'] = aorta
            df_help['filenameNuc'] = fileNuclei
            df_help['filenameKI67'] = fileRed
            df = pd.concat((df, df_help), axis=0, ignore_index=True)

    filename = "coordinatesKI67Intensity_" + original + "_" + age + "_" + aorta + "_TiledWithOutside.csv"
    df.to_csv(join(savePath,filename), index = False)




def fractionNorm(row):
    if row['medianOutsideNorm'] == 0.000000:
        return row['medianNorm']/(row['medianOutsideNorm'] + 0.00012)
    else:
        return row['medianNorm']/row['medianOutsideNorm']


def thresholding(row, threshold):
    if threshold == "thresholdIntensity":
        thresholdIntensity = thresholds['thresholdIntensity1'].loc[(thresholds['age'] == row['age']) & (thresholds['aorta'] == row['aorta'])].values[0]
        if row['medianNorm'] > thresholdIntensity:
            return 1
        else:
            return 0
    elif threshold == "thresholdFraction":
        thresholdFraction = thresholds['thresholdFraction1'].loc[(thresholds['age'] == row['age']) & (thresholds['aorta'] == row['aorta'])].values[0]
        if row['fractionNormNew'] > thresholdFraction :
            return 1
        else:
            return 0
    else:
        return 99999999999


def threshold(row):
    return np.min([row["thresholdIntensity"], row["thresholdFraction"]])


def xStitched(row, stitched_image):
    widthParam = ((row['imageNumber']-1)%12)
    widthStep = int(stitched_image.size[0]/12)
    return row['centroid-1'] + widthParam*widthStep

def yStitched(row, stitched_image):
    heightParam = int((row['imageNumber']-1)/12)
    heightStep = int(stitched_image.size[1]/12)
    return row['centroid-0'] + heightParam*heightStep
