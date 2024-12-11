import cv2
import numpy as np
from skimage import color
from skimage.segmentation import clear_border

def watershed(imagePath,distTransformMultiplyer):
    img1 = cv2.imread(imagePath)
    img = cv2.cvtColor(img1,cv2.COLOR_BGR2GRAY)
    ret1, thresh = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    kernel = np.ones((3,3),np.uint8)
    opening = cv2.morphologyEx(thresh,cv2.MORPH_OPEN,kernel, iterations = 2)
    opening = thresh
    opening = clear_border(opening)
    sure_bg = cv2.dilate(opening,kernel,iterations=10)
    opening8Bit = opening.astype(np.uint8)
    dist_transform = cv2.distanceTransform(opening8Bit,cv2.DIST_L2,3)
    ret2, sure_fg = cv2.threshold(dist_transform,distTransformMultiplyer*dist_transform.max(),255,0)
    sure_fg = np.uint8(sure_fg)
    unknown = cv2.subtract(sure_bg,sure_fg)
    ret3, markers = cv2.connectedComponents(sure_fg)
    markers = markers+10
    markers[unknown==255] = 0
    markers = cv2.watershed(img1,markers)
    img1[markers == -1] = [0,255,255]
    img2 = color.label2rgb(markers, bg_label=0)
    #connect touching clusters
    markersNew = np.copy(markers)
    markersNew[markersNew == 10] = 0
    minus1 = np.where(markersNew == -1)
    for i in range(len(minus1[0])):
        minus1y = minus1[0][i]
        minus1x = minus1[1][i]
        if (minus1y != 0) & (minus1y != (markersNew.shape[0] -1)) & (minus1x != 0) & (minus1x != (markersNew.shape[1]-1)):
            numbers = np.unique([markersNew[minus1y - 1,minus1x], markersNew[minus1y + 1,minus1x], markersNew[minus1y,minus1x-1], markersNew[minus1y,minus1x +1]])
            uniqueMasks  = set(numbers) - set([-1, 0])
            masks = list(uniqueMasks)
            if len(uniqueMasks) > 1:
                masks = list(uniqueMasks)
                markersNew = np.where(markersNew == masks[1], masks[0], markersNew)
            if 0 not in numbers:
                markersNew[minus1y, minus1x] = masks[0]

    return markersNew