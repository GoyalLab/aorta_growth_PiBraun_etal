import numpy as np
import cv2


def nucleiCheck(cytoPath, nucleiPath, borderpixel, accuracy):
    #input: path cytoplasm, path nuclei; pixel tube for border pixel
    nucleiSegmentation = np.load(nucleiPath, allow_pickle=True).item()
    cytoSegmentation = np.load(cytoPath, allow_pickle=True).item()
    nucleiMasks = nucleiSegmentation['masks']
    cytoMasks = cytoSegmentation['masks']
    cytoOutlines = cytoSegmentation['outlines']
    imgCyto = cytoSegmentation['img']
    excludeX = np.arange(imgCyto.shape[1] - borderpixel + 1 ,imgCyto.shape[1] +1)
    downValues = np.arange(0,borderpixel)
    excludeX = np.concatenate((downValues, excludeX))
    excludeY = np.arange(imgCyto.shape[0] - borderpixel +1 ,imgCyto.shape[0] +1)
    excludeY = np.concatenate((downValues, excludeY))
    # nucleicheck
    uniqueCells = np.unique(cytoMasks)
    for c in uniqueCells:
        if c != 0:
            cellCoord = np.where(cytoMasks==c)
            checkY = np.isin(cellCoord[0],excludeY)
            checkX = np.isin(cellCoord[1],excludeX)
            if np.any(checkY):
                cytoMasks[cellCoord[0], cellCoord[1]] = 0
                cytoOutlines[cellCoord[0], cellCoord[1]] = 0
            elif np.any(checkX):
                cytoMasks[cellCoord[0], cellCoord[1]] = 0
                cytoOutlines[cellCoord[0], cellCoord[1]] = 0
            else:
                nuclei = nucleiMasks[cellCoord[0], cellCoord[1]]
                uniqueNuclei = np.unique(nuclei)
                if len(np.unique(nuclei)) <=1:
                    cytoMasks[cellCoord[0], cellCoord[1]] = 0
                    cytoOutlines[cellCoord[0], cellCoord[1]] = 0
                elif len(np.unique(nuclei)) > 3:
                    cytoMasks[cellCoord[0], cellCoord[1]] = 0
                    cytoOutlines[cellCoord[0], cellCoord[1]] = 0
                else:
                    coverage = 0.0
                    for n in uniqueNuclei:
                        if n!= 0:
                            nucleiCoord = np.where(nucleiMasks==n)
                            #How many pixels does this nucleus have:
                            cytoNumber = cytoMasks[nucleiCoord[0], nucleiCoord[1]]
                            #occurence of nucleus pixel in cyto
                            cytoNumberNuclei = np.bincount(cytoNumber)[c]
                            # max percentage of how many pixels of nuc are in cellSize
                            if (cytoNumberNuclei/len(cytoNumber) > coverage):
                                coverage = cytoNumberNuclei/len(cytoNumber)

                    if (coverage < accuracy):
                        cytoMasks[cellCoord[0], cellCoord[1]] = 0
                        cytoOutlines[cellCoord[0], cellCoord[1]] = 0
    return cytoOutlines, cytoMasks, imgCyto


def nucleiCheckTA(cytoPath, nucleiPath, borderpixel, accuracy):
    be = 0
    ze = 0
    he = 0
    re = 0
    #input: path cytoplasm, path nuclei; pixel tube for border pixel
    cytoSegm = cv2.imread(cytoPath)
    SegmCopy = cytoSegm.copy()
    SegmCopy = SegmCopy[:,:,0]
    SegmCopy = np.where(SegmCopy==0, 255,0)
    SegmCopy = SegmCopy.astype('uint8')
    ret3, markers = cv2.connectedComponents(SegmCopy, connectivity = 4)
    nucleiSegmentation = np.load(nucleiPath, allow_pickle=True).item()
    cytoMasks = markers
    nucleiMasks = nucleiSegmentation['masks']
    excludeX = np.arange(cytoSegm.shape[1] - borderpixel + 1 ,cytoSegm.shape[1] +1)
    downValues = np.arange(0,borderpixel)
    excludeX = np.concatenate((downValues, excludeX))
    excludeY = np.arange(cytoSegm.shape[0] - borderpixel +1 ,cytoSegm.shape[0] +1)
    excludeY = np.concatenate((downValues, excludeY))
    # nucleicheck
    uniqueCells = np.unique(cytoMasks)
    for c in uniqueCells:
        if c != 0:
            cellCoord = np.where(cytoMasks==c)
            checkY = np.isin(cellCoord[0],excludeY)
            checkX = np.isin(cellCoord[1],excludeX)
            if np.any(checkY):
                cytoMasks[cellCoord[0], cellCoord[1]] = 0
                be += 1
            elif np.any(checkX):
                cytoMasks[cellCoord[0], cellCoord[1]] = 0
                be += 1
            else:
                nuclei = nucleiMasks[cellCoord[0], cellCoord[1]]
                uniqueNuclei = np.unique(nuclei)
                if len(np.unique(nuclei)) <=1:
                    cytoMasks[cellCoord[0], cellCoord[1]] = 0
                    ze += 1
                elif len(np.unique(nuclei)) > 3:
                    cytoMasks[cellCoord[0], cellCoord[1]] = 0
                    he += 1
                else:
                    coverage = 0.0
                    for n in uniqueNuclei:
                        if n!= 0:
                            nucleiCoord = np.where(nucleiMasks==n)
                            #How many pixels does this nucleus have:
                            cytoNumber = cytoMasks[nucleiCoord[0], nucleiCoord[1]]
                            #occurence of nucleus pixel in cyto
                            cytoNumberNuclei = np.bincount(cytoNumber)[c]
                            # max percentage of how many pixels of nuc are in cellSize
                            if (cytoNumberNuclei/len(cytoNumber) > coverage):
                                coverage = cytoNumberNuclei/len(cytoNumber)

                    if (coverage < accuracy):
                        cytoMasks[cellCoord[0], cellCoord[1]] = 0
                        re += 1
    return cytoMasks



def nucleiCheckTA1(cytoPath, nucleiPath, borderpixel, accuracy):
    be = 0
    ze = 0
    he = 0
    re = 0
    #input: path cytoplasm, path nuclei; pixel tube for border pixel
    cytoSegm = cv2.imread(cytoPath)
    SegmCopy = cytoSegm.copy()
    SegmCopy = SegmCopy[:,:,0]
    SegmCopy = np.where(SegmCopy==0, 255,0)
    SegmCopy = SegmCopy.astype('uint8')
    ret3, markers = cv2.connectedComponents(SegmCopy, connectivity = 4)
    nucleiSegmentation = np.load(nucleiPath, allow_pickle=True).item()
    cytoMasks = markers
    nucleiMasks = nucleiSegmentation['masks']
    excludeX = np.arange(cytoSegm.shape[1] - borderpixel + 1 ,cytoSegm.shape[1] +1)
    downValues = np.arange(0,borderpixel)
    excludeX = np.concatenate((downValues, excludeX))
    excludeY = np.arange(cytoSegm.shape[0] - borderpixel +1 ,cytoSegm.shape[0] +1)
    excludeY = np.concatenate((downValues, excludeY))
    # nucleicheck
    uniqueCells = np.unique(cytoMasks)
    for c in uniqueCells:
        if c != 0:
            cellCoord = np.where(cytoMasks==c)
            checkY = np.isin(cellCoord[0],excludeY)
            checkX = np.isin(cellCoord[1],excludeX)
            if np.any(checkY):
                cytoMasks[cellCoord[0], cellCoord[1]] = 0
                be += 1
            elif np.any(checkX):
                cytoMasks[cellCoord[0], cellCoord[1]] = 0
                be += 1
            else:
                nuclei = nucleiMasks[cellCoord[0], cellCoord[1]]
                uniqueNuclei = np.unique(nuclei)
                if len(np.unique(nuclei)) <=1:
                    cytoMasks[cellCoord[0], cellCoord[1]] = 0
                    ze += 1
                elif len(np.unique(nuclei)) >= 3:
                    cytoMasks[cellCoord[0], cellCoord[1]] = 0
                    he += 1
                else:
                    coverage = 0.0
                    cellCounter = 0
                    for n in uniqueNuclei:
                        if n!= 0:
                            nucleiCoord = np.where(nucleiMasks==n)
                            #How many pixels does this nucleus have:
                            cytoNumber = cytoMasks[nucleiCoord[0], nucleiCoord[1]]
                            #occurence of nucleus pixel in cyto
                            cytoNumberNuclei = np.bincount(cytoNumber)[c]
                            # max percentage of how many pixels of nuc are in cellSize
                            if (cytoNumberNuclei/len(cytoNumber) > coverage):
                                coverage = cytoNumberNuclei/len(cytoNumber)
                                cellCounter += 1

                    if ((coverage < accuracy) or (cellCounter > 1)):
                        cytoMasks[cellCoord[0], cellCoord[1]] = 0
                        re += 1
    return cytoMasks


def nucleiCheckTANOQC(cytoPath, nucleiPath, borderpixel, accuracy):
    be = 0
    ze = 0
    he = 0
    re = 0
    #input: path cytoplasm, path nuclei; pixel tube for border pixel
    cytoSegm = cv2.imread(cytoPath)
    SegmCopy = cytoSegm.copy()
    SegmCopy = SegmCopy[:,:,0]
    SegmCopy = np.where(SegmCopy==0, 255,0)
    SegmCopy = SegmCopy.astype('uint8')
    ret3, markers = cv2.connectedComponents(SegmCopy, connectivity = 4)
    nucleiSegmentation = np.load(nucleiPath, allow_pickle=True).item()
    cytoMasks = markers
    return cytoMasks