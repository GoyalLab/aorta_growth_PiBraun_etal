import numpy as np
from os.path import join
from sklearn.neighbors import NearestNeighbors
from random import seed
from random import sample
import pandas as pd

def stitchedNucImage(nucleiPath, scalingFactor = 10000):
    for i in range(4):
    #build row block
        nuc1 =np.load(join(nucleiPath, "Cropped_IMG-" + str(4*i +1) + '_seg.npy'), allow_pickle=True).item()['masks']
        nuc2 =np.load(join(nucleiPath, "Cropped_IMG-" + str(4*i +2) + '_seg.npy'), allow_pickle=True).item()['masks']
        nuc3 =np.load(join(nucleiPath, "Cropped_IMG-" + str(4*i +3) + '_seg.npy'), allow_pickle=True).item()['masks']
        nuc4 =np.load(join(nucleiPath, "Cropped_IMG-" + str(4*i +4) + '_seg.npy'), allow_pickle=True).item()['masks']
        nuc1 = nuc1 + scalingFactor *(4*i +1)
        nuc1[nuc1 == (scalingFactor *(4*i +1))] = 0
        nuc2 = nuc2 + scalingFactor *(4*i +2)
        nuc2[nuc2 == (scalingFactor *(4*i +2))] = 0
        nuc3 = nuc3 + scalingFactor *(4*i +3)
        nuc3[nuc3 == (scalingFactor *(4*i +3))] = 0
        nuc4 = nuc4 + scalingFactor *(4*i +4)
        nuc4[nuc4 == (scalingFactor *(4*i +4))] = 0

        row = np.concatenate((nuc1, nuc2, nuc3, nuc4), axis = 1)
        if i == 0:
            nucImage = row
        else:
            nucImage = np.concatenate((nucImage, row), axis = 0)
    return nucImage

def centeroidPoints(arr):
    length = arr.shape[0]
    sum_y = np.sum(arr[:, 0])
    sum_x = np.sum(arr[:, 1])
    return sum_y/length, sum_x/length

def clusterArrays(cluster, df):
    singletons = []
    doublets = []
    higher = []
    marked = []
    for key in cluster:
        if len(cluster[key]) == 1:
            singletons.append(df[df['label'] == cluster[key][0]][['centroid-0', 'centroid-1']].to_numpy()[0])
            marked.append(df[df['label'] == cluster[key][0]][['centroid-0', 'centroid-1']].to_numpy()[0])
        elif len(cluster[key]) == 2:
            points = np.asarray([df[df['label'] == cluster[key][i]][['centroid-0', 'centroid-1']].to_numpy()[0] for i in range(len(cluster[key]))])
            centroid = centeroidPoints(points)
            doublets.append(centroid)
            marked.append(centroid)
        else:
            points = np.asarray([df[df['label'] == cluster[key][i]][['centroid-0', 'centroid-1']].to_numpy()[0] for i in range(len(cluster[key]))])
            centroid = centeroidPoints(points)
            higher.append(centroid)
            marked.append(centroid)

    return np.asarray(marked), np.asarray(singletons), np.asarray(doublets), np.asarray(higher)

def kNN(array, neighbor):
    nbrs = NearestNeighbors(n_neighbors=neighbor, algorithm='ball_tree').fit(array)
    distances, indices = nbrs.kneighbors(array)
    return nbrs, distances, indices

def randomkNNDistances(nucIDs, len, df, neighbors=6):
    selectedRandom = sample(nucIDs, len)
    dfRandom = df[df['label'].isin(selectedRandom)]
    nucCoordRandom = dfRandom[['centroid-0', 'centroid-1']].to_numpy()
    nbrsRandom = NearestNeighbors(n_neighbors=neighbors, algorithm='ball_tree').fit(nucCoordRandom)
    distancesRand, indicesRandom = nbrsRandom.kneighbors(nucCoordRandom)
    return distancesRand

def randomClusters(markedArray, singletonsArray, doubletsArray, higherArray, nucIDs, df, simulations = 1000, seedNumber = 42, higher = True):
    seed(seedNumber)
    markedRandom = []
    singletonRandom = []
    doubletRandom = []
    higherRandom = []
    for _ in range(simulations):
        markedRandom.append(np.asarray(randomkNNDistances(nucIDs,len(markedArray), df)))
        singletonRandom.append(np.asarray(randomkNNDistances(nucIDs,len(singletonsArray), df)))
        doubletRandom.append(np.asarray((randomkNNDistances(nucIDs, len(doubletsArray), df))))
        if higher:
            higherRandom.append(np.asarray(randomkNNDistances(nucIDs, len(higherArray), df)))
    return np.asarray(markedRandom), np.asarray(singletonRandom), np.asarray(doubletRandom), np.asarray(higherRandom)

def getClusterSizes(df):
    dic = {}
    for index, row in df.iterrows():
        aorta = row["aorta"]
        for i in range(2,16):
            cluster = np.full(row[i], i-1)
            if i==2:
                clusterSize = cluster
            else:
                clusterSize = np.concatenate([clusterSize, cluster])
        dic.update({aorta: clusterSize})
    return dic


