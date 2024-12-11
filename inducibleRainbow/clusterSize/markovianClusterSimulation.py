import random
import numpy as np


def markovianCluster(markedClusters, dividedCells, simulations):
    random.seed(42)
    start = list(range(1, markedClusters + 1))
    resultsRandomCluster = []
    for i in range(simulations):
        randomCluster = start.copy()
        for i in range(dividedCells):
            randomCluster.append(random.randrange(1, markedClusters + 1))
        clusterSizesRandomCluster = []
        for i in range(1,markedClusters +1):
            clusterSizesRandomCluster.append(randomCluster.count(i))
        resultsRandomCluster.append(clusterSizesRandomCluster)
    return resultsRandomCluster

def markovianDivision(markedClusters, dividedCells, simulations):
    random.seed(42)
    start = list(range(1, markedClusters + 1))
    resultsRandomDivision = []
    for i in range(simulations):
        randomDivision = start.copy()
        for j in range(dividedCells):
            randomDivision.append(random.choice(randomDivision))
        clusterSizesRandomDivision = []
        for i in range(1,markedClusters + 1):
            clusterSizesRandomDivision.append(randomDivision.count(i))
        resultsRandomDivision.append(clusterSizesRandomDivision)
    return resultsRandomDivision


def markovianClusterInfer(clusterSizesStart, clusterSizesEnd, simulations):
    c = 1
    for i in range(len(clusterSizesStart)):
        if i == 0:
            clusterStart = np.full(clusterSizesStart[i], c)
        else:
            clusterStart = np.concatenate([clusterStart, np.full(clusterSizesStart[i], c)])
        c += 1

    markedClusters = int(len(np.unique(clusterStart)))

    dividedCells = 0
    for i in range(len(clusterSizesEnd)):
        dividedCells += clusterSizesEnd[i] -1

    random.seed(42)
    resultsRandomCluster = []
    for i in range(simulations):
        randomCluster = clusterStart.copy()
        for i in range(dividedCells):
            randomCluster = np.append(randomCluster,random.randrange(1, markedClusters + 1))
        unique, counts = np.unique(randomCluster, return_counts=True)
        resultsRandomCluster.append(counts)
    return resultsRandomCluster


def markovianDivisionInfer(clusterSizesStart, clusterSizesEnd, simulations):
    c = 1
    for i in range(len(clusterSizesStart)):
        if i == 0:
            clusterStart = np.full(clusterSizesStart[i], c)
        else:
            clusterStart = np.concatenate([clusterStart, np.full(clusterSizesStart[i], c)])
        c += 1

    dividedCells = 0
    for i in range(len(clusterSizesEnd)):
        dividedCells += clusterSizesEnd[i] -1

    random.seed(42)
    resultsRandomDivision = []
    for i in range(simulations):
        randomDivision = clusterStart.copy()
        for j in range(dividedCells):
            randomDivision = np.append(randomDivision, random.choice(randomDivision))
        unique, counts = np.unique(randomDivision, return_counts=True)
        resultsRandomDivision.append(counts)
    return resultsRandomDivision