import numpy as np
from dictances import bhattacharyya
from scipy.special import rel_entr
from sklearn.metrics import mutual_info_score
import sys

def clusterSizeDF(row, withSingletons):
    #2 with singletons, 3 for divided
    if withSingletons == True:
        start = 2
    else:
        start = 3
    clusterSizes = np.full(row[start], 1)
    for i in range(start+1,17):
        cluster = np.full(row[i], i-1)
        clusterSizes = np.append(clusterSizes, cluster)
    return np.asarray(clusterSizes)

def stats(clusterSizes):
    dividend = 0
    divisor = 0
    shannonlog2 = 0
    shannonlogN = 0
    simpson = 0
    mean = np.mean(clusterSizes)
    n = len(clusterSizes)
    hooverdivid = 0
    hooverdiv = 0

    for i in range(len(clusterSizes)):
        divisor += clusterSizes[i]

        shannonlog2 += clusterSizes[i]/n * np.log2(clusterSizes[i]/n)
        shannonlogN += clusterSizes[i]/n * np.log(clusterSizes[i]/n)

        simpson += (clusterSizes[i]/n)**2

        hooverdivid += np.abs(clusterSizes[i] - mean)
        hooverdiv += clusterSizes[i]

        for j in range(len(clusterSizes)):
            dividend += np.abs(clusterSizes[i] - clusterSizes[j])
    gini = dividend/(2*len(clusterSizes)*divisor)
    hoover = hooverdivid/(2*hooverdiv)
    shannonlog2 = -1* shannonlog2
    shannonlogN = -1* shannonlogN
    shannonEquiLog2 = shannonlog2/np.log2(n)
    shannonEquiLogN = shannonlogN/np.log(n)
    return gini, shannonlog2, shannonlogN, shannonEquiLog2, shannonEquiLogN, simpson, hoover

def getClusterSizes(dataPath, savePath):
    cluster = np.load(dataPath, allow_pickle=True).item()
    clusterSizes = []
    for key in cluster:
        clusterSizes.append(len(cluster[key]))
    df = pd.DataFrame({"C5RB-044": clusterSizes})
    df.to_csv(savePath)

def kl_divergenceLog2(p, q):
    pAna = p.copy()
    qAna = q.copy()
    for i in range(len(p)):
        if pAna[i] == 0:
            pAna[i] = sys.float_info.epsilon
        if qAna[i] == 0:
            qAna[i] = sys.float_info.epsilon
    return sum(pAna[i] * np.log2(pAna[i]/qAna[i]) for i in range(len(pAna)))

def kl_divergenceLog(p, q):
    pAna = p.copy()
    qAna = q.copy()
    for i in range(len(p)):
        if pAna[i] == 0:
            pAna[i] = sys.float_info.epsilon
        if qAna[i] == 0:
            qAna[i] = sys.float_info.epsilon
    return sum(pAna[i] * log(pAna[i]/qAna[i]) for i in range(len(pAna)))

def distanceScores(p,q):
    p = p/sum(p)
    q = q/sum(q)
    jensenShannon = scipy.spatial.distance.jensenshannon(p,q)
    wasserstein = wasserstein_distance(p,q)
    kullbackLeibler = kl_divergenceLog2(p, q)
    return {'JensenShannon': jensenShannon, 'Wasserstein': wasserstein, 'KullbackLeibler': kullbackLeibler}