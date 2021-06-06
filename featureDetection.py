#!/usr/bin/python

import os, sys
import numpy as np, pandas as pd
from pyteomics import mzxml
from datetime import datetime
from peaks import *
from utils import *


def filterFeatures(f):
    # Input arguement, f = array of features
    # A feature may contain multiple peaks from one scan
    # In this case, one with the largest intensity is chosen
    gMinRt, gMaxRt = 0, 0  # Global minimum and maximum RT over all features
    for i in range(len(f)):
        if len(f[i].ms1) != len(list(set(f[i].ms1))):
            temp = {}
            for j in range(len(f[i].ms1)):
                if f[i].ms1[j] in temp:
                    currIntensity = f[i].intensity[j]
                    if currIntensity > temp[f[i].ms1[j]]["intensity"]:
                        temp[f[i].ms1[j]]["intensity"] = currIntensity
                        temp[f[i].ms1[j]]["index"] = j
                else:
                    temp[f[i].ms1[j]] = {}
                    temp[f[i].ms1[j]]["intensity"] = f[i].intensity[j]
                    temp[f[i].ms1[j]]["index"] = j
            uInd = []
            for key in sorted(temp.keys()):
                uInd.append(temp[key]["index"])
            f[i].mz = [f[i].mz[u] for u in uInd]
            f[i].intensity = [f[i].intensity[u] for u in uInd]
            f[i].ms1 = [f[i].ms1[u] for u in uInd]
            f[i].rt = [f[i].rt[u] for u in uInd]
            f[i].index = [f[i].index[u] for u in uInd]

    return f


def organizeFeatures(f, noise, params):
    n = 0
    ms1ToFeatures = {}
    gMinRt, gMaxRt = 0, 0
    for i in range(len(f)):
        # Global minimum and maximum RT values need to be obtained (for later use of "percentage of true feature")
        if i == 0:
            gMinRt = min(f[i].rt)
            gMaxRt = max(f[i].rt)
        else:
            if min(f[i].rt) < gMinRt:
                gMinRt = min(f[i].rt)
            if max(f[i].rt) > gMaxRt:
                gMaxRt = max(f[i].rt)

        # 1. mz: mean m/z of a feauture = weighted average of m/z and intensity
        # 2. intensity: intensity of a feature (maximum intensity among the peaks consist of the feature)
        # 3. z: charge of the feature, set to 1 now, but modified later
        # 4. rt: RT of the representative peak (i.e. strongest peak) of a feature
        # 5. minRt and maxRt
        # 6. ms1: MS1 scan number of the representative peak of a feature
        # 7. minMs1 and maxMs1
        # 8. snRatio: signal-to-noise ratio of the feature
        # 9. pctTf: percentage of the true feature

        mz = np.sum(np.multiply(f[i].mz, f[i].intensity)) / np.sum(f[i].intensity)
        intensity = max(f[i].intensity)
        z = 1
        # isotope = 0  # Will be used later
        ind = np.argmax(f[i].intensity)
        rt = f[i].rt[ind]
        minRt = min(f[i].rt)
        maxRt = max(f[i].rt)
        ms1 = f[i].ms1[ind]
        minMs1 = min(list(map(int, f[i].ms1)))
        maxMs1 = max(list(map(int, f[i].ms1)))
        if ms1 in noise:
            noiseLevel = noise[ms1]
        else:
            noiseLevel = 500
        snRatio = intensity / noiseLevel
        featureIntensityThreshold = noiseLevel * float(params["signal_noise_ratio"])
        if intensity >= featureIntensityThreshold:
            pctTF = (maxRt - minRt)  * 100
            # Organize features in a structured numpy array form
            dataType = "f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8"
            if n == 0:
                res = np.array([(mz, intensity, z, rt, minRt, maxRt, ms1, minMs1, maxMs1, snRatio, pctTF)], dtype=dataType)
                n += 1
            else:
                res = np.append(res, np.array([(mz, intensity, z, rt, minRt, maxRt, ms1, minMs1, maxMs1, snRatio, pctTF)], dtype=res.dtype))
        else:
            continue
    res.dtype.names = ("mz", "intensity", "z", "RT", "minRT", "maxRT", "MS1", "minMS1", "maxMS1", "SNratio", "PercentageTF")
    res = pd.DataFrame(res) # Convert the features to a pandas DataFrame

    return res


def detectFeatures(inputFile, paramFile):
    print("  Feature detection")
    now = datetime.now()
    nowString = now.strftime("%Y/%m/%d %H:%M:%S")
    print("  " + nowString)
    reader = mzxml.read(inputFile)

    # Parameters and initialization
    params = getParams(paramFile)
    gap = int(params["skipping_scans"])

    # Obtain the array of MS1 spectra
    ms1 = getMs1(reader, params)
    nMs1 = len(ms1)

    #####################
    # Feature detection #
    #####################
    f = []
    cache = []
    noise = {}
    progress = progressBar(nMs1)
    for i in range(nMs1):
        progress.increment()
        # Append spectra to "cache" array considering the "gap"
        # When the i-th (MS1) spectrum is considered with the parameter, gap = 1,
        # "cache" array will contain
        # [(i-2)-th spectrum
        #  (i-1)-th spectrum
        #    (i)-th spectrum
        #  (i+1)-th spectrum
        #  (i+2)-th spectrum]
        minInd = max(0, i - gap - 1)
        maxInd = min(nMs1 - 1, i + gap + 1)
        if i == 0:
            for j in range(maxInd + 1):
                spec = detectPeaks(ms1[j], params)
                spec["index"] = j
                cache.append(spec)
        else:
            for j in range(oldMinInd, minInd):
                cache.pop(0)  # Remove the first element in cache
            for j in range(oldMaxInd + 1, maxInd + 1):
                spec = detectPeaks(ms1[j], params)
                spec["index"] = j
                cache.append(spec)

        ##################
        # Peak reduction #
        ##################
        # In each MS1 spectrum, each peak is examined whether it is reproducible in adjacent MS1 scans (considering "gap")
        # In addition to the reduction, the noise level in each MS1 scan is also estimated
        cache, noise = reducePeaks(i, minInd, maxInd, cache, noise, params)

        ##########################################
        # Peak merging (i.e. feature generation) #
        ##########################################
        # Using the retained MS1 peaks, pre-features are formed (i.e., a set of MS1 peaks over multiple MS1 scans at a specific m/z)
        f, cache = mergePeaks(i, minInd, cache, f, params)

        oldMinInd = minInd
        oldMaxInd = maxInd

    #########################
    # Filtering of features #
    #########################
    f = [feature for feature in f if feature is not None]   # Remove empty features
    f = filterFeatures(f)

    #######################
    # Output organization #
    #######################
    f = organizeFeatures(f, noise, params)  # Output format of features (using pandas DataFrame)

    return f


if __name__ == "__main__":
    # Initialization
    paramFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Dev\OOP\Dataset\jumpm_positive.params"
    inputFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Dev\OOP\Dataset\IROA_c18_target1.mzXML"
    features = detectFeatures(inputFile, paramFile)
    print()

