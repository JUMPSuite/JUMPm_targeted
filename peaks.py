import sys, numpy as np


class feature:
    def __init__(self, mz=None, intensity=None, rt=None, minRt=None, maxRt=None, ms1=None, minMs1=None, maxMs1=None, index=None):
        self.mz = mz
        self.intensity = intensity
        self.rt = rt
        self.minRt = minRt
        self.maxRt = maxRt
        self.ms1 = ms1
        self.minMs1 = minMs1
        self.maxMs1 = maxMs1
        self.index = index


def detectPeaks(spec, params):
    # Parameters
    if params["data_acquisition_mode"] == "1":
        isCentroided = 1
    elif params["data_acquisition_mode"] == "2":
        isCentroided = 0
    else:
        print("Please set the proper 'data_acquisition_mode' parameter")
        sys.exit("")
    try:
        intensityThreshold = float(params["min_peak_intensity"])
    except:
        intensityThreshold = 0

    # m/z and intensity arrays from a spectrum object
    mzArray = spec["m/z array"]
    intensityArray = spec["intensity array"]
    # Check whether there's any "zero" intensity or not
    # If there's any "zero" intensity, it indicates that the spectrum is NOT centroided
    idxZero = sum(intensityArray == 0)
    if idxZero > 0:
        isCentroided = 0
    nPeaks = len(mzArray)
    newMzArray = np.array([])
    newIntensityArray = np.array([])

    # Detect peaks (i.e. centroidization of MS1 spectrum)
    if isCentroided == 0:  # i.e. Profile mode MS1
        for i in range(2, nPeaks - 2):
            if intensityArray[i] > 0:
                # Consider 2 points before and after the point of interest x, i.e. 5 point window
                b2, b1, x, a1, a2 = intensityArray[(i - 2):(i + 3)]
                if x >= intensityThreshold:
                    if isMax(b2, b1, x, a1, a2):
                        # If x is the local maximum in a 5-point window, lower and upper bounds for a peak will be explored
                        # Refer Figure 1a and b in the paper, Cox and Mann, Nature Biotech. 2008; 26: 1367-22
                        minInd = findMinPeakIndex(i, intensityArray)
                        maxInd = findMaxPeakIndex(i, intensityArray)
                        if (maxInd - minInd) > 2:
                            newMz, newIntensity = findPeakCenter(minInd, i, maxInd, mzArray, intensityArray)
                            newMzArray = np.append(newMzArray, newMz)
                            newIntensityArray = np.append(newIntensityArray, newIntensity)
    else:  # i.e. Centroid mode MS1
        idx = np.where(intensityArray >= intensityThreshold)[0]
        newMzArray = mzArray[idx]
        newIntensityArray = intensityArray[idx]

    # Update "spec" object
    spec["m/z array"] = newMzArray
    spec["intensity array"] = newIntensityArray
    return spec


def isMax(b2, b1, x, a1, a2):
    if x > b1 and x > a1:
        return True
    if x > b2 and x == b1 and x > a1:
        return True
    if x > b1 and x == a1 and x > a2:
        return True
    return False


def findMinPeakIndex(ind, array):
    while ind > 0 and array[ind] != 0 and array[ind - 1] <= array[ind]:
        ind -= 1
    return ind + 1


def findMaxPeakIndex(ind, array):
    count = len(array)
    while ind < count and array[ind] != 0 and array[ind + 1] <= array[ind]:
        ind += 1
    return ind - 1


def findPeakCenter(minInd, centerInd, maxInd, mz, intensity):
    # Find the center of a peak composed of five data points
    centerMz = 0
    centerIntensity = 0
    for i in range(minInd, maxInd + 1):
        if intensity[i] >= centerIntensity:
            centerIntensity = intensity[i]  # Take the maximum intensity as a peak intensity

    # There"s a plateau, bu others are zeros
    if minInd == maxInd:
        centerMz = mz[maxInd]
        return centerMz, centerIntensity

    # Right-angled triangle-shaped peak
    if minInd == centerInd:
        centerMz = estimate2(mz[centerInd], mz[centerInd + 1], intensity[centerInd], intensity[centerInd + 1])
        return centerMz, centerIntensity

    # Left-angled triangle-shaped peak
    if maxInd == centerInd:
        centerMz = estimate2(mz[centerInd - 1], mz[centerInd], intensity[centerInd - 1], intensity[centerInd])
        return centerMz, centerIntensity

    # Typical bell(triangle)-shaped peak
    centerMz = estimate3(mz[centerInd - 1], mz[centerInd], mz[centerInd + 1], intensity[centerInd - 1],
                         intensity[centerInd], intensity[centerInd + 1])
    return centerMz, centerIntensity


def estimate2(m1, m2, i1, i2):
    centerVal = (m1 * i1 + m2 * i2) / (i1 + i2)  # Intensity-weighted average of m/z
    return centerVal


def estimate3(m1, m2, m3, i1, i2, i3):
    l1 = np.log(i1)
    l2 = np.log(i2)
    l3 = np.log(i3)
    centerVal = 0.5 * ((l1 - l2) * (m3 ** 2 - m1 ** 2) - (l1 - l3) * (m2 ** 2 - m1 ** 2)) / (
            (l1 - l2) * (m3 - m1) - (l1 - l3) * (m2 - m1))
    return centerVal


def getClosest(spec, mz, tol):
    ind = np.argmin(abs(spec["m/z array"] - mz))
    diff = abs(mz - spec["m/z array"][ind]) / mz * 1e6
    if diff < tol:
        return 1, ind
    else:
        return 0, None


def searchNeighbors(currInd, minInd, maxInd, specArray, mz, params, direction):
    if direction == "backward":
        startInd = currInd - 1
        endInd = minInd - 1
        increment = -1
    elif direction == "forward":
        startInd = currInd + 1
        endInd = maxInd + 1
        increment = 1
    tol = float(params["mz_tolerance_peak_matching"])
    scanWindow = int(params["skipping_scans"]) + 1
    match, nTry, ind = 0, 0, None
    for i in range(startInd, endInd, increment):
        neighbor = specArray[i - minInd]
        if neighbor["m/z array"].size == 0:
            continue
        else:
            match, ind = getClosest(neighbor, mz, tol)
        if match == 1:
            break
        nTry += 1
        if nTry > scanWindow:
            break
    return match, ind


def reducePeaks(currInd, minInd, maxInd, cache, noise, params):
    # Input arguments
    #   currInd = (MS1 scan) index of the current spectrum
    #   minInd, maxInd = minimum/maximum index of MS1 scan in the "cache" array
    #   cache = temporary array of MS1 scans
    #   noise = dictionary of the noise level (key = MS1 scan number, val = noise intensity level)
    #   params = dictionary of parameters
    currScan = cache[currInd - minInd]
    nPeaks = len(currScan["m/z array"])
    valid = np.array([])
    for i in range(nPeaks):
        mz = currScan["m/z array"][i]
        # Backward search first
        match, ind = searchNeighbors(currInd, minInd, maxInd, cache, mz, params, "backward")
        if match != 1:
            # If any peak close to "mz" is not found, then perform forward search
            match, ind = searchNeighbors(currInd, minInd, maxInd, cache, mz, params, "forward")
        if match == 1:
            valid = np.append(valid, i)  # Index of peaks to be retained in the current MS1 scan
    valid = valid.astype(int)

    # Noise estimation
    if len(valid) > 0:
        noisePeaks = np.setdiff1d(currScan["intensity array"], currScan["intensity array"][valid])
        if len(noisePeaks) > 0:
            noiseLevel = np.percentile(noisePeaks, 25)
        else:
            noiseLevel = 1  # Nominal noise level
    else:
        noiseLevel = 500  # Hard-coded default noise level = 500
    noise[currScan["num"]] = noiseLevel

    # Peak reduction
    currScan["m/z array"] = currScan["m/z array"][valid]
    currScan["intensity array"] = currScan["intensity array"][valid]
    return cache, noise


def extendFeature(featureIndexArray, features, cache, currInd, minInd, peakIndex):
    currScan = cache[currInd - minInd]
    featureIndexArray = list(set(featureIndexArray))  # Make the index unique
    # Case 1. When the current peak (with a given m/z) corresponds to multiple features
    if len(featureIndexArray) > 1:
        repInd = min(featureIndexArray)  # Set the representative index as the minimum of "featureInd"
        for f in featureIndexArray:
            # Merge to the lowest indexed feature and remove the "merged" features
            # e.g., Assume that the current peak corresponds well to features[10], features[11] and features[12]
            #       - features[11] and features[12] will be merged to features[10] (the lowest indexed feature)
            #         features[10].mz will be extended to include features[11].mz and features[12].mz (same for other attributes like intensity, etc.)
            #       - (temporary) array of spectra should be revised; all peaks previously indicating features[11] and features[12] should indicate features[10]
            #       - Finally, features[11] and features[12] will be changed to None (to avoid any confusion, the size of feature array will be kept)
            if f == repInd:
                continue

            # Extend the attributes of the representative feature
            features[repInd].mz.extend(features[f].mz)
            features[repInd].intensity.extend(features[f].intensity)
            features[repInd].ms1.extend(features[f].ms1)
            features[repInd].rt.extend(features[f].rt)
            features[repInd].index.extend(features[f].index)

            # Revise the array of MS1 spectra (i.e., cache array)
            for s in features[f].index:
                for t in range(len(cache)):
                    if cache[t]["index"] == s:
                        for u in range(len(cache[t]["featureIndex"])):
                            if cache[t]["featureIndex"][u] == f:
                                cache[t]["featureIndex"][u] = repInd
            features[f] = None  # Keep the size of feature array
    # Case 2. When the current peak corresponds to only one feature
    else:
        repInd = featureIndexArray[0]

    # Update the information of cache array
    if "featureIndex" in currScan:
        currScan["featureIndex"].append(repInd)
    else:
        currScan["featureIndex"] = [repInd]

    # Update the feature corresponding to the current peak by adding the information of the peak
    features[repInd].mz.append(currScan["m/z array"][peakIndex])
    features[repInd].intensity.append(currScan["intensity array"][peakIndex])
    features[repInd].ms1.append(int(currScan["num"]))
    features[repInd].rt.append(currScan["retentionTime"])
    features[repInd].index.append(currScan["index"])

    return features, cache


def generateFeature(features, cache, currInd, minInd, peakIndex):
    currScan = cache[currInd - minInd]
    nFeatures = len(features)
    if "featureIndex" in currScan:
        currScan["featureIndex"].append(nFeatures)
    else:
        currScan["featureIndex"] = [nFeatures]
    # Since there's no peak in adjacent MS1 scan(s), a new feature is defined
    f = feature()
    f.mz = [currScan["m/z array"][peakIndex]]
    f.intensity = [currScan["intensity array"][peakIndex]]
    f.rt = [currScan["retentionTime"]]
    f.ms1 = [int(currScan["num"])]
    f.index = [currInd]
    features.append(f)
    nFeatures += 1

    return features, cache


def mergePeaks(currInd, minInd, cache, fArray, params):
    # Input arguments
    #   currInd = (MS1 scan) index of the current spectrum
    #   minInd = minimum index of MS1 scan in the "cache" array
    #   cache = temporary array of MS1 scans
    #   fArray = array of features
    #   params = dictionary of parameters
    matchTol = float(params["mz_tolerance_peak_matching"])
    currScan = cache[currInd - minInd]  # Current MS1 scan having the "reduced" peaks
    nPeaks = len(currScan["m/z array"])
    for i in range(nPeaks):
        mz = currScan["m/z array"][i]
        flagExtendFeature = 0
        matchedFeatureIndex = []

        # Backward search (merging is performed by only backward search, i.e. previous scans)
        # This part is similar to "searchNeighbor" function, but there's a difference as follows,
        # - "searchNeighbor" function STOPS and RETURNS outputs when it finds a peak close to the given m/z in any neighbor MS1 scan
        # - The following part keeps looking for peaks close to the given m/z in neighbor MS1 scans
        for j in range(currInd - 1, minInd - 1, -1):
            prevScan = cache[j - minInd]
            if prevScan["m/z array"].size == 0:
                continue
            else:
                match, ind = getClosest(prevScan, mz, matchTol)
                # match = 1 means that the j-th peak in the (reduced) i-th scan can form a 3D-peak with
                # ind-th peak in the (reduced) previous scan, i.e., prevScan
                if match == 1:
                    matchedFeatureIndex.append(prevScan["featureIndex"][ind])
                    flagExtendFeature = 1

        if flagExtendFeature == 1:
            # Case 1. When there is/are peak(s) close to the given m/z in previous MS1 scan(s)
            #         In this case, the peak(s) in previous scans already formed a feature
            #         So, the current peak has to be added to the above feature (the feature will be "updated")
            fArray, cache = extendFeature(matchedFeatureIndex, fArray, cache, currInd, minInd, i)
        else:
            # Case 2. When there's no peak close to the given m/z in neighbor MS1 scans,
            #         a new feature needs to be generated from the current peak
            fArray, cache = generateFeature(fArray, cache, currInd, minInd, i)

    return fArray, cache
