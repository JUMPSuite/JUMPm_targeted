import re, sys, os, pickle, numpy as np, pandas as pd
from pyteomics import mass


def getParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line)  # Remove all whitespaces

            # Exception for "feature_files" parameter
            if "feature_files" in parameters and line.endswith(".feature"):
                parameters["feature_files"].append(line)
            elif "library" in parameters and line.endswith(".db"):
                parameters["library"].append(line)
            else:
                key = line.split('=')[0]
                val = line.split('=')[1]
                if key == "feature_files" or key == "library":
                    parameters[key] = [val]
                else:
                    parameters[key] = val
    return parameters


def getMs1(reader, params):
    try:
        firstScan = int(params["first_scan_extraction"])
    except:
        firstScan = 1
    try:
        lastScan = int(params["last_scan_extraction"])
    except:
        lastScan = 1000000
    ms1 = []
    n = 0

    # When targeted metabolites are given, minimum and maximum m/z to be considered are set
    # Note that the charge state is always assumed to be 1 (i.e., z = 1)
    proton = 1.007276466812
    isTargeted = 0
    minMz, maxMz = 1e4, 0
    for k, v in params.items():
        if k.startswith("Metabolite"):
            isTargeted = 1
            mz = mass.calculate_mass(formula=v) + proton    # Assume "positive" mode
            if mz < minMz:
                minMz = mz
            if mz > maxMz:
                maxMz = mz

    # To find isotopologues, minMz and maxMz should be a bit expanded
    minMz -= 5
    maxMz += 10

    with reader:
        for spec in reader:
            msLevel = spec["msLevel"]  # int type
            scanNum = spec["num"]  # str type
            if msLevel == 1 and firstScan <= int(scanNum) <= lastScan:
                spec["scanIndex"] = n

                # Reduce MS1 peaks when targeted
                if isTargeted == 1:
                    ind = (spec["m/z array"] >= minMz) & (spec["m/z array"] <= maxMz)
                    spec["m/z array"] = spec["m/z array"][ind]
                    spec["intensity array"] = spec["intensity array"][ind]

                # # Mass (m/z) shift according to the given parameter (unit of ppm)
                # if params["mass_shift"] is not None:
                #     spec["m/z array"] = spec["m/z array"] / (1 + float(params["mass_shift"]) / 1e6)
                ms1.append(spec)
                n += 1
            elif int(scanNum) > lastScan:
                break
    return ms1


def calcMS2Similarity(featSpec, libSpec):
    # Calculation of MS2 similarity between a feature and a library compound
    # Reference: Clustering millions of tandem mass spectra, J Proteome Res. 2008; 7: 113-22

    # Input arguments
    # featSpec (dictionary): MS2 spectrum of a feature (key = "mz", "intensity")
    # libSpec (dictionary): MS2 spectrum of a library compound (key = "mz", "intensity", "index" (ignorable))
    nPeaks = 30 # Default = 30 according to the above reference
    k = min(nPeaks, min(len(featSpec["mz"]), len(libSpec["mz"])))

    # Keep $k strongest peaks in both spectra
    # featDict[mz] = intensity
    # libDict[mz] = intensity
    featDict, libDict = {}, {}
    ind = np.argsort([-i for i in featSpec["intensity"]])
    for i in ind[0:k]:
        featDict[featSpec["mz"][i]] = featSpec["intensity"][i]
    ind = np.argsort([-i for i in libSpec["intensity"]])
    for i in ind[0:k]:
        libDict[libSpec["mz"][i]] = libSpec["intensity"][i]

    # Join two sets of m/z values and make a new set of unique m/z values
    # Duplicate masses are removed as follows
    # - We consider two peaks to have a similar mass if they are within 0.5 Da from each other)
    # - For those two peaks having similar mass, the lower one will be the unique one
    #   (e.g. One peak with m/z = 100 and the other peak with m/z = 100.4 -> they will be merged to m/z = 100)
    mzArray = list(featDict.keys()) + list(libDict.keys())
    mzArray = sorted(mzArray)
    mzDict = {}
    val = 0
    for mz in mzArray:
        if abs(mz - val) <= 0.5:
            mzDict[mz] = val
        else:
            mzDict[mz] = mz
            val = mz

    # Reduction of spectrum to a vector by assigning to each intensity to the unique m/z bins
    # And then, calculate the similarity; normalized dot-product
    s = {}
    for key, val in mzDict.items():
        s[val] = {}
        s[val]["feat"] = 0
        s[val]["lib"] = 0

    for key, val in mzDict.items():
        if key in featDict:
            s[val]["feat"] += np.sqrt(featDict[key])
        if key in libDict:
            s[val]["lib"] += np.sqrt(libDict[key])

    num, den1, den2 = 0, 0, 0    # numerator, denominator1, denominator2
    for mz in s.keys():
        num += s[mz]["feat"] * s[mz]["lib"]
        den1 += s[mz]["feat"] ** 2
        den2 += s[mz]["lib"] ** 2

    if den1 * den2 == 0:
        normDotProduct = 0
    else:
        normDotProduct = num / np.sqrt(den1 * den2)

    return normDotProduct


class progressBar:
    def __init__(self, total):
        self.total = total
        self.barLength = 20
        self.count = 0
        self.progress = 0
        self.block = 0
        self.status = ""

    def increment(self, nIncrement=None):
        if nIncrement == None:
            self.count += 1
        else:
            self.count = nIncrement
        self.progress = self.count / self.total
        self.block = int(round(self.barLength * self.progress))
        if self.progress == 1:
            self.status = "Done...\r\n"
        else:
            self.status = ""
        #         self.status = str(self.count) + "/" + str(self.total)
        text = "\r  Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block),
                                                     int(self.progress * 100), self.status)
        sys.stdout.write(text)
        sys.stdout.flush()
