import os, sys, logging, numpy as np, pandas as pd
from utils import *
from datetime import datetime
from pyteomics import mzxml
from featureDetection import *


def defineFeature(uid, monoMzs, tol, idx, spec):
    df = spec.loc[idx]
    mzArray, intensityArray, areaArray, rtArray = [], [], [], []
    for monoMz in monoMzs:
        lL = monoMz - monoMz * tol / 1e6
        uL = monoMz + monoMz * tol / 1e6
        # For each m/z value of isotopologues,
        # 1. Extract a (sub) dataframe containing peaks around the "mz"
        # 2. There may be multiple peaks within the tolerance in one scan. So the dataframe is grouped by the "intensity"
        #    and the strongest (maximum) peak is chosen
        subDf = df.loc[(df["m/z array"] >= lL) & (df["m/z array"] <= uL)]
        subDf = subDf.sort_values("intensity array", ascending=False).groupby("num", as_index=False).first()

        mz = sum(subDf["m/z array"] * subDf["intensity array"]) / sum(subDf["intensity array"])    # Weighted averaged m/z
        intensity = max(subDf["intensity array"])
        area = sum(subDf["intensity array"])
        rt = sum(subDf["retentionTime"] * subDf["intensity array"]) / sum(subDf["intensity array"])    # Weighted averaged RT
        mzArray.append(mz)
        intensityArray.append(intensity)
        areaArray.append(area)
        rtArray.append(rt)
    z = 1   # Assume that z is always 1 (as of now)
    rt = sum(rtArray) / len(rtArray)
    minRt, maxRt = min(df["retentionTime"]), max(df["retentionTime"])
    # ms1 = df["num"].astype(int).iloc[idx]
    minMs1, maxMs1 = min(df["num"].astype(int)), max(df["num"].astype(int))
    # noiseLevel = np.percentile(spec[spec["num"] == str(ms1)]["intensity array"], 25)
    # snRatio = intensity / noiseLevel
    gMaxRt, gMinRt = max(spec["retentionTime"]), min(spec["retentionTime"])
    pctTF = (maxRt - minRt) / (gMaxRt - gMinRt) * 100
    metabolite = uid # HMDB ID
    res = {"mz": mzArray,
           "intensity": intensityArray,
           "area": areaArray,
           "z": z,
           "RT": rt,
           "minRT": minRt,
           "maxRT": maxRt,
           # "MS1": ms1,
           "minMS1": minMs1,
           "maxMS1": maxMs1,
           # "snRatio": snRatio,
           "pctTF": pctTF,
           "metabolite": metabolite}
    return res


def detectPeaks(spec):
    # m/z and intensity arrays from a spectrum object
    mzArray = spec["m/z array"]
    intensityArray = spec["intensity array"]
    nPeaks = len(mzArray)
    newMzArray = np.array([])
    newIntensityArray = np.array([])

    # Detect peaks (i.e. centroidization of MS1 spectrum)
    for i in range(2, nPeaks - 2):
        if intensityArray[i] > 0:
            # Consider 2 points before and after the point of interest x, i.e. 5 point window
            b2, b1, x, a1, a2 = intensityArray[(i - 2):(i + 3)]
            if isMax(b2, b1, x, a1, a2):
                # If x is the local maximum in a 5-point window, lower and upper bounds for a peak will be explored
                # Refer Figure 1a and b in the paper, Cox and Mann, Nature Biotech. 2008; 26: 1367-22
                minInd = findMinPeakIndex(i, intensityArray)
                maxInd = findMaxPeakIndex(i, intensityArray)
                if (maxInd - minInd) > 2:
                    newMz, newIntensity = findPeakCenter(minInd, i, maxInd, mzArray, intensityArray)
                    newMzArray = np.append(newMzArray, newMz)
                    newIntensityArray = np.append(newIntensityArray, newIntensity)

    # Update "spec" object
    spec["m/z array"] = newMzArray
    spec["intensity array"] = newIntensityArray
    return spec


def getScanIndex(df, mz, tol):
    lL = mz - mz * tol / 1e6
    uL = mz + mz * tol / 1e6
    rows = (df["m/z array"] >= lL) & (df["m/z array"] <= uL)
    res = df.loc[rows].index
    return res

########################################
# DataFrame containing all MS1 spectra #
########################################

inputFile = "isotope_distribution.xlsx"
paramFile = "jumpm_targeted.params"
# mzxmlFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_QC3.mzXML"
mzxmlFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML"

# Input dataframe (from Surendhar's program)
# What is going to be a "key"? metabolite name? HMDB ID?
df = pd.read_excel(inputFile)
uids = df["idhmdb"].unique()

# Pyteomics reader to pandas DataFrame?
reader = mzxml.read(mzxmlFile)
spec = pd.DataFrame(reader)
spec.drop(columns=["peaksCount", "polarity", "scanType", "lowMz", "highMz",
                 "basePeakMz", "basePeakIntensity", "totIonCurrent", "id",
                 "collisionEnergy", "precursorMz"], inplace=True)
spec = spec[spec["msLevel"] == 1].reset_index(drop=True)
spec = spec.apply(pd.Series.explode)    # This operation extracts "m/z array" and "intensity array" into columns

# # Preliminary mass correction (so slow as of now)
# dfMzShift = pd.DataFrame()
# for uid in uids:
#     monoMzs = [float(mz.split(";")[0]) for mz in df[df["idhmdb"] == uid]["isotope_M/Z"]]
#     for monoMz in monoMzs:
#         lL = monoMz - 0.5   # Isolation window?
#         uL = monoMz + 0.5
#         rows = (spec["m/z array"] >= lL) & (spec["m/z array"] <= uL)
#         subDf = spec[rows].groupby("num").agg({"intensity array": "max", "m/z array": "first"})
#         subDf["m/z array"] = abs(subDf["m/z array"] - monoMz) / monoMz * 1e6    # PPM
#         dfMzShift = dfMzShift.append(subDf, ignore_index=True)

tol = 50
gap = 5     # Gap allowed in MS1 scan
n = 0       # Index of features
fArray = []
for uid in uids:    # For each given metabolite
    monoMzs = [float(mz.split(";")[0]) for mz in df[df["idhmdb"] == uid]["isotope_M/Z"]]
    # Let's consider MS1 scans where M0 and M1 coexist
    # Although M0, M1, ..., Mn should coexist in principle, some isotopologues are not easily detected

    scanIdx = []
    for monoMz in monoMzs:
        # lL = monoMz - monoMz * tol / 1e6
        # uL = monoMz + monoMz * tol / 1e6
        # rows = rows & ((spec["m/z array"] >= lL) & (spec["m/z array"] <= uL))
        if len(scanIdx) > 0:
            scanIdx = np.intersect1d(scanIdx, getScanIndex(spec, monoMz, tol))
        else:
            scanIdx = getScanIndex(spec, monoMz, tol)

    # Feature detection with considering the gap in MS1 scan
    if len(scanIdx) > 0:
        idxArray = []
        j = 0
        for i in range(1, len(scanIdx)):
            if scanIdx[i] - scanIdx[i - 1] > gap:
                if i - j > 1:  # Avoid "singleton" features
                    idxArray.append(scanIdx[j: i])
                j = i
            else:
                if i == len(scanIdx) - 1:  # Last element
                    idxArray.append(scanIdx[j: len(scanIdx)])
        for idx in idxArray:
            # Define feature(s)
            f = defineFeature(uid, monoMzs, tol, idx, spec)
            fArray.append(f)

    # if len(scanIdx) > 0:
    #     idxArray = []
    #     for i in range(1, len(scanIdx)):
    #         if scanIdx[i] - scanIdx[i - 1] > gap:
    #             if len(idxArray) == 0:
    #                 j = 0
    #             if i - j > 1:   # Avoid "singleton" features
    #                 idxArray.append(scanIdx[j: i])
    #             j = i
    #     for idx in idxArray:
    #         # Define feature(s)
    #         f = defineFeature(uid, idx, spec)
    #         fArray.append(f)


    # # Feature detection with considering the gap in MS1 scan
    # if sum(rows) > 0:
    #     idxArray = []
    #     for i in range(1, len(spec[rows].index)):
    #         if spec[rows].index[i] - spec[rows].index[i - 1] > gap:
    #             if len(idxArray) == 0:
    #                 j = 0
    #             if i - j > 1:   # Avoid "singleton" features
    #                 idxArray.append(spec[rows].index[j: i])
    #             j = i
    #     for idx in idxArray:
    #         # Define features
    #         f = defineFeature(uid, rows, idx, spec)
    #         fArray.append(f)
    print()





# inputFile = r"../../Datasets/hilic_neg/neg_ctrl3.mzXML"
# inputFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_ctrl3.mzXML"

# Pyteomics reader to pandas DataFrame?
# reader = mzxml.read(inputFile)
# df = pd.DataFrame(reader)
# df.drop(columns=["peaksCount", "polarity", "scanType", "lowMz", "highMz",
#                  "basePeakMz", "basePeakIntensity", "totIonCurrent", "id",
#                  "collisionEnergy", "precursorMz"], inplace=True)
# df = df[df["msLevel"] == 1]
# df = df.apply(pd.Series.explode)    # This operation extracts "m/z array" and "intensity array" into columns
# print()

# t = datetime.now()
# from featureDetection import detectFeatures
# df = detectFeatures(inputFile, "jumpm_targeted.params")
# print((datetime.now() - t).total_seconds())


# paramFile = "jumpm_targeted.params"
# params = getParams(paramFile)
# print(params)


###################################################
# From features detected by "featureDetection.py" #
###################################################
# import pickle
#
# inputFile = "isotope_distribution.xlsx"
# paramFile = "jumpm_targeted.params"
# mzxmlFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML"
#
# df = pd.read_excel(inputFile)
# features = detectFeatures(df, mzxmlFile, paramFile)
# # features = pickle.load(open("features.pickle", "rb"))
#
# # For each unique ID (as of now, HMDB ID),
# # 1. Look for the feature corresponding to M0 (monoisotopic peak of the metabolite)
# # 2. Look for the features corresponding its isotopologues (i.e., M1, M2, ..., Mn)
# tol = 50
# uids = df["idhmdb"].unique()
# for uid in uids:
#     mzs = [float(mz.split(";")[0]) for mz in df[df["idhmdb"] == uid]["isotope_M/Z"]]
#     idx = features["mz"] < 0   # Initial index (all False)
#     for mz in mzs:
#         lL = mz - mz * tol / 1e6
#         uL = mz + mz * tol / 1e6
#         idx = idx & ((features["mz"] >= lL) & (features["mz"] <= uL))
#     if sum(idx) > 0:
#         print()

