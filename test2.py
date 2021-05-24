import os, sys, logging, numpy as np, pandas as pd
from utils import *
from datetime import datetime
from pyteomics import mzxml
from featureDetection import *

# scanIdx = [64, 411, 1311, 1356, 1357, 1359, 1544, 1547, 1548]
# scanIdx = [1, 2, 3, 7, 20, 27, 69, 100, 103, 104, 105, 107, 338, 400]
# idxArray = []
# gap = 5
# j = 0
# for i in range(1, len(scanIdx)):
#     if scanIdx[i] - scanIdx[i - 1] > gap:
#         if i - j > 1:  # Avoid "singleton" features
#             idxArray.append(scanIdx[j: i])
#         j = i
#     else:
#         if i == len(scanIdx) - 1:   # Last element
#             idxArray.append(scanIdx[j: len(scanIdx)])
#
# print(idxArray)


def checkPeakMz(val, ref, tol):
    lL = ref - ref * tol / 1e6
    uL = ref + ref * tol / 1e6
    if lL <= val <= uL:
        return 1
    else:
        return -1


def findRtOverlap(min1, max1, min2, max2):
    if min1 < min2:
        # Examples
        # |----------------|    range of input1
        #       |----|          range of input2
        #          |----------| range of input2
        if max1 < max2:
            overlap = max1 - min2
        else:
            overlap = max2 - min2
    else:
        # Examples
        #      |--------|       range of input1
        # |----------|          range of input2
        # |-----------------|   range of input2
        if max1 < max2:
            overlap = max1 - min1
        else:
            overlap = max2 - min1

    return overlap



import pickle
[f, df] = pickle.load(open("tmpFeatures.pickle", "rb"))
tol = 50
f = f.sort_values(["metabolite", "mz"], ignore_index=True)  # Sort by HMDB ID and m/z value
uids = sorted(df["idhmdb"].unique())
fArray = []
for uid in uids:
    df1 = f[f["metabolite"] == uid] # DataFrame for the specific "uid
    idxArray = []
    used = []
    for i in range(df1.shape[0]):
        if i in used:
            continue
        idx = []
        for j in range(i + 1, df1.shape[0]):
            if j in used:
                continue
            res = findRtOverlap(df1.iloc[i]["minRT"], df1.iloc[i]["maxRT"], df1.iloc[j]["minRT"], df1.iloc[j]["maxRT"])
            if 0 < res < 600:    # Too much RT-overlap (> 10 minutes) is not meaningful
                idx.append(j)
                used.append(j)
        if len(idx) > 0:
            idx.append(i)
            idxArray.append(idx)
            used.append(i)

    # Arrange the feature representing isotopologues
    mzs = [float(mz.split(";")[0]) for mz in df[df["idhmdb"] == uid]["isotope_M/Z"]]
    for i in range(len(idxArray)):
        mzArray, intensityArray, rtArray = [], [], []
        df2 = df1.iloc[idxArray[i]].sort_values("mz")
        for j in range(len(mzs)):
            mz = mzs[j]
            row = df2["mz"].astype(int) == int(mz)
            if j == 0 and sum(row) == 0:    # If there's no feature corresponding to "monoisotopic" peak, then break (move to the next)
                break
            if sum(row) > 0:
                mzArray.append(df2[row]["mz"].values[0])
                intensityArray.append(df2[row]["intensity"].values[0])
                if j == 0:
                    rt = df2[row]["RT"].values[0]
                    minRt = df2[row]["minRT"].values[0]
                    maxRt = df2[row]["maxRT"].values[0]
                    minMs1 = df2[row]["minMS1"].values[0]
                    maxMs1 = df2[row]["maxMS1"].values[0]
            else:
                mzArray.append(0)
                intensityArray.append(0)
        if len(mzArray) > 0:
            newF = {"mz": mzArray, "intensity": intensityArray, "z": 1,
                    "RT": rt, "minRT": minRt, "maxRT": maxRt,
                    "minMS1": minMs1, "maxMS1": maxMs1, "metabolite": uid}
            fArray.append(newF)

print()
