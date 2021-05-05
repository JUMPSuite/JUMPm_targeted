#!/usr/bin/python

from pyteomics import mzxml
from utils import *


def groupMzValues(mz, ppm):
    mzdiff = np.diff(mz)
    if ppm > 0:
        res = mzdiff >= (mz[:-1] * ppm / 1e6)
    res = np.insert(res, 0, 0)
    res = np.cumsum(res)
    return res


# def estimateMzScattering(x):
#     # Kernel density estimation of "mzdiff" using Silverman's rule-of-thumb
#     sigma = np.std(x)
#     q75, q25 = np.percentile(x, [75 ,25])
#     iqr = q75 - q25
#     n = len(x)
#     bw = 0.9 * min(sigma, iqr / 1.34) * (n ** (-1 / 5)) # Bandwidth for KDE using rule-of-thumb
#     kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(x.reshape(-1, 1))
#     xScores = np.linspace(min(x) - 3 * bw, max(x) + 3 * bw, 512)
#     scores = np.exp(kde.score_samples(xScores.reshape(-1, 1)))
#     idx = np.where(np.diff(np.sign(np.diff(scores))) == 2)[0] + 1
#     if len(idx) > 1:
#         idx = idx[0]
#     return xScores[idx].item()


def mergeMs2(mzs, ints, mzGroups):
    mzArray, intArray = [], []
    for i in np.unique(mzGroups):
        idx = np.where(mzGroups == i)[0]
        mz = np.average(mzs[idx], weights=ints[idx])
        intensity = sum(ints[idx])
        mzArray.append(mz)
        intArray.append(intensity)
    spec = {"mz": np.array(mzArray), "intensity": np.array(intArray)}
    return spec


def intraConsolidation(ms2, tol):
    # input arguments
    # ms2: dictionary of MS2 spectra; key = scan number, val = {"mz": np.array, "intensity": np.array}
    # tol: m/z tolerance for merging MS2 peaks
    mzs = ms2["mz"]
    ints = ms2["intensity"]
    # Sort m/z values (and intensities accordingly)
    idx = np.argsort(mzs)
    mzs = mzs[idx]
    ints = ints[idx]
    mzGroups = groupMzValues(mzs, tol)
    spec = mergeMs2(mzs, ints, mzGroups)
    spec = simplifyMs2(spec)
    return spec


# def interConsolidation(specs, tol):
#     # input arguments
#     # specs: array of "feature-to-spectrum" [# features x # runs]
#     #        specs[i, j] = (intra-consolidated) MS2 spectrum of j-th run corresponding to i-th feature
#     # tol: m/z tolerance for merging MS2 peaks
#     specs = [i for i in specs if i is not None]  # Skip "None"
#     if len(specs) > 1:
#         mzs = np.array([])
#         ints = np.array([])
#         # Extract m/z (and intensity) values from individual MS2 spectrum and merge them to "mzs" (and "ints")
#         for s in specs:
#             mzs = np.append(mzs, s["mz"])
#             ints = np.append(ints, s["intensity"])
#         # Sort m/z values (and intensities accordingly)
#         idx = np.argsort(mzs)
#         mzs = mzs[idx]
#         ints = ints[idx]
#         mzGroups = groupMzValues(mzs, tol)
#         spec = mergeMs2(mzs, ints, mzGroups)
#     else:
#         spec = specs[0]
#     spec = simplifyMs2(spec)
#     return spec


def simplifyMs2(spec):
    # Simplification of the merged spectrum
    # 1. Limit the number of peaks in each .dta file to 100
    # 2. Divide m/z-range into 10 bins (e.g. 0~100, 100~200, etc.) and retain the 10 largest peaks in each bin
    if len(spec["mz"]) > 100:
        nBins = 10
        bins = np.linspace(min(spec["mz"]), max(spec["mz"]), (nBins + 1))
        filteredSpec = {"mz": [], "intensity": []}
        for i in range(nBins):
            ind = np.where((spec["mz"] >= bins[i]) & (spec["mz"] < bins[i + 1]))[0]
            if len(ind) > 10:
                # Select 10 highest intensity peaks
                ind10 = sorted(range(len(spec["intensity"][ind])), key=lambda j: spec["intensity"][ind][j],
                               reverse=True)[:10]
                ind = ind[ind10]
            filteredSpec["mz"] = np.append(filteredSpec["mz"], spec["mz"][ind])
            filteredSpec["intensity"] = np.append(filteredSpec["intensity"], spec["intensity"][ind])
        spec = filteredSpec
    # Sort the spectrum in ascending order of m/z
    ind = np.argsort(spec["mz"])
    spec["mz"] = spec["mz"][ind]
    spec["intensity"] = spec["intensity"][ind]
    return spec


def ms2ForFeatures(df, mzxmlFile, paramFile):
    print("  Identification of MS2 spectra for the features")
    print("  ==============================================")
    df["MS2"] = ""

    ######################################
    # Load parameters and initialization #
    ######################################
    params = getParams(paramFile)
    tolIsolation = float(params["isolation_window"]) / 2
    tolIntraMS2Consolidation = float(params["tol_intra_ms2_consolidation"])

    #################################################
    # Assignment of MS2 spectra to features         #
    # Consolidation of MS2 spectra for each feature #
    #################################################
    reader = mzxml.MzXML(mzxmlFile)
    print("  %s is being processed" % os.path.basename(mzxmlFile))
    print("  Looking for MS2 scan(s) responsible for each feature")
    ms2Array = []
    for i in range(df.shape[0]):
        minScan, maxScan = int(df.iloc[i]["minMS1"]), int(df.iloc[i]["maxMS1"])
        ms2Dict = {"mz": np.array([]), "intensity": np.array([])}
        for j in range(minScan, maxScan + 100):
            # progress.increment()
            spec = reader[str(j)]
            msLevel = spec["msLevel"]
            if msLevel == 1 & j > maxScan:
                break
            elif msLevel == 2:
                # Find MS2 scans which satisfy the following conditions

                # From the discussion around June 2020,
                # 1. In ReAdW-derived mzXML files, precursor m/z values are in two tags: "precursorMz" and "filterLine"
                # 2. Through Haiyan's manual inspection, the real precursor m/z value is closer to one in "filterLine" tag
                # 3. So, in this script, precursor m/z of MS2 scan is obtained from "filterLine" tag
                # 4. Note that it may be specific to ReAdW-derived mzXML files since MSConvert-derived mzXML files do not have "filterLine" tag
                # 4.1. In this case, maybe the use of mzML (instead of mzXML) would be a solution (to-do later)
                # precMz = spec["precursorMz"][0]["precursorMz"]  # Precursor m/z from "precursorMz" tag
                p = re.search("([0-9.]+)\\@", spec["filterLine"])
                precMz = float(p.group(1))
                if precMz - tolIsolation <= df.iloc[i]["mz"] <= precMz + tolIsolation:
                    ms2Dict["mz"] = np.append(ms2Dict["mz"], spec["m/z array"])
                    ms2Dict["intensity"] = np.append(ms2Dict["intensity"], spec["intensity array"])
        ms2Array.append(ms2Dict)

    print("  Merging MS2 spectra for each feature within a run (it may take a while)")
    for i in range(len(ms2Array)):
        ms2Array[i] = intraConsolidation(ms2Array[i], tolIntraMS2Consolidation)
    df["MS2"] = ms2Array

    return df
