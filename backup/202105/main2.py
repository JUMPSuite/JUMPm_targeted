import os, sys, logging, numpy as np, pandas as pd
from utils import *
from datetime import datetime
from pyteomics import mzxml
from featureDetection import *
import rpy2.robjects as ro
from rpy2.robjects.vectors import IntVector, FloatVector
ro.r['options'](warn=-1)


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
        if subDf.shape[0] > 0:
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


# # Too slow
# def mzCorrection(uids, inputDf, spec):
#     # Preliminary mass correction (so slow as of now)
#     t1 = datetime.now()
#     mzShifts = {}
#     for uid in uids:
#         monoMzs = [float(mz.split(";")[0]) for mz in inputDf[inputDf["idhmdb"] == uid]["isotope_M/Z"]]
#         monoMz = monoMzs[0] # Only consider the monoisotopic one
#         lL = monoMz - monoMz * 50 / 1e6  # +/- 50 ppm (pre-set)
#         uL = monoMz + monoMz * 50 / 1e6
#         rows = (spec["m/z array"] >= lL) & (spec["m/z array"] <= uL)
#         df = spec[rows].sort_values("intensity array", ascending=False).groupby("num", as_index=False).first()
#         df["m/z array"] = abs(df["m/z array"] - monoMz) / monoMz * 1e6  # PPM
#         mzShifts[monoMz] = df[["num", "m/z array"]]
#     t2 = (datetime.now() - t1).total_seconds()
#     print(t2)
#
#     return mzShifts


# # Too slow
# def mzCorrection(uids, inputDf, reader):
#     # Preliminary mass correction
#     refMzs = []
#     mzShifts = {}
#     for uid in uids:
#         monoMzs = [float(mz.split(";")[0]) for mz in inputDf[inputDf["idhmdb"] == uid]["isotope_M/Z"]]
#         monoMz = monoMzs[0]  # Only consider the monoisotopic one
#         refMzs.append(monoMz)
#         mzShifts[monoMz] = []
#
#     rt = []
#     n = 0
#     for spec in reader:
#
#         n += 1
#         if n % 100 == 0:
#             print(spec["num"])
#
#
#         if spec["msLevel"] != 1:
#             continue
#         rt.append(spec["retentionTime"])
#         for refMz in refMzs:
#             lL = refMz - refMz * 100 / 1e6  # Isolation window?
#             uL = refMz + refMz * 100 / 1e6
#             idx = (spec["m/z array"] >= lL) & (spec["m/z array"] <= uL)
#             if sum(idx) > 0:
#                 maxIdx = np.argmax(spec["intensity array"][idx])
#                 mz = spec["m/z array"][idx][maxIdx]
#                 mzShifts[refMz].append(abs(mz - refMz) / refMz * 1e6)
#             else:
#                 mzShifts[refMz].append(np.nan)
#     print()


def loess():
    rstring = """
    loess.as = function(x, y, degree = 1, criterion="aicc", family="gaussian", user.span=NULL, plot=FALSE, ...) {

        criterion <- match.arg(criterion)
        family <- match.arg(family)
        x <- as.matrix(x)

        if ((ncol(x) != 1) & (ncol(x) != 2)) stop("The predictor 'x' should be one or two dimensional!!")
        if (!is.numeric(x)) stop("argument 'x' must be numeric!")
        if (!is.numeric(y)) stop("argument 'y' must be numeric!")
        if (any(is.na(x))) stop("'x' contains missing values!")
        if (any(is.na(y))) stop("'y' contains missing values!")
        if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span))) 
            stop("argument 'user.span' must be a numerical number!")
        if(nrow(x) != length(y)) stop("'x' and 'y' have different lengths!")
        if(length(y) < 3) stop("not enough observations!")

        data.bind <- data.frame(x=x, y=y)
        if (ncol(x) == 1) {
            names(data.bind) <- c("x", "y")
        } else { names(data.bind) <- c("x1", "x2", "y") }

        opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){
            as.crit <- function (x) {
                span <- x$pars$span
                traceL <- x$trace.hat
                sigma2 <- sum(x$residuals^2 ) / (x$n-1)
                aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
                gcv <- x$n*sigma2 / (x$n-traceL)^2
                result <- list(span=span, aicc=aicc, gcv=gcv)
                return(result)
            }
            criterion <- match.arg(criterion)
            fn <- function(span) {
                mod <- update(model, span=span)
                as.crit(mod)[[criterion]]
            }
            result <- optimize(fn, span.range)
            return(list(span=result$minimum, criterion=result$objective))
        }

        control = loess.control(surface = "direct")
        if (ncol(x)==1) {
            if (is.null(user.span)) {
                fit0 <- loess(y ~ x, degree=degree, family=family, data=data.bind, control=control, ...)
                span1 <- opt.span(fit0, criterion=criterion)$span
            } else {
                span1 <- user.span
            }
            fit <- loess(y ~ x, degree=degree, span=span1, family=family, data=data.bind, control=control, ...)
        } else {
            if (is.null(user.span)) {
                fit0 <- loess(y ~ x1 + x2, degree=degree,family=family, data.bind, control=control, ...)
                span1 <- opt.span(fit0, criterion=criterion)$span
            } else {
                span1 <- user.span
            }
            fit <- loess(y ~ x1 + x2, degree=degree, span=span1, family=family, data=data.bind, control=control...)
        }
        return(fit)
    }
    """
    return ro.r(rstring)


def correctMz(input, refMz):
    rLoess = loess()
    rPredict = ro.r("predict")
    df = input.copy()

    # Pre-filtering of MS1 peaks based on an intensity threshold
    intensityThreshold = 10000
    df = df[df["intensity array"] > 10000]

    # Estimate the m/z-shifts around the "refMz" by investigating the peaks around "refMz" +/- 50 ppm
    lL = refMz - refMz * 50 / 1e6
    uL = refMz + refMz * 50 / 1e6
    rows = (df["m/z array"] >= lL) & (df["m/z array"] <= uL)
    # When there are multiple peaks within the tolerance in a MS1 scan, the strongest peak is considered
    subDf = df[rows].sort_values("intensity array", ascending=False).groupby("num", as_index=False).first()
    # LOESS modeling and m/z correction
    x = subDf["retentionTime"]
    y = (subDf["m/z array"] - refMz) / refMz * 1e6    # m/z-shifts measured in PPM
    mod = rLoess(FloatVector(x), FloatVector(y), 1, "aicc", "gaussian")
    df["m/z array"] = df["m/z array"] / (1 + np.array(rPredict(mod, FloatVector(df["retentionTime"]))) / 1e6)
    return df, mod



########################################
# DataFrame containing all MS1 spectra #
########################################

t = datetime.now()

inputFile = "isotope_distribution.xlsx"
paramFile = "jumpm_targeted.params"
# mzxmlFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_QC3.mzXML"
mzxmlFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML"

# Input dataframe (from Surendhar's program)
# What is going to be a "key"? metabolite name? HMDB ID?
inputDf = pd.read_excel(inputFile)
reader = mzxml.MzXML(mzxmlFile)

# Extract m/z values and unique IDs of given metabolites (from an excel file)
mzs = [float(mz.split(";")[0]) for mz in inputDf.groupby("idhmdb").min()["isotope_M/Z"]]
uids = inputDf.groupby("idhmdb", as_index=False).min()["idhmdb"]
gMinMz, gMaxMz = min(mzs) - 10, max(mzs) + 10    # minimum and maximum m/z values to be considered

# Sort m/z values and unique IDs accordingly
idx = np.argsort(mzs)
mzs = sorted(mzs)
uids = uids[idx]

# Extract MS1 spectra (m/z array and intensity array are reduced)
# and convert it to a pandas dataframe
print("  Extraction of MS1 spectra")
ms1 = []
for spec in reader:
    if spec["msLevel"] == 1:
        idx = (spec["m/z array"] >= gMinMz) & (spec["m/z array"] <= gMaxMz) & (spec["intensity array"] > 1000)
        spec["m/z array"] = spec["m/z array"][idx]
        spec["intensity array"] = spec["intensity array"][idx]
        ms1.append(spec)
df = pd.DataFrame(ms1)
df.drop(columns=["peaksCount", "polarity", "scanType", "lowMz", "highMz", "basePeakMz",
                 "basePeakIntensity", "totIonCurrent", "id", "msLevel", "filterLine"], inplace=True)
df = df.apply(pd.Series.explode)    # This operation extracts "m/z array" and "intensity array" into columns
print("  Done ...\n")

# Feature detection (with mass correction)
print("  Feature detection")
tol = 20
gap = 5     # Gap allowed in MS1 scan
n = 0       # Index of features
fArray = []
mods = {}
for uid in uids:    # For each given metabolite
    print("  Features for {} is being identified".format(uid))
    # m/z values of isotopologues of the metabolite
    mzs = [float(mz.split(";")[0]) for mz in inputDf[inputDf["idhmdb"] == uid]["isotope_M/Z"]]
    # Mass correction (using the monoisotopic m/z)
    dfCorrected, mod = correctMz(df, mzs[0])
    mods[uid] = mod
    # Feature detection
    print("    Identification of MS1 scans where isotopologues (of the metabolite) exist")
    scanIdx = []    # Scan indices where isotopologues coexist
    for i in range(len(mzs)):    # Let's consider only M0 and M1
        mz = mzs[i]
        si = getScanIndex(dfCorrected, mz, tol)   # Scan indices where the peak close to "mz" exists
        if i == 0:
            scanIdx = si
        else:
            scanIdx = np.intersect1d(scanIdx, si)
        if len(scanIdx) == 0:
            break
        # if len(scanIdx) > 0:
        #     scanIdx = np.intersect1d(scanIdx, si)
        #     # overlap = np.intersect1d(scanIdx, si)
        #     # # "overlap" indicates the scan indices where former (e.g., M0, M1 and M2) and current isotopologues (e.g., M3) coexist
        #     # if len(overlap) > 0:
        #     #     scanIdx = overlap
        #     # # Sometimes, M3 may not coexist with former isotologues and the overlap may be empty
        #     # # In this case, the previous "scanIdx" will be kept (M0 ~ M2 coexist, but M3 does not) as is
        # else:
        #     scanIdx = si

    # Feature detection with considering the gap in MS1 scan
    print("    Features are being formed and evaluated whether they form isotopologues")
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
            f = defineFeature(uid, mzs, tol, idx, dfCorrected)
            fArray.append(f)
    print("  Done ...")
print((datetime.now() - t).total_seconds())
print()
