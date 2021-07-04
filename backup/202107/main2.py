import os, numpy as np, pandas as pd
from pyteomics import mzxml
from datetime import datetime
from utils import *

def findPeak(spec, givenMz, tol):
    lL = givenMz - givenMz * tol / 1e6
    uL = givenMz + givenMz * tol / 1e6
    idx = (spec["m/z array"] >= lL) & (spec["m/z array"] <= uL)
    if sum(idx) > 0:    # A peak close to "givenMz" is found in "spec"
        maxIdx = np.argmax(spec["intensity array"][idx])
        mz = spec["m/z array"][idx][maxIdx]
        intensity = spec["intensity array"][idx][maxIdx]
    else:    # No peak is found around "givenMz"
        mz, intensity = 0, 0

    return mz, intensity


def findIsotopologue(mzxmlFile, infoDf, isRef, params):
    # Summarize the information of metabolites
    dictM0 = {}
    uids = np.array(infoDf[infoDf["isotopologues"] == "M0"]["idhmdb"])
    mzs = np.array(infoDf[infoDf["isotopologues"] == "M0"]["feature_m/z"])
    rts = np.array(infoDf[infoDf["isotopologues"] == "M0"]["feature_RT"])
    for i in range(len(uids)):
        dictM0[uids[i]] = {"mz": mzs[i], "rt": rts[i]}

    # Initialization
    delC = 1.003355
    tol = 5
    reader = mzxml.MzXML(mzxmlFile)
    df = pd.DataFrame(reader)
    res = {"id":[], "ms1":[], "rt":[], "mz":[], "intensity":[]}

    for uid in dictM0.keys():
        res["id"].append(uid)
        mzArray, intensityArray, ms1Array, rtArray = [], [], [], []

        # Extract m/z and RT informatio of the "uid"
        mz = dictM0[uid]["mz"]
        rt = dictM0[uid]["rt"]

        if isRef == 1:  # When the current run is a reference run
            idx = np.argmin(abs(df["retentionTime"] - rt))
            scanNum = df.iloc[idx]["num"]
        else:  # When the current run is not a reference run
            # Consider MS1 spectra within +/-2.5 min of the reference RT
            subDf = df[(df["msLevel"] == 1) & (df["retentionTime"] > rt - 2.5) & (df["retentionTime"] < rt + 2.5)]
            subDf = subDf.set_index(["retentionTime"]).apply(pd.Series.explode).reset_index()
            # Look for the MS1 peak whose m/z is withint the tolerance and intensity is the strongest among candidates
            subDf = subDf[(subDf["m/z array"] >= mz - mz * tol / 1e6) & (subDf["m/z array"] <= mz + mz * tol / 1e6)]
            subDf = subDf.sort_values(by="intensity array", ascending=False)
            scanNum = subDf.iloc[0]["num"]

        spec = reader[scanNum]
        rt = spec["retentionTime"]
        mz, intensity = findPeak(spec, mz, tol)
        mzArray.append(mz)
        intensityArray.append(intensity)
        ms1Array.append(int(scanNum))
        rtArray.append(rt)

        #####################################
        # Identification of M1, M2, ..., Mn #
        #####################################
        nIsotopologues = sum(infoDf["idhmdb"] == uid)
        for i in range(1, nIsotopologues):
            mz += delC  # Suppose that the tracer is 13C
            obsMz, obsIntensity = findPeak(spec, mz, tol)
            if obsMz > 0:
                mz = obsMz  # When an isotopologue is found, the next one will be searched from the current one
            else:
                scanNum, rt = 0, 0
            mzArray.append(obsMz)
            intensityArray.append(obsIntensity)
            ms1Array.append(int(scanNum))
            rtArray.append(rt)

        res["mz"].append(mzArray)
        res["intensity"].append(intensityArray)
        res["ms1"].append(ms1Array)
        res["rt"].append(rtArray)

    res = pd.DataFrame.from_dict(res)
    return res


def quantifyIsotopologues(infoDf, isoDf):
    # Input arguments
    # infoDf = a pandas dataframe containing the information of isotopologues, e.g., theoretical isotopic peak m/z and intensitry
    # isoDf = a pandas dataframe containing the identified isotopologues of given metabolites

    # Quantification of isotopologues
    intensityArray, pctArray = [], []
    cm = correctionMatrix(infoDf)   # Correction matrix derived from the theoretical information of isotopologues
    for i in range(isoDf.shape[0]):
        uid = isoDf.loc[i]["id"]
        intensity = isoDf.loc[i]["intensity"]
        correctedIntensity = np.dot(np.linalg.inv(cm[uid]), intensity)
        correctedIntensity[correctedIntensity < 0] = 0
        intensityArray.append(correctedIntensity)
        pctArray.append(correctedIntensity / sum(correctedIntensity) * 100)
    isoDf["correctedIntensity"] = intensityArray
    isoDf["labelingPct"] = pctArray

    return isoDf


def correctionMatrix(df):
    res = {}
    uids = df["idhmdb"].unique()
    for uid in uids:
        subDf = df[df["idhmdb"] == uid]
        n = subDf.shape[0]
        cm = np.zeros((n, n))
        for i in range(n):
            # Assume that "intensity" is already sorted and organized from M0 to Mn
            intensity = subDf.iloc[i]["isotope_intensity"].split(";")
            cm[i, :] = [float(val) / 100 for val in intensity]
        res[uid] = cm

    return res


def formatOutput(res, isoDf):
    intensityDf = pd.DataFrame()
    pctDf = pd.DataFrame()
    n = res.shape[0]
    mzArray = np.zeros(n)
    ms1Array = np.zeros(n)
    rtArray = np.zeros(n)

    for key, df in isoDf.items():
        # Explode columns of the dataframe
        df = df.set_index(["id"]).apply(pd.Series.explode).reset_index()

        # Add the corrected intensity and labeling percentage to the dataframes
        intensityDf[key.split(".")[0] + "_intensity"] = df["correctedIntensity"]
        pctDf[key.split(".")[0] + "_labelingPct"] = df["labelingPct"]

        # Add other information to the corresponding arrays
        if np.all(ms1Array == 0):
            mzArray = np.array(df["mz"])
            ms1Array = np.array(df["ms1"])
            rtArray = np.array(df["rt"])
        else:
            for i in range(n):
                mzArray[i] = str(mzArray[i]) + ";" + str(df.loc[i]["mz"])
                ms1Array[i] = str(ms1Array[i]) + ";" + str(df.loc[i]["ms1"])
                rtArray[i] = str(rtArray[i]) + ";" + str(df.loc[i]["rt"])

    res["observed_m/z"] = mzArray
    res["MS1scan"] = ms1Array
    res["RT"] = rtArray
    res = pd.concat([res, intensityDf], axis=1)
    res = pd.concat([res, pctDf], axis=1)

    return res


if __name__ == "__main__":

    startTime = datetime.now()
    startTimeString = startTime.strftime("%Y/%m/%d %H:%M:%S")
    print("  " + startTimeString)

    refDf = pd.read_csv("6_nolable_jumpm.csv")
    infoDf = pd.read_csv("isotope_distribution.txt", sep="\t")
    res = infoDf.copy()
    # Temporarily
    infoDf = infoDf.merge(refDf, left_on="name", right_on="name")
    mzxmlFiles = [
        r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\6_nolable.mzXML",
        r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML",
        r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\8_tracer.mzXML",
        r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\9_tracer.mzXML"]
    paramFile = "jumpm_targeted.params"
    params = getParams(paramFile)
    isoDf = {}

    for mzxmlFile in mzxmlFiles:
        if os.path.basename(mzxmlFile) == "6_nolable.mzXML":
            isRef = 1
        else:
            isRef = 0
        df = findIsotopologue(mzxmlFile, infoDf, isRef, params)
        df = quantifyIsotopologues(infoDf, df)
        isoDf[os.path.basename(mzxmlFile)] = df

    # Format the output dataframe
    res = formatOutput(res, isoDf)
    res.to_csv("tracer_result2.txt", sep="\t", index=False)

    elapsed = (datetime.now() - startTime).total_seconds()
    print("  Took {} seconds".format(int(elapsed)))
