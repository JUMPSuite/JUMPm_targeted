import os, numpy as np, pandas as pd
from pyteomics import mzxml
from datetime import datetime
from utils import *
# from isotopeCalculation import *
from isotopeCalcEmass import *


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
    ms1Scans = []
    ms1RTs = []
    for spec in reader:
        if spec["msLevel"] == 1:
            ms1Scans.append(spec["num"])
            ms1RTs.append(spec["retentionTime"])

    res = {"id":[], "ms1":[], "rt":[], "mz":[], "intensity":[], "pct": []}
    for uid in dictM0.keys():
        res["id"].append(uid)
        mzArray, intensityArray, ms1Array, rtArray = [], [], [], []

        # Extract m/z and RT information of the "uid"
        mz = dictM0[uid]["mz"]
        rt = dictM0[uid]["rt"]
        intensity = 0

        # Look for the monoisotopic peak of "uid" (i.e., M0)
        scanNum = 0
        if isRef == 1:  # When the current run is a reference run
            idx = np.argmin(abs(ms1RTs - rt))
            scanNum = ms1Scans[idx]
        else:  # When the current run is not a reference run
            idxes = [i for i, v in enumerate(ms1RTs) if (rt - 2.5) < v < (rt + 2.5)]
            if len(idxes) > 0:
                maxIntensity = 0
                for idx in idxes:
                    mzs = reader[ms1Scans[idx]]["m/z array"]
                    ints = reader[ms1Scans[idx]]["intensity array"]
                    j = np.argmin(abs(mzs - mz))
                    if mzs[j] >= (mz - mz * tol / 1e6) and mzs[j] <= (mz + mz * tol / 1e6) and ints[j] > maxIntensity:
                        maxIntensity = ints[j]
                        scanNum = ms1Scans[idx]

        if scanNum != 0:    # When there is M0 of "uid"
            # mz, rt and intensity are replaced with the observed ones
            spec = reader[scanNum]
            rt = spec["retentionTime"]
            mz, intensity = findPeak(spec, mz, tol)
        else:   # When there's no peak corresponding to M0 of "uid"
            pass    # mz and rt: from the feature of the reference run, intensity = 0

        mzArray.append(mz)
        intensityArray.append(intensity)
        ms1Array.append(int(scanNum))
        rtArray.append(rt)

        # Identification of M1, M2, ..., Mn
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
        res["pct"].append(intensityArray / sum(intensityArray) * 100)

    res = pd.DataFrame.from_dict(res)
    return res


def correctNaturalAbundance(df):
    # Input arguments
    # inputDf = a pandas dataframe containing the information of isotopologues and their quantity (uncorrected)

    # Quantification of isotopologues
    intensityArray, pctArray = [], []
    cm = correctionMatrix(df)   # Correction matrix derived from the theoretical information of isotopologues
    cols = [s for s in df.columns if s.endswith("intensity") and s != "isotope_intensity"]
    uids = df["idhmdb"].unique()
    for uid in uids:
        idx = df["idhmdb"] == uid
        for col in cols:
            intensity = df.loc[idx][col]
            correctedIntensity = np.dot(np.linalg.inv(cm[uid]), intensity)
            correctedIntensity[correctedIntensity < 0] = 0
            df.loc[idx, col] = correctedIntensity
            correctedPct = correctedIntensity / sum(correctedIntensity) * 100
            df.loc[idx, col.replace("intensity", "labelingPct")] = correctedPct

    return df


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

        # Note that the correction matrix (i.e., natural abundance matrix) should be arranged so that cm[i, j] represents
        # the fraction of the distribution of the j-th labeled species (i.e., isotopologues) corresponding to the i-th measured value (i.e., m/z)
        # Each row = m/z
        # Each column = isotopologue
        # The fraction should be scaled to be summed to 1 for each column
        res[uid] = cm.T

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
        intensityDf[key.split(".")[0] + "_intensity"] = df["intensity"]
        pctDf[key.split(".")[0] + "_labelingPct"] = df["pct"]

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


def getIsotopicDistributions(df, params):
    # Input arguments
    # df = a pandas dataframe containing the JUMPm result of the reference run
    #      "formula" is used for the input argument of calculating isotopic peaks of isotopologues
    #      "name" is used as a key for merging dataframes later
    # params = a dictionary of parameters

    res = pd.DataFrame()
    for i in range(df.shape[0]):
        formula = df.formula[i]
        if df.feature_ion[i][-1] == "+":
            charge = df.feature_z[i]
        elif df.feature_ion[i][-1] == "-":
            charge = -df.feature_z[i]
        res_i = getIsotopologues(formula, charge, params)
        res_i["name"] = df.name[i]
        res = res.append(res_i, ignore_index=True)

    return res


if __name__ == "__main__":
    startTime = datetime.now()
    startTimeString = startTime.strftime("%Y/%m/%d %H:%M:%S")
    print("  " + startTimeString)
    print("  JUMPm for targeted metabolites\n")

    # Input arguments
    # 1. A parameter file containing the following,
    #    - A table (.csv) containing the feature information of target metabolites in a reference run (sample)
    #      It includes the information of target metabolites, operation mode (negative or positive), charge and so on
    #    - Parameters for calculating theoretical natural abundances of target metabolitesList of target metabolites
    #    - Experimental conditions including a tracer, MS mode (pos or neg), etc.
    # 2.
    # 3. List of mzXML files

    paramFile = sys.argv[1]
    # paramFile = "jumpm_targeted.params"
    params = getParams(paramFile)

    # Mode 1, identification and quantification of isotopologues (of given target metabolites)
    if params["mode"] == "1":
        print("  Isotopologues of given target metabolites are identified and quantified")
        mzxmlFiles = sys.argv[2:]
        # mzxmlFiles = [
        #     r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\6_nolable.mzXML",
        #     r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML",
        #     r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\8_tracer.mzXML",
        #     r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\9_tracer.mzXML"]

        if len(mzxmlFiles) == 0:
            sys.exit("  You should specify mzXML files\n  e.g., jump -mpython -target jumpm_targeted.params test1.mzXML test2.mzXML ...")

        # Calculation of theoretical isotopic distributions (Surendhar's script)
        refInfoFile = params["ref_feature_information"]  # JUMPm result of the reference run
        refDf = pd.read_csv(refInfoFile)
        infoDf = getIsotopicDistributions(refDf, params)
        # infoDf = getIsotopicDistributions(paramFile, refInfoFile)
        infoDf = refDf.merge(infoDf, left_on="name", right_on="name")
        res = infoDf.copy()
        isoDf = {}
        for mzxmlFile in mzxmlFiles:
            print("  Working on {}".format(os.path.basename(mzxmlFile)))
            if os.path.basename(mzxmlFile) == params["ref_run"]:
                isRef = 1
            else:
                isRef = 0
            df = findIsotopologue(mzxmlFile, infoDf, isRef, params)
            isoDf[os.path.basename(mzxmlFile)] = df
        res = res[["idhmdb", "formula", "name", "feature_ion", "feature_z", "isotopologues", "isotope_m/z", "isotope_intensity"]]
        res = res.rename(columns={"feature_ion": "ion", "feature_z": "charge"})

        # Format the output dataframe
        res = formatOutput(res, isoDf)
        res.to_csv("tracer_result_emass.txt", sep="\t", index=False)

    # Mode 2, correction of natural abundances of isotopic peaks (of quantified isotopologues)
    elif params["mode"] == "2":
        print("  Quantity data of isotopologues is corrected by natural abundances")
        try:
            inputFile = params["quan_result"]
        except KeyError:
            sys.exit("  'quan_result' parameter should be correctly specified")
        try:
            df = pd.read_csv(inputFile, sep="\t")
        except FileNotFoundError:
            sys.exit("  'Please check 'quan_result' parameter whether the file path is correctly specified")
        res = correctNaturalAbundance(df)
        res.to_csv("tracer_corrected_result.txt", sep="\t", index=False)
    else:
        sys.exit("The parameter 'mode' should be properly set (either 1 or 2)")

    print()
    endTime = datetime.now()
    endTimeString = endTime.strftime("%Y/%m/%d %H:%M:%S")
    print("  " + endTimeString)
    elapsed = (endTime - startTime).total_seconds()
    print("  Finished in {} seconds".format(int(elapsed)))
