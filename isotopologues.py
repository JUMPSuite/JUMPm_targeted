import numpy as np, pandas as pd
from pyteomics import mzxml


def findIsotopologuePeak(features, spec, givenMz, tol1, tol2):
    lL = givenMz - givenMz * tol2 / 1e6
    uL = givenMz + givenMz * tol2 / 1e6

    # 1. Look for any peak(s) between lL and uL in the given MS1 spectrum (where the strongest M0 peak exists)
    idx = (spec["m/z array"] >= lL) & (spec["m/z array"] <= uL)
    # Case 1. When there's an isotoplogue peak in the MS1 spectrum
    rtShift = 0
    if sum(idx) > 0:
        ms1 = int(spec["num"])
        maxIdx = np.argmax(spec["intensity array"][idx])
        mz = spec["m/z array"][idx][maxIdx]
        intensity = spec["intensity array"][idx][maxIdx]
    # Case 2. When there's NO isotopologue peak in the MS1 spectrum
    else:
        mz, intensity, ms1, rt = findFeature(features, givenMz, tol1)
        if mz > 0:
            rtShift = rt - spec["retentionTime"]
            mz, intensity = 0, 0  # If the isotopologue is not found, no need to return mz and intensity

    return mz, intensity, ms1, rt, rtShift


def findFeature(features, givenMz, tol):
    # For the given m/z, look for the feature corresponding to the m/z
    lL = givenMz - givenMz * tol / 1e6
    uL = givenMz + givenMz * tol / 1e6
    rows = (features["mz"] >= lL) & (features["mz"] <= uL)

    # Once the feature is found, return the representative peak information (i.e., the strongest peak among peaks of the feature)
    if sum(rows) > 0:
        idx = np.argmax(features[rows]["intensity"])
        mz = features[rows].iloc[idx]["mz"]  # Observed m/z of the M0 peak
        intensity = features[rows].iloc[idx]["intensity"]  # Observed intensity of the M0 peak
        ms1 = int(features[rows].iloc[idx]["MS1"])
        rt = features[rows].iloc[idx]["RT"]
        return mz, intensity, ms1, rt
    else:
        return 0, 0, 0, 0

def findIsotopologues(features, mzxmlFile, df, tol1, tol2):
    # Input arguments
    # features = a pandas dataframe of all features
    # mzxmlFile = the mzXML file
    # df = a pandas dataframe of isotopologue information
    # tol1 = m/z tolerance (ppm) for finding the feature corresponding to M0
    # tol2 = m/z tolerance (ppm) for finding the isotopologue peak (from M(i) to M(i+1)) in a MS1 scan

    # Summarize the information of metabolites
    dictM0 = {}
    uids = np.array(df[df["isotopologues"] == "M0"]["idhmdb"])
    theoMzs = np.array(df[df["isotopologues"] == "M0"]["isotope_m/z"])
    for i in range(len(uids)):
        dictM0[uids[i]] = float(theoMzs[i].split(";")[0])

    # Looking for the features corresponding to M0 of metabolites
    reader = mzxml.MzXML(mzxmlFile)
    delC = 1.003355
    res = {"id": [], "mz": [], "intensity": [], "rt": [], "rtShift": []}
    for uid, theoMz in dictM0.items():
        res["id"].append(uid)
        mzArray, intensityArray, ms1Array, rtArray, rtShiftArray = [], [], [], [], []
        ########################
        # Identification of M0 #
        ########################
        # Look for the feature and its representative peak corresponding to M0 of the given metabolite
        mz, intensity, ms1, rt = findFeature(features, theoMz, tol1)


        # TO-DO: what happens if M0 is not found, i.e., mz == None?
        if mz == 0:
            print("M0 for {} is not found".format(uid))
            continue


        # Put M0 information to mzArray and intensityArray
        mzArray.append(mz)
        intensityArray.append(intensity)
        ms1Array.append(ms1)
        rtArray.append(rt)
        rtShiftArray.append(0)

        #####################################
        # Identification of M1, M2, ..., Mn #
        #####################################
        # At the "repMs1" scan, look for isotopologue peaks, i.e., M1, M2, ..., Mn
        spec = reader[str(int(ms1))]
        nIsotopologues = sum(df["idhmdb"] == uid)
        for i in range(1, nIsotopologues):
            mz += delC    # Suppose that the tracer is 13C
            obsMz, obsIntensity, obsMs1, obsRt, obsRtShift = findIsotopologuePeak(features, spec, mz, tol1, tol2)
            if obsMz > 0:
                mz = obsMz  # When an isotopologue is found, the next one will be searched from the current one

            mzArray.append(obsMz)
            intensityArray.append(obsIntensity)
            ms1Array.append(obsMs1)
            rtArray.append(obsRt)
            rtShiftArray.append(obsRtShift)

        res["mz"].append(mzArray)
        res["intensity"].append(intensityArray)
        res["rt"].append(rtArray)
        res["rtShift"].append(rtShiftArray)

    res = pd.DataFrame.from_dict(res)
    return res


def correctionMatrix(df):
    res = {}
    uids = df["idhmdb"].unique()
    for uid in uids:
        subDf = df[df["idhmdb"] == uid]
        n = subDf.shape[0]
        cm = np.zeros((n, n))
        for j in range(n):
            vec = subDf.iloc[j]["isotope_intensity"].split(";")
            for i in range(len(vec)):
                if i + j < n:
                    cm[i + j, j] = float(vec[i]) / 100
                else:
                    continue
        res[uid] = cm
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


if __name__ == "__main__":
    # import pickle
    # [res, isoDf] = pickle.load(open("tmp.pickle", "rb"))
    # intensityDf = pd.DataFrame()
    # pctDf = pd.DataFrame()
    # for key, df in isoDf.items():
    #     intensityArray = np.zeros(res.shape[0])
    #     pctArray = np.zeros(res.shape[0])
    #     for uid in df["id"].unique():
    #         rows = np.where(res["idhmdb"] == uid)
    #         intensityArray[rows] = list(df[df["id"] == uid]["correctedIntensity"])
    #         pctArray[rows] = list(df[df["id"] == uid]["labelingPct"])
    #
    #     intensityDf[key + "_intensity"] = intensityArray
    #     pctDf[key + "_labelingPct"] = pctArray
    #
    # res = pd.concat([res, intensityDf], axis=1)
    # res = pd.concat([res, pctDf], axis=1)
    # print()

    inputFile = "isotope_distribution.xlsx"
    mzxmlFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML"
    features = pd.read_csv("features_7_tracer.txt", sep="\t")
    df = pd.read_excel(inputFile)
    tol1, tol2 = 15, 5
    res = findIsotopologues(features, mzxmlFile, df, tol1, tol2)

