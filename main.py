import os, numpy as np, pandas as pd
from featureDetection import detectFeatures
from datetime import datetime
from isotopologues import *
# from featureToMS2 import ms2ForFeatures
# from librarySearch import searchLibrary


def formatOutput(res, isoDf):
    intensityDf = pd.DataFrame()
    pctDf = pd.DataFrame()
    n = res.shape[0]
    mzArray = np.zeros(n)
    ms1Array = np.zeros(n)
    rtArray = np.zeros(n)
    rtShiftArray = np.zeros(n)
    scoreArray = np.zeros(n)

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
            rtShiftArray = np.array(df["rtShift"])
            scoreArray = np.array(df["ms2Score"])
        else:
            for i in range(n):
                mzArray[i] = str(mzArray[i]) + ";" + str(df.loc[i]["mz"])
                ms1Array[i] = str(ms1Array[i]) + ";" + str(df.loc[i]["ms1"])
                rtArray[i] = str(rtArray[i]) + ";" + str(df.loc[i]["rt"])
                rtShiftArray[i] = str(rtShiftArray[i]) + ";" + str(df.loc[i]["rtShift"])
                scoreArray[i] = str(scoreArray[i]) + ";" + str(df.loc[i]["ms2Score"])

    res["observed_m/z"] = mzArray
    res["MS1scan"] = ms1Array
    res["RT"] = rtArray
    res["RTshift"] = rtShiftArray
    res["MS2score"] = scoreArray
    res = pd.concat([res, intensityDf], axis=1)
    res = pd.concat([res, pctDf], axis=1)

    return res



if __name__ == "__main__":
    now = datetime.now()
    nowString = now.strftime("%Y/%m/%d %H:%M:%S")
    print("  " + nowString)

    # Input processing
    infoDf = pd.read_csv("isotope_distribution.txt", sep="\t")  # Isotopologue m/z and intensity information
    refDf = pd.read_csv("6_nolable_jumpm.csv")  # Jump -m result of the reference (unlabeled) run
    paramFile = "jumpm_targeted.params"
    mzxmlFiles = [r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\6_nolable.mzXML",
                  r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML",
                  r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\8_tracer.mzXML",
                  r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\9_tracer.mzXML"]

    # Feature-based approach (last seen on 6/20/2021)

    # Input processing
    # Parameters are
    # 1. List of targeted metabolites (names and formulas)
    # 2. Feature detection-related parameters
    inputFile = "isotope_distribution.txt"
    paramFile = "jumpm_targeted.params"
    mzxmlFiles = [r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\6_nolable.mzXML",
                  r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML",
                  r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\8_tracer.mzXML",
                  r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\9_tracer.mzXML"]

    # Input dataframe (from Surendhar's program)
    # What is going to be a "key"? metabolite name? HMDB ID?
    infoDf = pd.read_csv(inputFile, sep="\t")
    isoDf = {}
    res = infoDf.copy()
    for mzxmlFile in mzxmlFiles:
        # 1. Feature detection
        features = detectFeatures(mzxmlFile, paramFile)
        # 2. Identification of isotopologues
        df = findIsotopologues(features, mzxmlFile, infoDf, 15, 5)
        # 3. Quantification of isotopologues
        df = quantifyIsotopologues(infoDf, df)
        isoDf[os.path.basename(mzxmlFile)] = df

    # Format the output dataframe
    res = formatOutput(res, isoDf)
    res.to_csv("tracer_result.txt", sep="\t", index=False)

