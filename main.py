import os, numpy as np, pandas as pd
from featureDetection import detectFeatures
from datetime import datetime
from isotopologues import *
# from featureToMS2 import ms2ForFeatures
# from librarySearch import searchLibrary


def formatOutput(res, isoDf):
    intensityDf = pd.DataFrame()
    pctDf = pd.DataFrame()
    for key, df in isoDf.items():
        intensityArray = np.zeros(res.shape[0])
        pctArray = np.zeros(res.shape[0])

        for uid in df["id"].unique():
            rows = np.where(res["idhmdb"] == uid)
            intensityArray[rows] = list(df[df["id"] == uid]["correctedIntensity"])
            pctArray[rows] = list(df[df["id"] == uid]["labelingPct"])

        intensityDf[key + "_intensity"] = intensityArray
        pctDf[key + "_labelingPct"] = pctArray

    res = pd.concat([res, intensityDf], axis=1)
    res = pd.concat([res, pctDf], axis=1)
    return res


if __name__ == "__main__":
    now = datetime.now()
    nowString = now.strftime("%Y/%m/%d %H:%M:%S")
    print("  " + nowString)

    # Input processing
    # Parameters are
    # 1. List of targeted metabolites (names and formulas)
    # 2. Feature detection-related parameters
    inputFile = "isotope_distribution.xlsx"
    paramFile = "jumpm_targeted.params"
    # mzxmlFiles = [r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\6_nolable.mzXML",
    #               r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML",
    #               r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\8_tracer.mzXML",
    #               r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\9_tracer.mzXML"]
    mzxmlFiles = [
        r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML",
        r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\8_tracer.mzXML"]

    # Input dataframe (from Surendhar's program)
    # What is going to be a "key"? metabolite name? HMDB ID?
    infoDf = pd.read_excel(inputFile)
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
    print()

# Feature detection for the given metabolites
# MS2 processing for the features
# MS2-based score calculation for the features representing the given metabolites
# featureDict = {}
# for mzxmlFile in mzxmlFiles:
#     # Feature detection (with given masses)
#     # The method is almost identical to one for Jump -m unlabeled pipeline
#     # Here, the feature detection is only performed for the peaks around given m/z values
#     # (Every m/z value of MS1 peaks is checked against the given m/z list)
#     features = detectFeatures(mzxmlFile, paramFile)
#
#     # # MS2 spectra for the identified features (need to be consolidated within each run)
#     # features = ms2ForFeatures(features, mzxmlFile, paramFile)   # MS2 spectrum for each feature is added
#     #
#     # # MS2-based score calculation
#     # # It is not identification. The metabolites are already given (name, formula, mass, etc.)
#     # # As a confirmation, they are evaluated by MS2-based scores.
#     # # RT-shifts are also obtained for our information.
#     # features = searchLibrary(features, paramFile)
#     # print()
#     #
#     # featureDict[mzxmlFile] = features
# print()






# Quantification




# Calculation of labeling percentages





