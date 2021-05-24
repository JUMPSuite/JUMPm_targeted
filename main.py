import pandas as pd
from featureDetection import detectFeatures
from featureToMS2 import ms2ForFeatures
from librarySearch import searchLibrary

# Input processing
# Parameters are
# 1. List of targeted metabolites (names and formulas)
# 2. Feature detection-related parameters
inputFile = "isotope_distribution.xlsx"
paramFile = "jumpm_targeted.params"
# mzxmlFiles = [r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_ctrl3.mzXML",
#               r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_HK1_4.mzXML",
#               r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_QC3.mzXML"]
mzxmlFiles = [r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\13Ctracer_rawdata\7_tracer.mzXML"]


# Input dataframe (from Surendhar's program)
# What is going to be a "key"? metabolite name? HMDB ID?
df = pd.read_excel(inputFile)

# Feature detection for the given metabolites
# MS2 processing for the features
# MS2-based score calculation for the features representing the given metabolites
featureDict = {}
for mzxmlFile in mzxmlFiles:
    # Feature detection (with given masses)
    # The method is almost identical to one for Jump -m unlabeled pipeline
    # Here, the feature detection is only performed for the peaks around given m/z values
    # (Every m/z value of MS1 peaks is checked against the given m/z list)
    features = detectFeatures(df, mzxmlFile, paramFile)

    # MS2 spectra for the identified features (need to be consolidated within each run)
    features = ms2ForFeatures(features, mzxmlFile, paramFile)   # MS2 spectrum for each feature is added

    # MS2-based score calculation
    # It is not identification. The metabolites are already given (name, formula, mass, etc.)
    # As a confirmation, they are evaluated by MS2-based scores.
    # RT-shifts are also obtained for our information.
    features = searchLibrary(features, paramFile)
    print()

    featureDict[mzxmlFile] = features
print()






# Quantification




# Calculation of labeling percentages





