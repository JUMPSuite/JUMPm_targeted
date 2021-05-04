from featureDetection import detectFeatures

# Input processing
# Parameters are
# 1. List of targeted metabolites (names and formulas)
# 2. Feature detection-related parameters
paramFile = "jumpm_targeted.params"
mzxmlFiles = [r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_ctrl3.mzXML",
              r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_HK1_4.mzXML",
              r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_QC3.mzXML"]

# Input dataframe (from Surendhar's program)
# What is going to be a "key"? metabolite name? HMDB ID?


# Feature detection (with given masses)
# The method is almost identical to one for Jump -m unlabeled pipeline
# Here, the feature detection is only performed for the peaks around given m/z values
# (Every m/z value of MS1 peaks is checked against the given m/z list)
features = {}
for mzxmlFile in mzxmlFiles:
    features = detectFeatures(mzxmlFile, paramFile)
    print()
print()


# MS2-based score calculation
# It is not identification. The metabolites are already given (name, formula, mass, etc.)
# As a confirmation, they are evaluated by MS2-based scores.
# RT-shifts are also obtained for our information.




# Quantification




# Calculation of labeling percentages





