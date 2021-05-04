import os, sys, logging, numpy as np, pandas as pd
from utils import *
from datetime import datetime
from pyteomics import mzxml

# inputFile = r"../../Datasets/hilic_neg/neg_ctrl3.mzXML"
inputFile = r"C:\Users\jcho\OneDrive - St. Jude Children's Research Hospital\UDrive\Research\Projects\7Metabolomics\Datasets\hilic_neg\neg_ctrl3.mzXML"

# Pyteomics reader to pandas DataFrame?
reader = mzxml.read(inputFile)
df = pd.DataFrame(reader)
df.drop(columns=["peaksCount", "polarity", "scanType", "lowMz", "highMz",
                 "basePeakMz", "basePeakIntensity", "totIonCurrent", "id",
                 "collisionEnergy", "precursorMz"], inplace=True)
df = df[df["msLevel"] == 1]
df = df.apply(pd.Series.explode)    # This operation extracts "m/z array" and "intensity array" into columns
print()

# t = datetime.now()
# from featureDetection import detectFeatures
# df = detectFeatures(inputFile, "jumpm_targeted.params")
# print((datetime.now() - t).total_seconds())


# paramFile = "jumpm_targeted.params"
# params = getParams(paramFile)
# print(params)