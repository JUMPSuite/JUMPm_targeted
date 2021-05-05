#!/usr/bin/python

import sqlite3
from utils import *


# def adductDictionary(params):
#     adduct = {}
#     for key, val in params.items():
#         if key.startswith("adduct"):
#             key = re.sub(r'adduct_', '', key)
#             adduct[key] = float(val)
#     return adduct


def calcMS2Similarity(featSpec, libSpec, params):
    # Calculation of MS2 similarity between a feature and a library compound
    # Reference: Clustering millions of tandem mass spectra, J Proteome Res. 2008; 7: 113-22

    # Input arguments
    # featSpec (dictionary): MS2 spectrum of a feature (key = "mz", "intensity")
    # libSpec (dictionary): MS2 spectrum of a library compound (key = "mz", "intensity", "index" (ignorable))
    nPeaks = int(params["num_peaks_ms2_similarity"])    # Default = 30 according to the above reference
    k = min(nPeaks, min(len(featSpec["mz"]), len(libSpec["mz"])))

    # Keep $k strongest peaks in both spectra
    # featDict[mz] = intensity
    # libDict[mz] = intensity
    featDict, libDict = {}, {}
    ind = np.argsort([-i for i in featSpec["intensity"]])
    for i in ind[0:k]:
        featDict[featSpec["mz"][i]] = featSpec["intensity"][i]
    ind = np.argsort([-i for i in libSpec["intensity"]])
    for i in ind[0:k]:
        libDict[libSpec["mz"][i]] = libSpec["intensity"][i]

    # Join two sets of m/z values and make a new set of unique m/z values
    # Duplicate masses are removed as follows
    # - We consider two peaks to have a similar mass if they are within 0.5 Da from each other)
    # - For those two peaks having similar mass, the lower one will be the unique one
    #   (e.g. One peak with m/z = 100 and the other peak with m/z = 100.4 -> they will be merged to m/z = 100)
    mzArray = list(featDict.keys()) + list(libDict.keys())
    mzArray = sorted(mzArray)
    mzDict = {}
    val = 0
    for mz in mzArray:
        if abs(mz - val) <= 0.5:
            mzDict[mz] = val
        else:
            mzDict[mz] = mz
            val = mz

    # Reduction of spectrum to a vector by assigning to each intensity to the unique m/z bins
    # And then, calculate the similarity; normalized dot-product
    s = {}
    for key, val in mzDict.items():
        s[val] = {}
        s[val]["feat"] = 0
        s[val]["lib"] = 0

    for key, val in mzDict.items():
        if key in featDict:
            s[val]["feat"] += np.sqrt(featDict[key])
        if key in libDict:
            s[val]["lib"] += np.sqrt(libDict[key])

    num, den1, den2 = 0, 0, 0    # numerator, denominator1, denominator2
    for mz in s.keys():
        num += s[mz]["feat"] * s[mz]["lib"]
        den1 += s[mz]["feat"] ** 2
        den2 += s[mz]["lib"] ** 2

    if den1 * den2 == 0:
        normDotProduct = 0
    else:
        normDotProduct = num / np.sqrt(den1 * den2)
    return normDotProduct


def queryLibrary(mz, m, z, conn, tol):
    # Standard library search
    # Retrieve library compounds satisfying conditions and calculate MS2-based and RT-based similarity (if exist)

    # Query is made using m/z, neutral mass and charge (if available) information of each feature
    # Neutral mass is the main variable of query
    # m/z is used to prevent searching adduct compounds when 'no adduct' feature is queried (and vice versa)
    # So, a larger tolerance (2 * tol) is applied to m/z-based query not to affect neutral mass-based query
    if z == 0:
        sqlQuery = r"SELECT * FROM library WHERE abs(((?) - mass) / mass * 1e6) < (?) " \
                   r"AND abs(((?) - precursor_mz) / precursor_mz * 1e6) < (?)"
        df = pd.read_sql_query(sqlQuery, conn, params=(m, tol, mz, 2 * tol))
        # # Adduct search
        # dfAdduct = pd.DataFrame()
        # for k, v in adducts.items():
        #     dfAdduct = dfAdduct.append(pd.read_sql_query(sqlQuery, conn, params=(m - v, tol, mz, 2 * tol)), ignore_index=True)
        # df = df.append(dfAdduct, ignore_index=True)
    else:
        sqlQuery = r"SELECT * FROM library WHERE abs(((?) - mass) / mass * 1e6) < (?) " \
                   r"AND abs(((?) - precursor_mz) / precursor_mz * 1e6) < (?) AND charge = (?)"
        df = pd.read_sql_query(sqlQuery, conn, params=(m, tol, mz, 2 * tol, int(z)))
        # # Adduct search
        # dfAdduct = pd.DataFrame()
        # for k, v in adducts.items():
        #     dfAdduct = dfAdduct.append(pd.read_sql_query(sqlQuery, conn, params=(m - v, tol, mz, 2 * tol, int(z))), ignore_index=True)
        # df = df.append(dfAdduct, ignore_index=True)

    return df


def searchLibrary(df, paramFile):
    print("  Calculating MS2-based scores for the features")

    ##################################
    # Load parameters and initialize #
    ##################################
    try:
        params = getParams(paramFile)
    except:
        sys.exit("Parameter file cannot be found or cannot be loaded")
    condition = params["LC_column"].lower()
    if params["mode"] == "1":
        condition = condition + "p"
    elif params["mode"] == "-1":
        condition = condition + "n"
    else:
        sys.exit("'mode' parameter should be either 1 or -1")
    proton = 1.007276466812
    matchMzTol = float(params["library_mass_tolerance"])  # Unit of ppm
    nFeatures = df.shape[0]
    # While full["feature_RT"] has the unit of minute, the library compounds have RTs in the unit of second
    # So, within this function, full["feature_RT"] needs to be converted to the unit of second
    df["RT"] = df["RT"] * 60

    ##########################
    # Perform library search #
    ##########################
    libFile = params["library"][0]
    try:
        conn = sqlite3.connect(libFile)
    except:
        sys.exit("Library file cannot be found or cannot be loaded.")
    simMs2Array = []
    for i in range(nFeatures):
        # Feature information
        fZ = df["z"].iloc[i]
        fSpec = df["MS2"].iloc[i]
        if np.isnan(fZ) or fSpec is None:  # When MS2 spectrum of the feature is not defined, skip it
            continue
        fMz = df["mz"].iloc[i]
        if params["mode"] == "1":  # Positive mode
            fMass = fZ * (fMz - proton)
        elif params["mode"] == "-1":  # Negative mode
            fMass = fZ * (fMz + proton)

        # Retrieve library compounds of which neutral masses are similar to feature mass
        simMs2 = -1
        dfQuery = queryLibrary(fMz, fMass, fZ, conn, matchMzTol)
        if dfQuery.empty:
            print("  There's no library compound corresponding to a feature at m/z = {}".format(fMz))
        elif dfQuery.shape[0] > 1:
            print("  There are more than one compound corresponding to a feature at m/z = {}".format(fMz))
        else:
            # Unique library compound is found for the i-th feature
            uid = dfQuery["id"].iloc[0]
            uid = uid.replace("##Decoy_", "")
            sqlQuery = r"SELECT * FROM {}".format(uid)
            try:
                libSpec = pd.read_sql_query(sqlQuery, conn)
            except:
                continue
            if not libSpec.empty:
                # Calculate the score based on MS2 spectrum
                libSpec = libSpec.to_dict(orient="list")
                simMs2 = calcMS2Similarity(fSpec, libSpec, params)
        simMs2Array.append(simMs2)

    df["MS2_score"] = simMs2Array
    return df
