import re, sys, os, pickle, numpy as np, pandas as pd


def getParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line)  # Remove all whitespaces

            # Exception for "feature_files" parameter
            if "feature_files" in parameters and line.endswith(".feature"):
                parameters["feature_files"].append(line)
            elif "library" in parameters and line.endswith(".db"):
                parameters["library"].append(line)
            else:
                key = line.split('=')[0]
                val = line.split('=')[1]
                if key == "feature_files" or key == "library":
                    parameters[key] = [val]
                else:
                    parameters[key] = val
    return parameters


def getMs1(reader, params):
    try:
        firstScan = int(params["first_scan_extraction"])
    except:
        firstScan = 1
    try:
        lastScan = int(params["last_scan_extraction"])
    except:
        lastScan = 1000000
    ms1 = []
    ind = 0
    with reader:
        for spec in reader:
            msLevel = spec["msLevel"]  # int type
            scanNum = spec["num"]  # str type
            if msLevel == 1 and firstScan <= int(scanNum) <= lastScan:
                spec["scanIndex"] = ind
                ms1.append(spec)
                ind += 1
            elif int(scanNum) > lastScan:
                break
    return ms1


def readFeatures(featureFile):
    with open(featureFile, 'r') as file:
        features = file.readlines()
    return features


class progressBar:
    def __init__(self, total):
        self.total = total
        self.barLength = 20
        self.count = 0
        self.progress = 0
        self.block = 0
        self.status = ""

    def increment(self, nIncrement=None):
        if nIncrement == None:
            self.count += 1
        else:
            self.count = nIncrement
        self.progress = self.count / self.total
        self.block = int(round(self.barLength * self.progress))
        if self.progress == 1:
            self.status = "Done...\r\n"
        else:
            self.status = ""
        #         self.status = str(self.count) + "/" + str(self.total)
        text = "\r  Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block),
                                                     int(self.progress * 100), self.status)
        sys.stdout.write(text)
        sys.stdout.flush()
