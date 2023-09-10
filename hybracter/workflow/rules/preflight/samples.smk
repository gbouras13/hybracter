"""
Function for parsing the config and identifying samples and read files
"""

from itertools import chain

"""
long
"""

def samplesFromCsvLong(csvFile):
    """
    Read samples and files from a CSV Long Read Only
    3 cols
    1 = sample 
    2 = Long read
    3 = MinChromLength
    """
    outDict = {}
    with open(csvFile,'r') as csv:
        for line in csv:
            l = line.strip().split(',')
            if len(l) == 3:
                outDict[l[0]] = {}
                if os.path.isfile(l[1]) and l[2].isnumeric():
                    outDict[l[0]]['LR'] = l[1]
                    outDict[l[0]]['MinChromLength'] = l[2]
                else:
                    sys.stderr.write("\n"
                                     f"    FATAL: Error parsing {csvFile}. {l[1]} \n"
                                     f"    does not exist or  {l[2]} is not an integer. \n"
                                     "    Check formatting, and that \n" 
                                     "    file names and file paths are correct.\n"
                                     "\n")
                    sys.exit(1)
    return outDict

"""
short
"""

def samplesFromCsvShort(csvFile):
    """
    Read samples and files from a CSV Hybrid
    3 cols
    1 = sample 
    2 = Long read
    3 = MinChromLength
    4 = R1 Short
    5 = R2 Short
    """
    outDict = {}
    with open(csvFile,'r') as csv:
        for line in csv:
            l = line.strip().split(',')
            if len(l) == 5:
                outDict[l[0]] = {}
                if os.path.isfile(l[1])  and l[2].isnumeric() and os.path.isfile(l[3]) and os.path.isfile(l[4]):
                    outDict[l[0]]['LR'] = l[1]
                    outDict[l[0]]['MinChromLength'] = l[2]
                    outDict[l[0]]['R1'] = l[3]
                    outDict[l[0]]['R2'] = l[4]
                else:
                    sys.stderr.write("\n"
                                     f"    FATAL: Error parsing {csvFile}. One of \n"
                                     f"    {l[1]} or \n"
                                     f"    {l[3]} or \n"
                                     f"    {l[4]} \n"
                                     f"    does not exist or  {l[2]} is not an integer. \n"
                                     "    Check formatting, and that \n" 
                                     "    file names and file paths are correct.\n"
                                     "\n")
                    sys.exit(1)
    return outDict


def parseSamples(csvfile, long_flag):
    if os.path.isfile(csvfile) and long_flag is True:
        sampleDict = samplesFromCsvLong(csvfile)
    elif os.path.isfile(csvfile) and long_flag is False:
        sampleDict = samplesFromCsvShort(csvfile)
    else:
        sys.stderr.write("\n"
                        f"    FATAL: something is wrong. Likely {csvfile} is neither a file nor directory.\n"
                        "\n")
        sys.exit(1)

    return sampleDict