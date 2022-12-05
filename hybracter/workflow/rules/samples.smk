"""
Function for parsing the 'Assemblies' config and identifying samples and read files
"""

from itertools import chain

# def samplesFromDirectory(dir):
#     """Parse samples from a directory"""
#     outDict = {}
#     # https://stackoverflow.com/questions/11860476/how-to-unnest-a-nested-list
#     samples= glob_wildcards(os.path.join(dir,'{sample}.fastq.gz'))
#     samples2 = chain(*samples)
#     for sample in samples2:
#         outDict[sample] = {}
#         fastq = os.path.join(dir,f'{sample}.fastq.gz')
#         if os.path.isfile(fastq):
#             outDict[sample]['fastq'] = fastq
#         else:
#             sys.stderr.write("\n"
#                              "    FATAL: Error globbing files."
#                              f"    {fastq} \n"
#                              "    does not exist. Ensure consistent formatting and file extensions."
#                              "\n")
#             sys.exit(1)
#     return outDict

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


def parseSamples(csvfile, LR_ONLY):
    # for reading from directory
    #if os.path.isdir(readFileDir):
    #   sampleDict = samplesFromDirectory(readFileDir)
    if os.path.isfile(csvfile) and LR_ONLY == True:
        sampleDict = samplesFromCsvLong(csvfile)
    elif os.path.isfile(csvfile) and LR_ONLY == False:
        sampleDict = samplesFromCsvShort(csvfile)
    else:
        sys.stderr.write("\n"
                         f"    FATAL: {csvfile} is neither a file nor directory.\n"
                         "\n")
        sys.exit(1)
    if len(sampleDict.keys()) == 0:
        sys.stderr.write("\n"
                         "    FATAL: We could not detect any samples at all.\n"
                         "\n")
        sys.exit(1)
    return sampleDict