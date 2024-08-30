"""
Function for parsing the config and identifying samples and read files
"""

from itertools import chain
import sys

"""
long
"""


def samplesFromCsvLong(csvFile, subsample_depth, datadir):
    """
    Read samples and files from a CSV Long Read Only
    3 cols
    1 = sample
    2 = Long read
    3 = MinChromLength
    """
    outDict = {}
    with open(csvFile, "r") as csv:
        for line in csv:
            l = line.strip().split(",")
            if len(l) == 3:
                outDict[l[0]] = {}
                # where datafir isn't specified
                if datadir is None:
                    long_fastq = l[1]
                else:
                    long_fastq = os.path.join(datadir, l[1])
                if os.path.isfile(long_fastq) and l[2].isnumeric():
                    outDict[l[0]]["LR"] = long_fastq
                    outDict[l[0]]["MinChromLength"] = l[2]
                    outDict[l[0]]["TargetBases"] = int(l[2]) * subsample_depth
                else:
                    sys.stderr.write(
                        "\n"
                        f"    FATAL: Error parsing {csvFile}. {long_fastq} \n"
                        f"    does not exist or  {l[2]} is not an integer. \n"
                        "    Check formatting, and that \n"
                        "    file names and file paths are correct.\n"
                        "\n"
                    )
                    sys.exit(1)
            else:
                sys.stderr.write(
                    "\n"
                    f"    FATAL: Error parsing {csvFile}. Line {l} \n"
                    f"    does not have 3 columns. \n"
                    f"    Please check the formatting of {csvFile}. \n"
                )
                sys.exit(1)
    return outDict


"""
short
"""


def samplesFromCsvShort(csvFile, subsample_depth, datadir):
    """
    Read samples and files from a CSV Hybrid
    5 cols
    1 = sample
    2 = Long read
    3 = MinChromLength
    4 = R1 Short
    5 = R2 Short
    """
    outDict = {}
    with open(csvFile, "r") as csv:
        for line in csv:
            l = line.strip().split(",")
            if len(l) == 5:
                outDict[l[0]] = {}
                # where datafir isn't specified
                if datadir is None:
                    long_fastq = l[1]
                    r1_fastq = l[3]
                    r2_fastq = l[4]
                else:
                    long_fastq = os.path.join(datadir, l[1])
                    r1_fastq = os.path.join(datadir, l[3])
                    r2_fastq = os.path.join(datadir, l[4])

                if (
                    os.path.isfile(long_fastq)
                    and l[2].isnumeric()
                    and os.path.isfile(r1_fastq)
                    and os.path.isfile(r2_fastq)
                ):
                    outDict[l[0]]["LR"] = long_fastq
                    outDict[l[0]]["MinChromLength"] = l[2]
                    outDict[l[0]]["R1"] = r1_fastq
                    outDict[l[0]]["R2"] = r2_fastq
                    outDict[l[0]]["TargetBases"] = int(l[2]) * subsample_depth
                else:
                    sys.stderr.write(
                        "\n"
                        f"    FATAL: Error parsing {csvFile}. One of \n"
                        f"    {long_fastq} or \n"
                        f"    {r1_fastq} or \n"
                        f"    {r2_fastq} \n"
                        f"    does not exist or  {l[2]} is not an integer. \n"
                        "    Check formatting, and that \n"
                        "    file names and file paths are correct.\n"
                        "\n"
                    )
                    sys.exit(1)
            else:
                sys.stderr.write(
                    "\n"
                    f"    FATAL: Error parsing {csvFile}. Line {l} \n"
                    f"    does not have 5 columns. \n"
                    f"    Please check the formatting of {csvFile}. \n"
                )
                sys.exit(1)
    return outDict


def parseSamples(csvfile, long_flag, subsample_depth, datadir):
    if os.path.isfile(csvfile) and long_flag is True:
        sampleDict = samplesFromCsvLong(csvfile, subsample_depth, datadir)
    elif os.path.isfile(csvfile) and long_flag is False:
        sampleDict = samplesFromCsvShort(csvfile, subsample_depth, datadir)
    else:
        sys.stderr.write(
            "\n"
            f"    FATAL: something is wrong. Likely {csvfile} is neither a file nor directory.\n"
            "\n"
        )
        sys.exit(1)

    # checks for dupes

    SAMPLES = list(sampleDict.keys())

    # Check for duplicates
    has_duplicates = len(SAMPLES) != len(set(SAMPLES))

    # error out if dupes
    if has_duplicates is True:
        sys.stderr.write(
            f"Duplicates found in the SAMPLES list in column 1 of {csvfile}.\n"
            f"Please check {csvfile} and give each sample a unique name!"
        )
        sys.exit(1)

    return sampleDict
