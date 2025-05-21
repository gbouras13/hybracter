"""
Function for parsing the config and identifying samples and read files
"""

from itertools import chain
import sys

"""
check_datadir function
"""


def check_datadir(datadir: str):
    # if only one (hybracter long)
    if "," not in datadir:
        return datadir, None
    else:
        # if more than one (hybracter hybrid)
        parts = datadir.split(",")

        if len(parts) > 2:
            sys.exit(
                f"Error: You must specify a maximum of two comma separated input directories with --datadir, the first being where the long-read FASTQs are, the second the short-read.\n You have specified --datadir {datadir}.\n Please check this."
            )

        # If exactly one comma, assign the two parts
        datadirlong, datadirshort = parts
        return datadirlong, datadirshort


"""
long
"""


def samplesFromCsvLong(csvFile, subsample_depth, datadir, min_depth, auto):
    """
    Read samples and files from a CSV Long Read Only
    3 cols
    1 = sample
    2 = Long read
    3 = MinChromLength
    IF auto - only need 2

    chrom_size, target_bases and min_depth handled with functions if --auto
    """
    outDict = {}
    with open(csvFile, "r") as csv:
        for line in csv:
            l = line.strip().split(",")
            if auto:
                if len(l) == 2:
                    outDict[l[0]] = {}
                    # where datafir isn't specified
                    if datadir is None:
                        long_fastq = l[1]
                    else:
                        datadirlong, datadirshort = check_datadir(datadir)
                        long_fastq = os.path.join(datadirlong, l[1])
                    if os.path.isfile(long_fastq):
                        outDict[l[0]]["LR"] = long_fastq
                    else:
                        sys.stderr.write(
                            "\n"
                            f"    FATAL: Error parsing {csvFile}. {long_fastq} \n"
                            f"    does not exist. \n"
                            "    Check formatting, and that \n"
                            "    file names and file paths are correct.\n"
                            "\n"
                        )
                        sys.exit(1)
                else:
                    sys.stderr.write(
                        "\n"
                        f"    FATAL: Error parsing {csvFile}. Line {l} \n"
                        f"    does not have 2 columns. \n"
                        f"    Please check the formatting of {csvFile}. \n"
                    )
                    sys.exit(1)

            else:
                if len(l) == 3:
                    outDict[l[0]] = {}
                    # where datafir isn't specified
                    if datadir is None:
                        long_fastq = l[1]
                    else:
                        datadirlong, datadirshort = check_datadir(datadir)
                        long_fastq = os.path.join(datadirlong, l[1])
                    if os.path.isfile(long_fastq) and l[2].isnumeric():
                        outDict[l[0]]["LR"] = long_fastq
                        outDict[l[0]]["MinChromLength"] = l[2]
                        outDict[l[0]]["TargetBases"] = int(l[2]) * subsample_depth
                        outDict[l[0]]["MinBases"] = int(l[2]) * min_depth
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


def samplesFromCsvShort(csvFile, subsample_depth, datadir, min_depth, auto):
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
            if auto:
                if len(l) == 4:
                    outDict[l[0]] = {}
                    # where datafir isn't specified
                    if datadir is None:
                        long_fastq = l[1]
                        r1_fastq = l[2]
                        r2_fastq = l[3]
                    else:
                        # if the user supplies a datadir
                        datadirlong, datadirshort = check_datadir(datadir)
                        long_fastq = os.path.join(datadirlong, l[1])
                        # all fastqs in 1 dir
                        if datadirshort is None:
                            r1_fastq = os.path.join(datadirlong, l[2])
                            r2_fastq = os.path.join(datadirlong, l[3])
                        # separate dirs
                        else:
                            r1_fastq = os.path.join(datadirshort, l[2])
                            r2_fastq = os.path.join(datadirshort, l[3])

                    if (
                        os.path.isfile(long_fastq)
                        and os.path.isfile(r1_fastq)
                        and os.path.isfile(r2_fastq)
                    ):
                        outDict[l[0]]["LR"] = long_fastq
                        outDict[l[0]]["R1"] = r1_fastq
                        outDict[l[0]]["R2"] = r2_fastq
                    else:
                        sys.stderr.write(
                            "\n"
                            f"    FATAL: Error parsing {csvFile}. One of \n"
                            f"    {long_fastq} or \n"
                            f"    {r1_fastq} or \n"
                            f"    {r2_fastq} \n"
                            f"    does not exist. \n"
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

            else:
                if len(l) == 5:
                    outDict[l[0]] = {}
                    # where datafir isn't specified
                    if datadir is None:
                        long_fastq = l[1]
                        r1_fastq = l[3]
                        r2_fastq = l[4]
                    else:
                        # if the user supplies a datadir
                        datadirlong, datadirshort = check_datadir(datadir)
                        long_fastq = os.path.join(datadirlong, l[1])
                        # all fastqs in 1 dir
                        if datadirshort is None:
                            r1_fastq = os.path.join(datadirlong, l[3])
                            r2_fastq = os.path.join(datadirlong, l[4])
                        # separate dirs
                        else:
                            r1_fastq = os.path.join(datadirshort, l[3])
                            r2_fastq = os.path.join(datadirshort, l[4])
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
                        outDict[l[0]]["MinBases"] = int(l[2]) * min_depth
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

def samplesFromCsvAutomatic(csvFile, subsample_depth, datadir, min_depth, auto):
    """
    Read samples and files from a CSV in automatic mode
    Accepts mixed format with either:
    - 3 cols (sample, long_read, MinChromLength) for long-read only (if auto=False)
    - 2 cols (sample, long_read) for long-read only (if auto=True)
    - 5 cols (sample, long_read, MinChromLength, short_read_1, short_read_2) for hybrid (if auto=False)
    - 4 cols (sample, long_read, short_read_1, short_read_2) for hybrid (if auto=True)
    Empty values for short reads indicate long-read only samples
    """
    outDict = {}
    with open(csvFile, "r") as csv:
        for line in csv:
            l = line.strip().split(",")
            # Initialize sample with default None values for R1/R2
            outDict[l[0]] = {
                "R1": None,
                "R2": None
            }

            # Handle paths based on datadir
            if datadir is None:
                long_fastq = l[1]
            else:
                datadirlong, datadirshort = check_datadir(datadir)
                long_fastq = os.path.join(datadirlong, l[1])

            # Validate and add long read path
            if not os.path.isfile(long_fastq):
                sys.stderr.write(f"\nFATAL: Long read file {long_fastq} does not exist\n")
                sys.exit(1)
            
            outDict[l[0]]["LR"] = long_fastq

            # Handle MinChromLength based on auto parameter
            if not auto:
                chrom_length_idx = 2
                if not l[chrom_length_idx].isnumeric():
                    sys.stderr.write(f"\nFATAL: MinChromLength {l[chrom_length_idx]} is not a number\n")
                    sys.exit(1)
                outDict[l[0]]["MinChromLength"] = l[chrom_length_idx]
                outDict[l[0]]["TargetBases"] = int(l[chrom_length_idx]) * subsample_depth
                outDict[l[0]]["MinBases"] = int(l[chrom_length_idx]) * min_depth

            # Handle short reads if present and not empty
            if len(l) >= 4:  # Check if we have enough columns for short reads
                r1_idx = 2 if auto else 3
                r2_idx = 3 if auto else 4
                
                # Only process if both R1 and R2 are non-empty
                if l[r1_idx].strip() and l[r2_idx].strip():
                    if datadir is None:
                        r1_fastq = l[r1_idx]
                        r2_fastq = l[r2_idx]
                    else:
                        if datadirshort is None:
                            r1_fastq = os.path.join(datadirlong, l[r1_idx])
                            r2_fastq = os.path.join(datadirlong, l[r2_idx])
                        else:
                            r1_fastq = os.path.join(datadirshort, l[r1_idx])
                            r2_fastq = os.path.join(datadirshort, l[r2_idx])

                    if not (os.path.isfile(r1_fastq) and os.path.isfile(r2_fastq)):
                        sys.stderr.write(
                            f"\nFATAL: Short read files {r1_fastq} or {r2_fastq} do not exist\n"
                        )
                        sys.exit(1)

                    outDict[l[0]]["R1"] = r1_fastq
                    outDict[l[0]]["R2"] = r2_fastq

    return outDict


def parseSamples(csvfile, long_flag, subsample_depth, datadir, min_depth, auto, automatic_flag=False):
    """
    Parse samples from CSV file
    automatic_flag: If True, use automatic mode to detect sample types
    """
    if not os.path.isfile(csvfile):
        sys.stderr.write(f"\nFATAL: {csvfile} is not a file\n")
        sys.exit(1)
        
    if automatic_flag:
        return samplesFromCsvAutomatic(csvfile, subsample_depth, datadir, min_depth, auto)
    elif long_flag:
        return samplesFromCsvLong(csvfile, subsample_depth, datadir, min_depth, auto)
    else:
        return samplesFromCsvShort(csvfile, subsample_depth, datadir, min_depth, auto)

# helper function to see if a sample should be hybrid or long in the automatic mode
def filter_samples_by_workflow_type(samples_dict):
    """
    Filter samples into hybrid and long-read-only based on presence of short read data
    Input: Dictionary from parseSamples
    Returns: Dict with 'hybrid' and 'long' sample lists
    """
    hybrid_samples = []
    long_samples = []

    for sample_id, sample_data in samples_dict.items():
        # Check if R1 and R2 exist and are not None
        if sample_data.get("R1") is not None and sample_data.get("R2") is not None:
            hybrid_samples.append(sample_id)
        else:
            long_samples.append(sample_id)

    return {
        'hybrid': hybrid_samples,
        'long': long_samples
    }

