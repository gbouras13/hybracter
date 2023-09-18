"""
Defines all functions used in hybracter
"""


# define functions
# get long reads
def get_input_lr_fastqs(wildcards):
    return dictReads[wildcards.sample]["LR"]


# get min chrom length (define chrom size)
def getMinChromLength(wildcards):
    return dictReads[wildcards.sample]["MinChromLength"]


def get_input_r1(wildcards):
    return dictReads[wildcards.sample]["R1"]


def get_input_r2(wildcards):
    return dictReads[wildcards.sample]["R2"]
