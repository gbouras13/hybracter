
# input function for the rule aggregate
def aggregate_long_read_polish_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[0].open() as f:
        if f.read().strip() == "C":
            return os.path.join(MEDAKA_RD_2,"{sample}", "consensus.fasta")
        else:
            return os.path.join(MEDAKA_INCOMPLETE,"{sample}", "consensus.fasta")


### from the long_read_polishing 


rule aggregate_long_read_polish_input:
    input:
        aggregate_long_read_polish_input
    output:
        os.path.join(AGGREGATE_LONG_READ_POLISH,"{sample}.txt")
    shell:
        """
        touch {output}
        """

rule aggr_long_read_polish:
    """Aggregate."""
    input:
        expand(os.path.join(AGGREGATE_LONG_READ_POLISH,"{sample}.txt"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_long_read_polish.txt")
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    shell:
        """
        touch {output[0]}
        """


######### short reads ############

# input function for the rule aggregate
def aggregate_short_read_polish_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[0].open() as f:
        if f.read().strip() == "C":
            return os.path.join(POLYPOLISH,"{sample}.fasta")
        else:
            return os.path.join(POLYPOLISH_INCOMPLETE,"{sample}.fasta")


### from the short_read_polishing 


rule aggregate_short_read_polish_input:
    input:
        aggregate_short_read_polish_input
    output:
        os.path.join(AGGREGATE_SHORT_READ_POLISH,"{sample}.txt")
    shell:
        """
        touch {output}
        """


rule aggr_long_read_polish:
    """Aggregate."""
    input:
        expand(os.path.join(AGGREGATE_SHORT_READ_POLISH,"{sample}.txt"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_short_read_polish.txt")
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    shell:
        """
        touch {output[0]}
        """

####### polca ##########


# input function for the rule aggregate polca
def aggregate_polca_polish_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[0].open() as f:
        if f.read().strip() == "C":
            return os.path.join(POLCA,"{sample}.fasta.PolcaCorrected.fa")
        else:
            return os.path.join(POLCA_INCOMPLETE,"{sample}.fasta.PolcaCorrected.fa")


### from the short_read_polishing 
rule aggregate_polca_polish_input:
    input:
        aggregate_polca_polish_input
    output:
        os.path.join(AGGREGATE_POLCA_POLISH,"{sample}.txt")
    shell:
        """
        touch {output}
        """


rule aggr_polca:
    """Aggregate."""
    input:
        expand(os.path.join(AGGREGATE_POLCA_POLISH,"{sample}.txt"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_polca.txt")
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    shell:
        """
        touch {output[0]}
        """
