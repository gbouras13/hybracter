"""

Stores all the aggregation rules and checkpoints post assembly - due to the split between complete and incomplete

"""

"""
long read polishing
"""

# input function for the rule aggregate
def aggregate_long_read_polish_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[0].open() as f:
        if f.read().strip() == "C":
            return os.path.join(dir.out.medaka_rd_2,"{sample}", "consensus.fasta")
        else: # if incomplete  
                return os.path.join(dir.out.medaka_incomplete,"{sample}.fasta")

### from the long_read_polishing 

rule aggregate_long_read_polish_input:
    input:
        aggregate_long_read_polish_input
    output:
        os.path.join(dir.out.aggr_lr_polish ,"{sample}.txt")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        echo {input[0]}
        touch {output}
        """

rule aggr_long_read_polish_flag:
    """Aggregate."""
    input:
        expand(os.path.join(dir.out.aggr_lr_polish ,"{sample}.txt"), sample = SAMPLES)
    output:
        flag = os.path.join(dir.out.flags, "aggr_long_read_polish.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """


"""
Short read polishing
"""

# input function for the rule aggregate
def aggregate_short_read_polish_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[0].open() as f:
        if f.read().strip() == "C":
            return os.path.join(dir.out.polypolish,"{sample}.fasta")
        else:
            return os.path.join(dir.out.polypolish_incomplete,"{sample}.fasta")


### from the short_read_polishing 


rule aggregate_short_read_polish_input:
    input:
        aggregate_short_read_polish_input
    output:
        os.path.join(dir.out.aggr_sr_polish ,"{sample}.txt")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output}
        """


rule aggr_short_read_polish_flag:
    """Aggregate."""
    input:
        expand(os.path.join(dir.out.aggr_sr_polish,"{sample}.txt"), sample = SAMPLES)
    output:
        flag = os.path.join(dir.out.flags, "aggr_short_read_polish.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
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
            return os.path.join(dir.out.polca,"{sample}", "{sample}.fasta.PolcaCorrected.fa")
        else:
            return os.path.join(dir.out.polca_incomplete,"{sample}.fasta.PolcaCorrected.fa")


### from the short_read_polishing 
rule aggregate_polca_polish_input:
    input:
        aggregate_polca_polish_input
    output:
        os.path.join(dir.out.aggr_polca_polish,"{sample}.txt")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output}
        """


rule aggr_polca_flag:
    """Aggregate."""
    input:
        expand(os.path.join(dir.out.aggr_polca_polish,"{sample}.txt"), sample = SAMPLES)
    output:
        flag = os.path.join(dir.out.flags, "aggr_polca.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """
