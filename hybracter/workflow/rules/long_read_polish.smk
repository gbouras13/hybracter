rule medaka_round_1:
    input:
        os.path.join(CHROMOSOME_PRE_POLISH,"{sample}.fasta"),
        get_input_lr_fastqs
    output:
        directory(os.path.join(MEDAKA_RD_1,"{sample}")),
        os.path.join(MEDAKA_RD_1,"{sample}", "consensus.fasta")
    conda:
        os.path.join('..', 'envs','medaka.yaml')
    params:
        MEDAKA_MODEL
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    threads:
        BigJobCpu
    shell:
        """
        medaka_consensus -i {input[1]} -d {input[0]} -o {output[0]} -m {params[0]}  -t {threads}
        """

rule dnaapler:
    input:
        os.path.join(MEDAKA_RD_1,"{sample}", "consensus.fasta")
    output:
        os.path.join(DNAAPLER, "{sample}", "{sample}_reoriented.fasta")
    threads:
        1
    conda:
        os.path.join('..', 'envs','dnaapler.yaml')
    params:
        os.path.join(DNAAPLER, "{sample}")
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    threads:
        BigJobCpu
    shell:
        """
        dnaapler chromosome -i {input[0]} -o {params[0]} -p {wildcards.sample} -t {threads} -f
        """

rule medaka_round_2:
    input:
        os.path.join(DNAAPLER, "{sample}", "{sample}_reoriented.fasta"),
        get_input_lr_fastqs
    output:
        directory(os.path.join(MEDAKA_RD_2,"{sample}")),
        os.path.join(MEDAKA_RD_2,"{sample}", "consensus.fasta")
    conda:
        os.path.join('..', 'envs','medaka.yaml')
    params:
        MEDAKA_MODEL
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    threads:
        BigJobCpu
    shell:
        """
        medaka_consensus -i {input[1]} -d {input[0]} -o {output[0]} -m {params[0]}  -t {threads}
        """

rule dnaapler_no_polish:
    input:
        os.path.join(CHROMOSOME_PRE_POLISH,"{sample}.fasta")
    output:
        os.path.join(DNAAPLER_NO_POLISH, "{sample}", "{sample}_reoriented.fasta")
    threads:
        1
    conda:
        os.path.join('..', 'envs','dnaapler.yaml')
    params:
        os.path.join(DNAAPLER_NO_POLISH, "{sample}")
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    threads:
        BigJobCpu
    shell:
        """
        dnaapler chromosome -i {input[0]} -o {params[0]} -p {wildcards.sample} -t {threads} -f
        """