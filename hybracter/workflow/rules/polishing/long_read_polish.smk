rule medaka_round_1:
    input:
        os.path.join(dir.out.chrom_pre_polish,"{sample}.fasta"),
        get_input_lr_fastqs
    output:
        directory(os.path.join(dir.out.medaka_rd_1 ,"{sample}")),
        os.path.join(dir.out.medaka_rd_1 ,"{sample}", "consensus.fasta")
    conda:
        os.path.join(dir.env,'medaka.yaml')
    params:
        MEDAKA_MODEL
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time
    threads:
        config.resources.big.cpu
    shell:
        """
        medaka_consensus -i {input[1]} -d {input[0]} -o {output[0]} -m {params[0]}  -t {threads}
        """

rule dnaapler:
    input:
        os.path.join(dir.out.medaka_rd_1 ,"{sample}", "consensus.fasta")
    output:
        os.path.join(dir.out.dnaapler , "{sample}", "{sample}_reoriented.fasta")
    conda:
        os.path.join(dir.env,'dnaapler.yaml')
    params:
        os.path.join(dir.out.dnaapler, "{sample}")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        dnaapler chromosome -i {input[0]} -o {params[0]} -p {wildcards.sample} -t {threads} -f
        """

rule medaka_round_2:
    input:
        os.path.join(dir.out.dnaapler , "{sample}", "{sample}_reoriented.fasta"),
        get_input_lr_fastqs
    output:
        directory(os.path.join(dir.out.medaka_rd_2,"{sample}")),
        os.path.join(dir.out.medaka_rd_2,"{sample}", "consensus.fasta")
    conda:
        os.path.join(dir.env,'medaka.yaml')
    params:
        MEDAKA_MODEL
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time
    threads:
        config.resources.big.cpu
    shell:
        """
        medaka_consensus -i {input[1]} -d {input[0]} -o {output[0]} -m {params[0]}  -t {threads}
        """

