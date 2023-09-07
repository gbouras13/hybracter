rule medaka_incomplete:
    input:
        os.path.join(INCOMPLETE_PRE_POLISH,"{sample}.fasta"),
        get_input_lr_fastqs
    output:
        directory(os.path.join(MEDAKA_INCOMPLETE,"{sample}")),
        os.path.join(MEDAKA_INCOMPLETE,"{sample}", "consensus.fasta")
    conda:
        os.path.join('..', 'envs','medaka.yaml')
    params:
        MEDAKA_MODEL
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time
    threads:
        config.resources.big.cpus
    shell:
        """
        medaka_consensus -i {input[1]} -d {input[0]} -o {output[0]} -m {params[0]}  -t {threads}
        """


