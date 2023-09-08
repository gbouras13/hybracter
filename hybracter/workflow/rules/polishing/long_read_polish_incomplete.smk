rule medaka_incomplete:
    input:
        os.path.join(dir.out.incomp_pre_polish,"{sample}.fasta"),
        get_input_lr_fastqs
    output:
        directory(os.path.join(dir.out.medaka_incomplete,"{sample}")),
        os.path.join(dir.out.medaka_incomplete,"{sample}", "consensus.fasta")
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


