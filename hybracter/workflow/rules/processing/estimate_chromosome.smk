"""
@richardstoeckl with a nice implementation
https://github.com/gbouras13/hybracter/issues/90
"""
rule sanitiseInputReads:
    input:
        fastq=get_input_lr_fastqs(sample=wildcards.sample)
    output:
        fastq=temp(os.path.join(dir.out.qc, "{sample}_sanitised.fastq.gz"))
    conda:
        os.path.join(dir.env, "seqkit.yaml")
    threads:
        1
    shell:
        """
        seqkit replace -p "\t" -r '_' {input.fastq} -o {output.fastq}
        """

checkpoint kmc:
    input:
        fastq = rules.sanitiseInputReads.output.fastq,
    output:
        kmcLOG = os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt"),
        dir = temp(directory(os.path.join(dir.out.kmc,"{sample}"))),
    conda:
        os.path.join(dir.env, "kmc.yaml")
    threads:
        config.resources.med.cpu
    shell:
        "kmc -k21 -ci10 -t{threads} -v -fq {input.fastq} {wildcards.sample} {output.dir} > {output.kmcLOG} 2>&1"


