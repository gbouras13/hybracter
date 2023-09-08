
rule polca:
    input:
        os.path.join(dir.out.polypolish_incomplete,"{sample}.fasta")
    output:
        os.path.join(dir.out.polca,"{sample}.fasta.PolcaCorrected.fa")
    params:
        os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"),
        os.path.join(dir.out.fastp,"{sample}_2.fastq.gz"),
        dir.out.polca
    conda:
        os.path.join(dir.env,'polca.yaml')
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time
    threads:
        config.resources.big.cpu
    shell:
        """
        # real struggle running polca honestly
        cd {params[2]}
        polca.sh -a ../../../../{input[0]}  -r "../../../../{params[0]} ../../../../{params[1]}" -t {threads}
        """


rule polca_incomplete:
    input:
        os.path.join(dir.out.polypolish_incomplete,"{sample}.fasta")
    output:
        os.path.join(dir.out.polca_incomplete,"{sample}.fasta.PolcaCorrected.fa")
    params:
        os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"),
        os.path.join(dir.out.fastp,"{sample}_2.fastq.gz"),
        dir.out.polca_incomplete
    conda:
        os.path.join(dir.env,'polca.yaml')
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time
    threads:
        config.resources.big.cpu
    shell:
        """
        # real struggle running polca honestly
        cd {params[2]}
        polca.sh -a ../../../../{input[0]}  -r "../../../../{params[0]} ../../../../{params[1]}" -t {threads}
        """






