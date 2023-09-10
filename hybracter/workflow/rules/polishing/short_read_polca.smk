
rule polca:
    input:
        polypolish_fasta = os.path.join(dir.out.polypolish,"{sample}.fasta")
    output:
        fasta = os.path.join(dir.out.polca,"{sample}", "{sample}.fasta.PolcaCorrected.fa"),
        version = os.path.join(dir.out.versions, "{sample}", "polca_complete_masurca.version"),
        polca_input_fasta = os.path.join(dir.out.polca, "{sample}", "{sample}.fasta")
    params:
        #r1 = os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"),
        #r2 = os.path.join(dir.out.fastp,"{sample}_2.fastq.gz"),
        dir = os.path.join(dir.out.polca, "{sample}"),
        reads = ' '.join(['"'+ os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"), os.path.join(dir.out.fastp,"{sample}_2.fastq.gz"+'"')])
    conda:
        os.path.join(dir.env,'polca.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        # real struggle running polca honestly
        cp {input.polypolish_fasta} {output.polca_input_fasta}
        cd {params.dir}
        polca.sh -a {output.input_fasta}  -r {params.reads} -t {threads}
        masurca --version > {output.version}
        """


rule polca_incomplete:
    input:
        polypolish_fasta = os.path.join(dir.out.polypolish_incomplete,"{sample}.fasta")
    output:
        fasta = os.path.join(dir.out.polca_incomplete,"{sample}", "{sample}.fasta.PolcaCorrected.fa"),
        version = os.path.join(dir.out.versions, "{sample}", "polca_complete_masurca.version"),
        polca_input_fasta = os.path.join(dir.out.polca_incomplete, "{sample}", "{sample}.fasta")
    params:
        #r1 = os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"),
        #r2 = os.path.join(dir.out.fastp,"{sample}_2.fastq.gz"),
        dir = os.path.join(dir.out.polca_incomplete, "{sample}"),
        reads = ' '.join(['"../../../../'+ os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"), '../../../../'+os.path.join(dir.out.fastp,"{sample}_2.fastq.gz"+'"')])
    conda:
        os.path.join(dir.env,'polca.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        # real struggle running polca honestly
        cp {input.polypolish_fasta} {output.polca_input_fasta}
        cd {params.dir}
        polca.sh -a {output.polca_input_fasta}  -r {params.reads} -t {threads}
        masurca --version > {output.version}
        """






