
rule polca:
    input:
        polypolish_fasta = os.path.join(dir.out.polypolish,"{sample}.fasta")
    output:
        fasta = os.path.join(dir.out.polca,"{sample}", "{sample}.fasta.PolcaCorrected.fa"),
        version = os.path.join(dir.out.versions, "{sample}", "polca_complete_masurca.version"),
        polca_input_fasta = os.path.join(dir.out.polca, "{sample}", "{sample}.fasta"),
    params:
        polca_input_fasta = "{sample}.fasta",
        dir = os.path.join(dir.out.polca, "{sample}"),
        reads = ' '.join(['"../../../../../'+ os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"), '../../../../../'+os.path.join(dir.out.fastp,"{sample}_2.fastq.gz"+'"')]),
        version = os.path.join(dir.out.versions, "{sample}", "polca_complete_masurca.version"),
        copy_fasta = os.path.join(dir.out.intermediate_assemblies, "{sample}", "{sample}_polca.fasta")
    conda:
        os.path.join(dir.env,'polca.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "polca", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "polca", "{sample}.log")
    shell:
        """
        # real struggle running polca honestly
        # need these workaround with '../../' etc
        cp {input.polypolish_fasta} {output.polca_input_fasta}
        cd {params.dir}
        polca.sh -a {params.polca_input_fasta}  -r {params.reads} -t {threads} 
        masurca --version > {params.version}
        """

rule copy_polca:
    input:
        fasta = os.path.join(dir.out.polca,"{sample}", "{sample}.fasta.PolcaCorrected.fa")
    output:
        copy_fasta = os.path.join(dir.out.intermediate_assemblies, "{sample}", "{sample}_polca.fasta")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        cp {input.fasta} {output.copy_fasta}
        """

rule polca_incomplete:
    input:
        polypolish_fasta = os.path.join(dir.out.polypolish_incomplete,"{sample}.fasta")
    output:
        fasta = os.path.join(dir.out.polca_incomplete,"{sample}", "{sample}.fasta.PolcaCorrected.fa"),
        version = os.path.join(dir.out.versions, "{sample}", "polca_complete_masurca.version"),
        polca_input_fasta = os.path.join(dir.out.polca_incomplete, "{sample}", "{sample}.fasta")
    params:
        polca_input_fasta = "{sample}.fasta",
        dir = os.path.join(dir.out.polca_incomplete, "{sample}"),
        reads = ' '.join(['"../../../../../'+ os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"), '../../../../../'+os.path.join(dir.out.fastp,"{sample}_2.fastq.gz"+'"')]),
        version = os.path.join(dir.out.versions, "{sample}", "polca_complete_masurca.version"),
        copy_fasta = os.path.join(dir.out.intermediate_assemblies, "{sample}", "{sample}_polca.fasta")
    conda:
        os.path.join(dir.env,'polca.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "polca_incomplete", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "polca_incomplete", "{sample}.log")
    shell:
        """
        # real struggle running polca honestly
        # need these workaround with '../../' etc
        CURR_DIR=PWD
        mkdir -p {params.intermediate_dir}
        cp {input.polypolish_fasta} {output.polca_input_fasta}
        cd {params.dir}
        polca.sh -a {params.polca_input_fasta}  -r {params.reads} -t {threads} 
        cd $CURR_DIR
        masurca --version > {params.version}
        cp {output.fasta} {params.copy_fasta} 

        """

rule copy_polca_incomplete:
    input:
        fasta = os.path.join(dir.out.polca_incomplete,"{sample}", "{sample}.fasta.PolcaCorrected.fa")
    output:
        copy_fasta = os.path.join(dir.out.intermediate_assemblies, "{sample}", "{sample}_polca.fasta")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        cp {input.fasta} {output.copy_fasta}
        """
