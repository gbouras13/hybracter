
rule polca:
    input:
        polypolish_fasta=os.path.join(dir.out.polypolish, "{sample}.fasta"),
    output:
        fasta=os.path.join(
            dir.out.polca, "{sample}", "{sample}.fasta.PolcaCorrected.fa"
        ),
        version=os.path.join(
            dir.out.versions, "{sample}", "polca_complete_masurca.version"
        ),
        polca_input_fasta=os.path.join(dir.out.polca, "{sample}", "{sample}.fasta"),
    params:
        polca_input_fasta="{sample}.fasta",
        dir=os.path.join(dir.out.polca, "{sample}"),
        reads=" ".join(
            [
                '"../../../../../'
                + os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
                "../../../../../"
                + os.path.join(dir.out.fastp, "{sample}_2.fastq.gz" + '"'),
            ]
        ),
        version=os.path.join(
            dir.out.versions, "{sample}", "polca_complete_masurca.version"
        ),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_polca.fasta"
        ),
    conda:
        os.path.join(dir.env, "polca.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "polca", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "polca", "{sample}.log"),
    shell:
        """
        CURR_DIR=$(pwd)
        cp {input.polypolish_fasta} {output.polca_input_fasta}
        cd {params.dir}
        if polca.sh -a {params.polca_input_fasta}  -r {params.reads} -t {threads}  ; then
            echo "POLCA SUCCEEDED"
            cd $CURR_DIR
        else
            echo "POLCA DID NOT DETECT ANY VARIANTS vs POLYPOLISH. COPYING THE POLYPOLISH FASTA."
            cd $CURR_DIR
            cp {input.polypolish_fasta} {output.fasta}
        fi  
        masurca --version > {params.version}
        cp {output.fasta} {params.copy_fasta} 
        """


rule compare_assemblies_polca_vs_polypolish:
    """
    compare assemblies 
    """
    input:
        reference=os.path.join(dir.out.polypolish, "{sample}.fasta"),
        assembly=os.path.join(
            dir.out.polca, "{sample}", "{sample}.fasta.PolcaCorrected.fa"
        ),
    output:
        diffs=os.path.join(dir.out.differences, "{sample}", "polca_vs_polypolish.txt"),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "compare_assemblies.py")


rule polca_incomplete:
    input:
        polypolish_fasta=os.path.join(dir.out.polypolish_incomplete, "{sample}.fasta"),
    output:
        fasta=os.path.join(
            dir.out.polca_incomplete, "{sample}", "{sample}.fasta.PolcaCorrected.fa"
        ),
        version=os.path.join(
            dir.out.versions, "{sample}", "polca_complete_masurca.version"
        ),
        polca_input_fasta=os.path.join(
            dir.out.polca_incomplete, "{sample}", "{sample}.fasta"
        ),
    params:
        polca_input_fasta="{sample}.fasta",
        dir=os.path.join(dir.out.polca_incomplete, "{sample}"),
        reads=" ".join(
            [
                '"../../../../../'
                + os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
                "../../../../../"
                + os.path.join(dir.out.fastp, "{sample}_2.fastq.gz" + '"'),
            ]
        ),
        version=os.path.join(
            dir.out.versions, "{sample}", "polca_complete_masurca.version"
        ),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_polca.fasta"
        ),
    conda:
        os.path.join(dir.env, "polca.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "polca_incomplete", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "polca_incomplete", "{sample}.log"),
    shell:
        """
        CURR_DIR=$(pwd)
        cp {input.polypolish_fasta} {output.polca_input_fasta}
        cd {params.dir}
        if polca.sh -a {params.polca_input_fasta}  -r {params.reads} -t {threads}  ; then
            echo "POLCA SUCCEEDED"
            cd $CURR_DIR
        else
            echo "POLCA DID NOT DETECT ANY VARIANTS vs POLYPOLISH. COPYING THE POLYPOLISH FASTA."
            cd $CURR_DIR
            cp {input.polypolish_fasta} {output.fasta}
        fi  
        masurca --version > {params.version}
        cp {output.fasta} {params.copy_fasta} 
        """
