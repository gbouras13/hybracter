rule compile_c_program:
    output:
        os.path.join(dir.env, "ALE-master", "src", "ALE")
    params: 
        rule_dir = dir.env 
    conda:
        os.path.join(dir.env,'make.yaml')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.med.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        cd {params.rule_dir}
        wget https://github.com/sc932/ALE/archive/master.zip
        unzip master.zip 
        cd ALE-master/src
        make
        """
