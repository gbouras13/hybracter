# Advanced Configuration and `hybracter`

The `hybracter` configuration file contains settings related to resource memory, cpu and time consumption for various stages of the pipeline. There are also specific settings for optional host depletion, but it is unlikely you will want to changes these!

The default settings ensure the resource intensive parts of `hybracter` (namely, the assembly, reorientation and polishing steps) will use available threads when run on a single isolate on a standard 4-core 8-thread laptop.

For example the default settings are:

```yaml
resources:
    big:                   #### For larger jobs like assembly and polishing (Medaka/Polypolish/pypolca)
        cpu: 16               # number of threads used
        mem: 32000            # Mem for the job in MB
        time: '23:59:00'      # time in hh:mm:ss
    med:                   #### For medium jobs like reorientation (dnaapler)
        cpu: 8    
        mem: 16000    
        time: '08:00:00'  
    sml:                   #### For small jobs such as python scripting 
        cpu: 1             #### It is not recommended to change this
        mem: 4000    
        time: '00:00:05'  
```

Users can copy and modify the default config file to suit their own needs by running `hybracter config`. 

```bash
hecatomb config
```

This will copy the default config file to the `hybracter_out` directory as `config.yaml`. You can modify this as you please to suit your system, and you can specify it using the `--configfile` parameter.

```bash
hybracter hybrid --configfile my_modified_config.yaml
```

For example, if you have 16 threads and 64 GB available on your system and you have multiple isolates, you may want to reduce the big job cpu from 16 to 8 or 4 and the memory from 32000 MB to 16000MB. This is because it will probably to be faster for your overall sample set for hybracter to run Flye twice simultaneously with 8 threads than once on 16 threads.

You could achieve this as follows:

```yaml
resources:
    big:                   #### For larger jobs like assembly and polishing (Medaka/Polypolish/pypolca)
        cpu: 8               # number of threads used
        mem: 16000            # Mem for the job in MB
        time: '23:59:00'      # time in hh:mm:ss
    med:                   #### For medium jobs like reorientation (dnaapler)
        cpu: 4    
        mem: 8000    
        time: '08:00:00'  
    sml:                   #### For small jobs such as python scripting 
        cpu: 1             #### It is not recommended to change this
        mem: 4000    
        time: '00:00:05'  
```