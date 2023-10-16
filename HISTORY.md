# History

## v0.2.0 (16 October 2023)

* Replacement of POLCA by [pypolca](https://github.com/gbouras13/pypolca) as it will be easier to integrate, install and maintain, and allows for POLCA use on MacOS.
* Adds `--dnaapler_custom_db`
* @simone-pignotti fixed some errors with fastp redirection and threads.
* Adds various small improvements (mem for PBS use, various other params and paths, some of the conda envs).

## v0.1.2 (6 October 2023)

* Fixes [bugs](https://github.com/gbouras13/hybracter/issues/13) with bwa index creation and typos in some output files.
* Thanks  @simone-pignotti for detecting and fixing [it](https://github.com/gbouras13/hybracter/pull/14)


## v0.1.1 (5 October 2023)

* Fixes a small [bug](https://github.com/gbouras13/hybracter/issues/9) with samples.smk
* Thanks @npbhavya and @simone-pignotti for detecting it.

## v0.1.0 (28 September 2023)

* Initial release.

```
 _           _                    _            
| |__  _   _| |__  _ __ __ _  ___| |_ ___ _ __ 
| '_ \| | | | '_ \| '__/ _` |/ __| __/ _ \ '__|
| | | | |_| | |_) | | | (_| | (__| ||  __/ |   
|_| |_|\__, |_.__/|_|  \__,_|\___|\__\___|_|   
       |___/


Usage: hybracter [OPTIONS] COMMAND [ARGS]...

  For more options, run: hybracter command --help

Options:
  -h, --help  Show this message and exit.

Commands:
  install        Downloads and installs the plassembler database
  hybrid         Run hybracter with hybrid long and paired end short reads
  hybrid-single  Run hybracter hybrid on 1 isolate
  long           Run hybracter with only long reads
  long-single    Run hybracter long on 1 isolate
  test-hybrid    Test hybracter hybrid
  test-long      Test hybracter long
  config         Copy the system default config file
  citation       Print the citation(s) for hybracter
  version        Print the version for hybracter
```