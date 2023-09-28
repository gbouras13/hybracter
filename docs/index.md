# `hybracter`

`hybracter` is modern hybrid (and long-only) bacterial assembly pipeline for many isolates using Snakemake and [Snaketool](https://github.com/beardymcjohnface/Snaketool).

## Overview

`hybracter` is designed for assembling many bacterial isolate genomes using the embarassingly parallel power of HPC and Snakemake profiles. It is designed for applications where you have a number of isolates with Oxford Nanopore Technologies (ONT) long reads and optionally matched paired-end short reads for polishing.

`hybracter` is desined to straddle the fine line between being as fully feature-rich as possible with as much information as you need to decide upon the best assembly, while also being a one-line automated program. Perfect for lazy people like me. 

`hybracter` is largely based off Ryan Wick's [magnificent tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial) and associated [paper](https://doi.org/10.1371/journal.pcbi.1010905). `hybracter` differs in that it adds some additional steps regarding targetted plasmid assembly with [plassembler](https://github.com/gbouras13/plassembler), contig reorientation with [dnaapler](https://github.com/gbouras13/dnaapler) and extra polishing and statistical summaries.

![Image](hybracter.png)


