#!/bin/bash

conda activate badread

badread simulate  --reference assembly.fasta --quantity 60x   --seed 43 --length 10000,10000| gzip > test_long_reads.fastq.gz

conda activate iss

iss generate --genomes assembly.fasta --cpus 8 --model novaseq --compress  --output test_short_reads --n_reads 0.05M


