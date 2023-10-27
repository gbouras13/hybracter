#!/bin/bash

conda activate badread
# Badread v0.3.0 nanopore2023

badread simulate  --reference assembly.fasta --quantity 60x   --seed 43 --length 10000,10000| gzip > test_long_reads_tmp.fastq.gz
badread simulate --reference C222_plasmid.fasta --quantity 120x --length 800,500 | gzip > C222_reads.fastq.gz

# to get more C222 long reads to test hybracter long now plassembler should work :)

cat C222_reads.fastq.gz  test_long_reads_tmp.fastq.gz > test_long_reads.fastq.gz

conda activate iss

iss generate --genomes assembly.fasta --cpus 8 --model novaseq --compress  --output test_short_reads --n_reads 0.05M


