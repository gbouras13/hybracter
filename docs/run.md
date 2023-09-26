# Running `hybracter`

## Available Commands

* `hybracter hybrid`: Assemble multiple genomes from isolates that have long-reads and paired-end short reads.
* `hybracter hybrid-single`: Assembles a single genome from an isolate with long-reads and paired-end short reads. It takes similar parameters to [Unicycler](https://github.com/rrwick/Unicycler).
* `hybracter long`: Assemble multiple genomes from isolates that have long-reads only.
* `hybracter long-single`: Assembles a single genome from an isolate with long-reads only.
* `hybracter install`: Downloads and installs the required `plassembler` database.

## Input csv

`hybracter hybrid` and `hybracter long` require an input csv file to be specified with `--input`. No other inputs are required.

* This file requires no headers.
* Other than the reads, `hybracter` requires a value for a lower bound the minimum chromosome length for each isolate in base pairs. It must be an integer.
* `hybracter` will denote contigs about this value as chromosome(s) and if it can recover a chromosome, it will denote the isolate as complete.
* In practice, I suggest choosing 90% of the estimated chromosome size for this value.
* e.g. for _S. aureus_, I'd choose 2500000, _E. coli_, 4000000, _P. aeruginosa_ 5500000.

#### `hybracter hybrid`

* `hybracter hybrid` requires an input csv file with 5 columns. 
* Each row is a sample.
* Column 1 is the sample name you want for this isolate. 
* Column 2 is the long read fastq file.
* Column 3 is the minimum chromosome length for that sample.
* Column 4 is the R1 short read fastq file
* Column 5 is the R2 short read fastq file.

e.g.

```
s_aureus_sample1,sample1_long_read.fastq.gz,2500000,sample1_SR_R1.fastq.gz,sample1_SR_R2.fastq.gz
p_aeruginosa_sample2,sample2_long_read.fastq.gz,5500000,sample2_SR_R1.fastq.gz,sample2_SR_R2.fastq.gz
```

#### `hybracter long`

`hybracter long` also requires an input csv with no headers, but only 3 columns.

* `hybracter long` requires an input csv file with 3 columns. 
* Each row is a sample.
* Column 1 is the sample name you want for this isolate. 
* Column 2 is the long read fastq file.
* Column 3 is the minimum chromosome length for that sample.

e.g.

```
s_aureus_sample1,sample1_long_read.fastq.gz,2500000
p_aeruginosa_sample2,sample2_long_read.fastq.gz,5500000
```

## Usage

#### `hybracter install`

You will first need to install the `hybracter` databases.

```
hybracter install
```

Alternatively, can also specify a particular directory to store them - you will need to specify this with `-d <databases directory>` when you run `hybracter`.

```
hybracter install -d  <databases directory>
```

Once that is done, run `hybracter hybrid` or `hybracter long` as follows.

#### `hybracter hybrid`

```
hybracter hybrid -i <input.csv> -o <output_dir> -t <threads> 
```

* `hybracter hybrid` requires only a CSV file specified with `-i` or `--input`
* `--no_polca` will turn off POLCA polishing. **Must be used on MacOS.**
* Use `--min_length` to specify the minimum long-read length for Filtlong.
* Use `--min_quality` to specify the minimum long-read quality for Filtlong.
* You can specify a FASTA file containing contaminants with `--contaminants`. All long reads that map to contaminants will be filtered out.
  * You can specify Escherichia phage lambda (a common contaminant in Nanopore library preparation) using `--contaminants lambda`.
* `--skip_qc` will skip all read QC steps.
* You can change the `--medakaModel` (all available options are listed in `hybracter hybrid -h`)
* You can change the `--flyeModel` (all available options are listed in `hybracter hybrid -h`)


#### `hybracter hybrid-single`

```
hybracter hybrid-single -l <longread FASTQ> -1 <R1 short reads FASTQ> -2 <R2 short reads FASTQ> -s <sample name> -c <chromosome size> -o <output_dir> -t <threads>  [other arguments]
```

#### `hybracter long`

```
hybracter long -i <input.csv> -o <output_dir> -t <threads> [other arguments]
```

##### Other Arguments 

* `hybracter long` requires only a CSV file specified with `-i` or `--input`
* Use `--min_length` to specify the minimum long-read length for Filtlong.
* Use `--min_quality` to specify the minimum long-read quality for Filtlong.
* You can specify a FASTA file containing contaminants with `--contaminants`. All long reads that map to contaminants will be filtered.
  * You can specify Escherichia phage lambda (a common contaminant in Nanopore library preparation) using `--contaminants lambda`.
* `--skip_qc` will skip all read QC steps.
* You can change the `--medakaModel` (all available options are listed in `hybracter long -h`)
* You can change the `--flyeModel` (all available options are listed in `hybracter long -h`)


#### `hybracter long-single`

```
hybracter long-single -l <longread FASTQ> -s <sample name> -c <chromosome size>  -o <output_dir> -t <threads>  [other arguments]
```