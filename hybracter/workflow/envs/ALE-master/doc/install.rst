.. ALE documentation master file, created by
   sphinx-quickstart on Fri Dec 16 21:11:32 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============

Getting Started:
   #. `Installation`_
   #. `Running ALE`_
   #. `Plotting the Output`_

Installation
------------

Download the latest source::

   $ git clone git@github.com:sc932/ALE.git

Enter the directory and run **make**::

   $ cd ALE/src
   $ make

To build this documentation (optional) run::

   $ cd ../doc
   $ make html

This documentation can now browsed from *../doc/_build/html/index.html*

Running ALE
-----------

After installation we can run ALE from the ALE src directory::

   $ ./ALE
   Usage: ./ALE [-options] readSorted.[s|b]am assembly.fasta[.gz] ALEoutput.ale
      Options:
      -h : print out help

From this point if you have a bam/sam file and an assembly you can run ALE directly. What follows is an example of how to create these synthetically from scratch and produce figures similar to that of the paper.

ALE Synthetic Example From Scratch
++++++++++++++++++++++++++++++++++

We will make some synthetic reads from the first part of E.Coli K12

First we download the genome for Escherichia_coli_K_12_substr__DH10B::

   $ cd example
   $ wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid20079/CP000948.fna

Then we extract the first 350k bases::

   $ head -50001 CP000948.fna > Ecoli_first350k.fna

Now we generate a set of 2 million synthetic reads (for more info see documentation)::

   $ ../synthReadGen -ip 1.0 -nr 2000000 -ps 10 -b Ecoli_first350k.fna Ecoli_first350k.fastq

At this point we can run **bowtie** to get an initial mapping::

   $ bowtie-build Ecoli_first350k.fna Ecoli_first350k
   $ bowtie -t -I 0 -X 300 --fr -a -l 10 -v 1 -e 300 -S --threads 2 Ecoli_first350k -1 part1_Ecoli_first350k.fastq  -2 part2_Ecoli_first350k.fastq Ecoli_first350k.map.sam

Now we have the required .sam and .fna files to run ALE::

   $ ../ALE Ecoli_first350k.map.sam Ecoli_first350k.fna Ecoli_first350k.ale

This results in a file Ecoli_first350k.ale in the following format::

   $ head -6 Ecoli_first350k.ale
   # Reference: gi|170079663|ref|NC_010473.1| 350000
   # contig position depth ln(depthLike) ln(placeLike) ln(kmerLike) ln(totalLike)
   0 1.000000 -60.000000 0.194888 -5.760798 -65.565910
   0 1 3.000000 -60.000000 0.466271 -5.608334 -65.142063
   0 2 5.000000 -60.000000 0.010585 -5.541655 -65.531071
   0 3 12.000000 -60.000000 -0.057731 -5.380759 -65.438491

We can use this information in its raw format or plot it using *plotter3.py*

For a more complete example see :mod:`image_maker`

Plotting the Output
-------------------

The authors recommend using IGV to view the output.

http://www.broadinstitute.org/igv/

Just import the assembly, bam and ALE scores. You can convert the .ale file to a set of .wig files with ale2wiggle.py and IGV can read those directly.  Depending on your genome size you may want to convert the .wig files to the BigWig format.

Here we show how to use the built in debugging plotter. Note, this plotter is no longer in active development and is not supported. We recommend using IGV.

We can invoke the plotter by running::

   $ ./plotter3.py ALEoutput.ale

Which results in output similar to the following figure (link to figure)

.. figure:: example/Ecoli_first352k.ale.pdf.png
   :align:  center

For a full list of options please see the :mod:`plotter3` or run::

  $ ./plotter3.py -h
  Usage: ./ALE_plotter.py [-options] <inputfile.ale>

  where basic options are:
    -h      : show brief help on version and full usage
    -nosave : do not save the figure as a pdf (instead plot to screen)

  parameter options accepting <f>loats and <i>ntegers and <s>trings (default):
    -s <i>   : the starting position to plot (for all contigs, ie a single insert length)
    -e <i>   : the ending position of the plot
    -pt <s>  : plot type 't'otal 'k'mer 'p'lacement 'd'epth (-pt dpkt)
    -dsw <i> : depth smoothing window, averaging over position (-dsw 10000)
    -psw <i> : placement smoothing window (-psw 1000)
    -ksw <i> : kmer smoothing window (-ksw 1000)
    -t <f>   : threshold percentage, see paper (-t 0.99999)
    -pt <f>  : plot threshold, only plot if more than % of errors (-pt 0.0)
    -st <i>  : number of standard deviations to engage threshold (-st 5)
    -fn <s>  : figure name (default: contig name)
    -mps <i> : minimum plot size in bp (-mps 20000)
    -sc <s>  : plot only a specific contig (ie -sc contigName213)
    -pmo     : plot meta information only (off)
    -dpm     : don't plot meta information at all (off)

