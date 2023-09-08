.. ALE documentation master file, created by
   sphinx-quickstart on Fri Dec 16 21:11:32 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

synthReadGen
============

Jump to:
   #. `Compiling synthReadGen`_
   #. `Running synthReadGen`_

Compiling synthReadGen
----------------------

By running **make** in the src/ directory *synthReadGen* should be automatically compiled::
   
   $ cd src
   $ make

Alternately, you can compile it manually with::

   $ cc -g -O3 synthReadGen.c -o synthReadGen -lz -lm -Isamtools-0.1.18 -Lsamtools-0.1.18

Running synthReadGen
--------------------

The usage can be found by running::

   $./synthReadGen -h
   Usage: ./synthReadGen [options] <inputFile> <outputFile>
   
   Options: <i>nt <f>loat [default]
     -h      : print out this help
     -id <i> : set distribution used for insert length
               [1 = normal], 2 = poisson
     -ld <i> : set distribution used for read length
               [1 = normal], 2 = poisson
     -im <f> : inward insert length mean [200.0]
     -om <f> : outward insert length mean [500.0]
     -is <f> : inward insert length std dev [10.0]
     -os <f> : outward insert length std dev [15.0]
     -ip <f> : probability for an inward read [0.5]
     -er <c> : illumina error char [^]
     -nr <i> : number of reads to make [1000]
     -rl <x> : read length mean [85.0]
     -rs <x> : read length sigma [7.0]
     -ps <x> : no error for first x bases in a read [0]
     -b      : outputs two fastq files for bowtie mapping [off]

