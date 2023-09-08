.. ALE documentation master file, created by
   sphinx-quickstart on Fri Dec 16 21:11:32 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

plotter3.py
===========

Jump to:
   #. `Requirements for plotter3.py`_
   #. `Running plotter3.py`_
   #. `plotter3.py functions and classes`_

Plotting
--------

The authors recommend using IGV to view the output.

http://www.broadinstitute.org/igv/

Just import the assembly, bam and ALE scores. You can convert the .ale file to a set of .wig files with ale2wiggle.py and IGV can read those directly.  Depending on your genome size you may want to convert the .wig files to the BigWig format.

Requirements for plotter4.py
----------------------------

.. automodule:: plotter3

Running plotter3.py
-------------------

Here we show how to use the built in debugging plotter. Note, this plotter is no longer in active development and is not supported. We recommend using IGV.

We can invoke the plotter by running::

   $ ./plotter3.py ALEoutput.ale

Which results in output similar to the following figure (link to figure)

.. figure:: example/Ecoli_first350k.ale.pdf.png
   :align:  center

For a full list of options please see the documentation/source below or run::

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

plotter3.py functions and classes
---------------------------------

.. autofunction:: plotter3.read_in_info

.. autoclass:: plotter3.Contig
   :members:

