.. ALE documentation master file, created by
   sphinx-quickstart on Fri Dec 16 21:11:32 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

artificial_errors.py
====================

Jump to:
   #. `Running artificial_errors.py`_
   #. `artificial_errors.py functions and classes`_

Running artificial_errors.py
----------------------------

We can invoke the artificial error maker by running::

   $ ./artificial_errors.py [-options] <inputfile.fna>

This will create a new file errors_<inputfile.fna> that has the transformations requested in [-options] (performed left to right). If no options are given the errors_<inputfile.fna> will be identical to <inputfile.fna>.

Options are::

    $ ./artificial_errors.py -h
    Usage: ./artificial_errors.py [-options] <inputfile.fna>

    where basic options are:
      -h      : show brief help on version and full usage

    parameter options accepting <i>ntegers and <s>trings (default):
      Note: transformations will be made left to right
      -ase <i> <i> : add substitution error at <location> for <length>
      -ade <i> <i> : add deletion error at <location> for <length>
      -aie <i> <i> : add insertion error at <location> for <length>
      -inv <i> <i> : add inversion error at <location> for <length>
      -cip <i> <i> : copy part of the assembly at <location> for <length>
      -trp <i>     : transpose assembly around <pivot>
      -ab  <i>     : add a break (split into 2 contigs) at <location>
      -o   <s>     : output file name (error_ + inputfile.fna)

artificial_errors.py functions and classes
------------------------------------------

.. automodule:: artificial_errors
   :members:



