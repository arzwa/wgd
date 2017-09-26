.. wgd documentation master file, created by
   sphinx-quickstart on Mon Apr 10 10:31:22 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

wgd: Whole genome duplication analysis in Python
************************************************

This Python package and corresponding command line interface (CLI) were developed for various analyses related
to whole genome duplications (WGDs). Here the Python API Package is documented as well as the various command line
utilities bundled in the ``wgd`` CLI.

To install ``wgd``, clone the repository available at https://github.com/arzwa/wgd, navigate into the repo, run
``pip install .``.

Make sure you have the following external software tools installed:

1. PAML (http://abacus.gene.ucl.ac.uk/software/paml.html)
2. MUSCLE (http://www.drive5.com/muscle/)
3. Optional: ``mcl`` (https://micans.org/mcl/index.html)
4. Optional: Blastp (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)


Python package
==============

Contents:

.. toctree::
   :maxdepth: 2

   Markov clustering (mcl) <mcl>
   Codeml wrapper <codeml>
   Alignment tools <alignment>
   Ks distribution analysis <ks>
   Positive selection screening <pos>

Command line tools
==================

Note that besides using the python library directly, you can also use the provided command line utilities.
Upon cloning/downloading, run `wgd` at the command to see the various available commands.
These are written in python using the ``click`` library and require the ``wgd`` module and its dependencies to be
installed correctly.

* :ref:`cli`

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

