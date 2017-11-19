.. wgd documentation master file, created by
   sphinx-quickstart on Mon Apr 10 10:31:22 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

wgd: Whole genome duplication analysis in Python
************************************************

This Python package and corresponding command line interface (CLI) were developed for various analyses related
to whole genome duplications (WGDs). Here the Python API is documented as well as the various command line
utilities bundled in the ``wgd`` CLI.

To install ``wgd``, clone the repository available at https://github.com/arzwa/wgd, navigate into the repo, run
``pip install .``.

External software
=================

For full functionality, make sure you have the following external software tools installed:

1. ``codeml`` from the PAML package (http://abacus.gene.ucl.ac.uk/software/paml.html)
2. ``muscle`` for multiple sequence alignment (http://www.drive5.com/muscle/)
3. ``mcl`` (https://micans.org/mcl/index.html)
4. ``blastp`` and ``makeblastdb`` from the Blast+ suite (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
5. ``i-adhore`` from the I-ADHoRe 3.0 suite (http://bioinformatics.psb.ugent.be/beg/tools/i-adhore30)

Of course, you don't need all of these tools installed for each feature of ``wgd``.
Please put these tools in your PATH environment variable as well, it makes life easier.

Command line tools
==================

The command line tools are the main functionality of the ``wgd`` package. You can find tools for the following
analyses:

1. All-*versus*-all Blastp analysis and MCL clustering
2. Whole paranome |Ks| (and |Ka| and ω) distribution construction
3. One-versus-one orthologs |Ks| (and |Ka| and ω;) distribution construction
4. Mixture modeling of |Ks| distributions and WGD-specific paralog extraction
5. Visualization of (multiple) |Ks| distributions
6. Intragenomic co-linearity/synteny analysis and anchor based |Ks| distribution construction
7. Co-linearity dotplot construction

All relative information can be find here:

.. toctree::
   :maxdepth: 1

   Command line interface <cli>


Python package
==============

For those interested in the underlying structure of ``wgd``, here you can find the full documentation of the API.

Contents:

.. toctree::
   :maxdepth: 2

   Blast & Markov clustering (MCL) <blast_mcl>
   Codeml anlysis <codeml>
   Alignment tools <alignment>
   Ks distribution analysis <ks>
   Mixture modeling <mix>
   Co-linearity analysis <syn>
   Utilities <utils>
   Visualization <viz>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. |Ks| replace:: K\ :sub:`S`
.. |Ka| replace:: K\ :sub:`A`