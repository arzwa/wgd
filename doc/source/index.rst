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

``wgd`` requires the following third party executables (preferably these should also be in the `PATH` environment variable):

For ``wgd blast``:
- ``blast`` (Altschul _et al._ 2008), from which it uses the ``blastp`` and ``makeblastdb`` commands,
``sudo apt-get install ncbi-blast+`` will often suffice for installation
- ``mcl`` (https://micans.org/mcl/index.html). Get MCL using your package
manager ``sudo apt-get install mcl`` or download it at the provided link.

For ``wgd ks``:
- ``muscle`` (Edgar 2004) for multiple sequence alignment. MUSCLE can be
retrieved from http://www.drive5.com/muscle/.
- ``codeml`` from the PAML software package (Yang 1997). PAML can be downloaded
from the following link: http://abacus.gene.ucl.ac.uk/software/paml.html
- Depending on which weighting method you choose for node-weighting the Ks distribution, you might also
need ``FastTree`` (http://www.microbesonline.org/fasttree/) (Price `et al.` 2010) or ``phyml``
(http://www.atgc-montpellier.fr/phyml/) (Guindon & Gascuel 2003) (not necessary however).

For ``wgd syn``
- ``i-adhore`` from the I-ADHoRe 3.0 suite (http://bioinformatics.psb.ugent.be/beg/tools/i-adhore30) (Proost `et al.` 2012)

Of course, you don't need all of these tools installed for each feature of ``wgd``.

Command line tools
==================

The command line tools are the main functionality of the ``wgd`` package. You can find tools for the following
analyses:

1. All-*versus*-all Blastp analysis and MCL clustering
2. Whole paranome |Ks| (and |Ka| and ω) distribution construction
3. One-versus-one orthologs |Ks| (and |Ka| and ω;) distribution construction
4. Mixture modeling of |Ks| distributions and WGD-specific paralog extraction
5. Interactive visualization of (multiple) |Ks| distributions and kernel density estimates thereof
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
   Phylogenetics tools <phy>
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