[![Documentation Status](https://readthedocs.org/projects/wgd/badge/?version=latest)](http://wgd.readthedocs.io/en/latest/?badge=latest)

Copyright (C) 2018 Arthur Zwaenepoel

# TO DO

- It is not necessary to strip gaps, as gaps are treated as missing data!

- So I think the best approach is just run codeml on a full codon alignment,
save the pairwise alignment length for ech pair, filter out those that are 
too short and proceed as normal. This will be equivalent as the pairwise 
method but way faster!

- We could also use codeml with tree! So first infer a tree with codonphyml,
then run codeml on the full alignment + tree, and finally collect ds values 
and weights. This seems to me the most consistent approach. As pairwise 
estimates are actually inconsistent when sequences are related in a tree 
structure!

- Am I using the best model for Ks estimation? Gamma rates? If using full 
alignment, do I need to use F3x4?

# Whole genome duplication analysis

Python package and command line interface (CLI) for the analysis
of whole genome duplications (WGDs). Tested with Python3.5 & Python3.6
on Linux Ubuntu. If you don't have python or pip installed a simple
`sudo apt-get install python3-pip` will do on Ubuntu, Debian or Linux Mint.

To install: clone the repo, navigate to it and install it with pip

    $ git clone https://github.com/arzwa/wgd.git
    $ cd wgd
    $ pip install .

For the command line interface, upon installation run

    $ wgd

to get a list of the available commands. To get usage instructions for
a command (e.g. `ks`) run

    $ wgd ks --help

To use as a Python package as well as to find additional documentation
for the CLI, please consult the docs at http://wgd.readthedocs.io/en/latest/

## Third party software

`wgd` requires the following third party executables (preferably these should also be in the `PATH` environment variable):

For `wgd blast`:
- `blast` (Altschul _et al._ 2008), from which it uses the `blastp` and `makeblastdb` commands,
`sudo apt-get install ncbi-blast+` will often suffice for installation
- `mcl` (https://micans.org/mcl/index.html). Get MCL using your package
manager ``sudo apt-get install mcl`` or download it at the provided link.

For `wgd ks`:
- `muscle` (Edgar 2004) for multiple sequence alignment. MUSCLE can be
retrieved from http://www.drive5.com/muscle/.
- `codeml` from the PAML software package (Yang 1997). PAML can be downloaded
from the following link: http://abacus.gene.ucl.ac.uk/software/paml.html
- Depending on which weighting method you choose for node-weighting the Ks distribution, you might also
need `FastTree` (http://www.microbesonline.org/fasttree/) (Price _et al._ 2010) or ``phyml``
(http://www.atgc-montpellier.fr/phyml/) (Guindon & Gascuel 2003) (not necessary however).

For `wgd syn`
- ``i-adhore`` from the I-ADHoRe 3.0 suite (http://bioinformatics.psb.ugent.be/beg/tools/i-adhore30) (Proost _et al._ 2012)

Of course, you don't need all of these tools installed for each feature of ``wgd``. Note that **`wgd` is still under development**, if you use `wgd` in your work, please cite this repository (when the first release is ready, this will have a DOI).
