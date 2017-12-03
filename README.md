# Whole genome duplication analysis

Python package and command line interface (CLI) for the analysis
of whole genome duplications (WGDs). Tested with Python3.5 & Python3.6
on Linux Ubuntu. If you don't have python or pip installed a simple
`sudo apt-get install python3-pip` will do on Ubuntu, Debian or Linux Mint.

To install: clone the repo, navigate to it and install it with pip

    $ git clone https://gitlab.psb.ugent.be/arzwa/wgd.git
    $ cd wgd
    $ pip install .

For the command line interface, upon installation run

    $ wgd

to get a list of the available commands. To get usage instructions for
a command (e.g. `ks`) run

    $ wgd ks --help

To use as a Python package as well as to find additional documentation
for the CLI, please consult the **extensive documentation** in the `doc`
directory or at https://arzwa.github.io/wgd/index.html.

## Third party software

`wgd` requires the following third party executables (preferably these
should also be in the `PATH` environment variable):

- `blast` (Altschul _et al._ 2008), from which it uses the `blastp` and `makeblastdb` commands,
`sudo apt-get install ncbi-blast+` will often suffice for installation
- ``mcl`` (https://micans.org/mcl/index.html). Get MCL using your package
manager ``sudo apt-get install mcl`` or download it at the provided link.
- `muscle` (Edgar 2004) for multiple sequence alignment. MUSCLE can be
retrieved from http://www.drive5.com/muscle/.
- `codeml` from the PAML software package (Yang 1997). PAML can be downloaded
from the following link: http://abacus.gene.ucl.ac.uk/software/paml.html
- `i-adhore` (Proost _et al._ 2012) for co-linearity analysis. I-ADHoRe 3.0
can be obtained from http://bioinformatics.psb.ugent.be/beg/tools/i-adhore30






