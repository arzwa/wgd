[![Documentation Status](https://readthedocs.org/projects/wgd/badge/?version=latest)](http://wgd.readthedocs.io/en/latest/?badge=latest) [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2097)


Copyright (C) 2018 Arthur Zwaenepoel

VIB/UGent center for plant systems biology -
Bioinformatics & evolutionary genomics group http://bioinformatics.psb.ugent.be/beg/

# wgd - simple command line tools for the analysis of ancient whole genome duplications

## Installation

Python package and command line interface (CLI) for the analysis
of whole genome duplications (WGDs). Tested with Python3.5 & Python3.6
on Linux Ubuntu. If you don't have python or pip installed a simple
`sudo apt-get install python3-pip` should do.

To install: clone the repo, navigate to it and install it with pip

    $ git clone https://github.com/arzwa/wgd.git
    $ cd wgd
    $ pip install .

Note that depending on your python installation and whether you're in a
virtualenv, ``pip`` may default either to ``pip2`` or ``pip3``. If the
above installation step fails, please try to use ``pip3`` instead of
``pip``.

For the command line interface, upon installation run

    $ wgd

to get a list of the available commands. To get usage instructions for
a command (e.g. `ksd`) run

    $ wgd ksd --help

For **external software** requirements: please consult the relevant section
in the docs: https://wgd.readthedocs.io/en/latest/index.html#external-software

To use as a Python package as well as to find additional documentation
and examples for the CLI, please consult the docs at
http://wgd.readthedocs.io/en/latest/

## Singularity container

A singularity container is available for ``wgd``, allowing all to use
all tools in ``wgd`` except ``wgd syn``, without having to install all
required software on your system. To install Singularity follow
the instructions [here](https://www.sylabs.io/docs/)

If you have singulaity installed (and you're in the virtual machine when
running on Windows or Mac), you can run the following to get the container

    singularity pull --name wgd.simg shub://arzwa/wgd

Then you can use ``wgd`` as follows

    singularity exec wgd.simg wgd <command>

## Changes

- 08/01/2019: fixed ImportError for interactive histogram visualization

## Notes

**Bug tracking:** If the program crashes, exits unexpectedly or some
unexpected results are obtained, please run it again with the
``--verbosity debug`` flag *before* the subcommand of interest (*e.g.*
``wgd --verbosity debug ksd gf.mcl cds.fasta``). If the anomaly persists,
please open an issue on this GitHub site.

**Note on input data:** while the input data is rather straightforward
(a CDS fasta file will do for most analyses) it may be of interest that
the wgd suite was extensively tested with data from the PLAZA platform,
so for examples of the right input data formats (in particular CDS fasta
files for sequence data and GFF files for structural annotation), please
have a look [there](https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v4_dicots/download/).
It is generally advised not to include pipe characters (`|`) in your gene 
IDs, since these can have special meanings in certain parts of `wgd`.

**Note on virtualenv:** you can install wgd in a _virtual environment_
(using [`virtualenv`](https://virtualenv.pypa.io/en/stable/)). If you
would however encounter problems with running the executable directly
(e.g. `wgd --help` doesn't work) you can circumvent this by directly
calling the CLI, using `python3 ./wgd_cli.py --help` (assuming you are
currently in the directory where you cloned wgd).

## Citation
 
Please cite us at https://doi.org/10.1093/bioinformatics/bty915

```
Zwaenepoel, A., and Van de Peer, Y. wgd - simple command line tools for the analysis of ancient whole genome duplications. Bioinformatics., bty915, https://doi.org/10.1093/bioinformatics/bty915
```

For citation of the tools used in wgd, please consult the documentation at
https://wgd.readthedocs.io/en/latest/index.html#citation.

