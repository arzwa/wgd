[![Documentation Status](https://readthedocs.org/projects/wgd/badge/?version=latest)](http://wgd.readthedocs.io/en/latest/?badge=latest)

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

For the command line interface, upon installation run

    $ wgd

to get a list of the available commands. To get usage instructions for
a command (e.g. `ksd`) run

    $ wgd ksd --help

To use as a Python package as well as to find additional documentation
and examples for the CLI, please consult the docs at
http://wgd.readthedocs.io/en/latest/

## Notes

**Note on input data:** while the input data is rather straightforward
(a CDS fasta file will do for most analyses) it may be of interest that
the wgd suite was extensively tested with data from the PLAZA platform,
so for examples of the right input data formats (in particular CDS fasta
files for sequence data and GFF files for structural annotation), please
have a look [there](https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v4_dicots/download/).

**Note on virtualenv:** you can install wgd in a _virtual environment_
(using `virtualenv`). If you would however encounter problems with
running the executable directly (e.g. `wgd --help` doesn't work) you can
circumvent this by directly calling the CLI, using `python3 ./wgd_cli.py
--help` (assuming you are currently in the directory where you cloned
wgd).

## Citation

Until this package is described in a formal publication, please cite
this repository if you use the software in your research. For citation
of the tools used in wgd, please consult the documentation at
http://wgd.readthedocs.io/en/latest/.